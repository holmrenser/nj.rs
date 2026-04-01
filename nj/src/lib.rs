//! Neighbor-Joining phylogenetic tree inference library.
//!
//! # Data flow
//!
//! ```text
//! [FASTA / Python dict / JS object]
//!         │
//!         ▼
//!      NJConfig  (config.rs)
//!         │
//!         ▼
//!   detect_alphabet()  ──►  Alphabet::DNA | Alphabet::Protein
//!         │
//!         ▼
//!     MSA<DNA|Protein>  (msa.rs)
//!      ├── bootstrap() ──► bootstrap_clade_counts()
//!      └── into_dist::<Model>()
//!               │
//!               ▼
//!           DistMat  (dist.rs)
//!               │
//!               ▼
//!         neighbor_joining()  ──►  NJState::run()  (nj.rs)
//!               │
//!               ▼
//!           TreeNode  (tree.rs)
//!               │
//!               ▼
//!           to_newick()  ──►  Newick String
//! ```
//!
//! # Public API
//!
//! The single public entry point is [`nj`], which accepts an [`NJConfig`] and
//! returns a Newick string. Everything else is internal implementation detail
//! exposed only to the Python and WASM wrapper crates.
//!
//! # Model–alphabet compatibility
//!
//! | Model | DNA | Protein |
//! |-------|-----|---------|
//! | `PDiff` | ✓ | ✓ |
//! | `JukesCantor` | ✓ | — |
//! | `Kimura2P` | ✓ | — |
//! | `Poisson` | — | ✓ |
//!
//! Providing an incompatible model returns an `Err` from [`nj`].
pub mod alphabet;
pub mod config;
pub mod distance_matrix;
pub mod models;
pub mod msa;
pub mod nj;
pub mod tree;

use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::HashMap;

use crate::alphabet::{Alphabet, AlphabetEncoding, DNA, Protein};
use crate::config::SubstitutionModel;
pub use crate::config::{MSA, NJConfig, SequenceObject};
use crate::distance_matrix::DistMat;
use crate::models::{JukesCantor, Kimura2P, ModelCalculation, PDiff, Poisson};
use crate::tree::{NameOrSupport, TreeNode};

/// Fills `out` with the leaf indices of all taxa in the subtree rooted at `node`.
///
/// Bits in `out` are set to `true` for each leaf encountered. The bit position
/// is looked up in `idx` by the leaf's name label. Returns `Err` if a leaf
/// has no name label (should not occur for well-formed NJ trees).
fn bitset_of(
    node: &TreeNode,
    idx: &HashMap<String, usize>,
    out: &mut BitVec<u8, Lsb0>,
) -> Result<(), String> {
    match &node.children {
        None => match &node.label {
            Some(NameOrSupport::Name(name)) => {
                let i = idx[name];
                out.set(i, true);
                Ok(())
            }
            _ => Err("Leaf node without a name label".into()),
        },
        Some([l, r]) => {
            bitset_of(l, idx, out)?;
            bitset_of(r, idx, out)?;
            Ok(())
        }
    }
}

/// Recursively counts how many times each non-trivial clade appears in `tree`.
///
/// A clade is represented as a raw-byte encoding of a `BitVec` over the `n_taxa`
/// leaf indices. Only clades with `1 < size < n_taxa` (i.e. proper internal
/// clades) are counted. Each call increments the clade's entry in `counter` by 1.
/// Used by [`bootstrap_clade_counts`] to aggregate over bootstrap replicates.
fn count_clades(
    tree: &TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counter: &mut HashMap<Vec<u8>, usize>,
) -> Result<(), String> {
    if let Some([l, r]) = &tree.children {
        // compute clade bitvec
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(tree, idx, &mut bv).unwrap();

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            // unique by structure; no HashSet needed
            counter
                .entry(bv.as_raw_slice().to_vec())
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }

        // recursion
        count_clades(l, idx, n_taxa, counter)?;
        count_clades(r, idx, n_taxa, counter)?;
    }
    Ok(())
}

/// Performs bootstrap sampling and counts clades across bootstrap trees.
fn bootstrap_clade_counts<A: AlphabetEncoding, M: ModelCalculation<A>>(
    msa: &MSA<A>,
    n_bootstrap_samples: usize,
) -> Result<Option<HashMap<Vec<u8>, usize>>, String> {
    if n_bootstrap_samples == 0 {
        return Ok(None);
    }
    let idx_map: HashMap<String, usize> = msa.to_index_map();
    let mut counter = HashMap::new();
    for _ in 0..n_bootstrap_samples {
        let tree = msa
            .bootstrap()?
            .into_dist::<M>()
            .neighbor_joining()
            .expect("NJ bootstrap iteration failed");
        count_clades(&tree, &idx_map, msa.n_sequences, &mut counter)?;
    }
    Ok(Some(counter))
}

/// Annotates internal nodes with bootstrap support values from `counts`.
///
/// For each internal node, computes its clade `BitVec`, looks up the count in
/// `counts`, and assigns a [`NameOrSupport::Support`] label if a matching
/// entry is found. Nodes whose clade was never observed in bootstrap replicates
/// receive no label.
fn add_bootstrap_to_tree(
    node: &mut TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counts: &HashMap<Vec<u8>, usize>,
) {
    if node.children.is_some() {
        // compute clade bitvec
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(node, idx, &mut bv).unwrap();

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            if let Some(c) = counts.get(&bv.as_raw_slice().to_vec()) {
                node.label = Some(NameOrSupport::Support(*c));
            }
        }

        // recursion
        if let Some([l, r]) = &mut node.children {
            add_bootstrap_to_tree(l, idx, n_taxa, counts);
            add_bootstrap_to_tree(r, idx, n_taxa, counts);
        }
    }
}

/// Heuristically detects whether the MSA contains DNA or protein sequences.
///
/// Returns [`Alphabet::DNA`] unless any sequence contains a byte that is not
/// in `{A, C, G, T, U, N, -}` (case-insensitive), in which case
/// [`Alphabet::Protein`] is returned. This covers all 20 standard amino acids
/// since letters like `D`, `E`, `F`, `H`, `I`, `K`, `L`, `M`, `P`, `Q`,
/// `R`, `S`, `V`, `W`, `Y` cannot appear in a DNA alignment.
fn detect_alphabet(msa: &[SequenceObject]) -> Result<Alphabet, String> {
    // Simple heuristic: if any character is > A,C,G,T,N, assume protein
    let mut is_protein = false;

    for seq in msa {
        for c in seq.sequence.bytes() {
            match c.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'-' => { /* still possible DNA */ }
                _ => {
                    is_protein = true;
                    break;
                }
            }
        }
        if is_protein {
            break;
        }
    }

    Ok(if is_protein {
        Alphabet::Protein
    } else {
        Alphabet::DNA
    })
}

/// Runs NJ with model `M` on alphabet `A` and returns a Newick string.
///
/// If `n_bootstrap_samples > 0`, generates that many bootstrap replicates,
/// collects clade counts via [`bootstrap_clade_counts`], runs NJ on the
/// original distances, and annotates the tree before serialising to Newick.
fn run_nj<A, M>(msa: MSA<A>, n_bootstrap_samples: usize) -> Result<String, String>
where
    A: AlphabetEncoding,
    M: ModelCalculation<A>,
{
    // bootstrap_clade_counts should be generic over A,M too (not shown here)
    let clade_counts = bootstrap_clade_counts::<A, M>(&msa, n_bootstrap_samples)?;

    let mut main_tree = msa.into_dist::<M>().neighbor_joining()?;
    let newick = match clade_counts {
        Some(counts) => {
            let main_idx_map: HashMap<String, usize> = msa.to_index_map();
            add_bootstrap_to_tree(&mut main_tree, &main_idx_map, msa.n_sequences, &counts);
            main_tree.to_newick()
        }
        None => main_tree.to_newick(),
    };
    Ok(newick)
}

/// Infers a phylogenetic tree from an aligned MSA and returns a Newick string.
///
/// This is the single public entry point for the library. The alphabet is
/// auto-detected from the sequences; `conf.substitution_model` must be
/// compatible with the detected alphabet (see the module-level compatibility
/// table). Returns `Err` for an empty MSA, an incompatible model, or any
/// internal NJ failure.
pub fn nj(conf: NJConfig) -> Result<String, String> {
    if conf.msa.is_empty() {
        return Err("Input MSA is empty".into());
    }
    let alphabet = detect_alphabet(&conf.msa)?;
    match alphabet {
        Alphabet::DNA => {
            // build MSA specialized to DNA (pre-encodes using DNA::encode)
            let msa =
                MSA::<DNA>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));

            match conf.substitution_model {
                SubstitutionModel::PDiff => run_nj::<DNA, PDiff>(msa, conf.n_bootstrap_samples),
                SubstitutionModel::JukesCantor => {
                    run_nj::<DNA, JukesCantor>(msa, conf.n_bootstrap_samples)
                }
                SubstitutionModel::Kimura2P => {
                    run_nj::<DNA, Kimura2P>(msa, conf.n_bootstrap_samples)
                }
                // Poisson is a protein model — either disallow here or handle by error:
                SubstitutionModel::Poisson => {
                    Err("Poisson is a protein model; cannot use with DNA".into())
                }
            }
        }

        Alphabet::Protein => {
            let msa =
                MSA::<Protein>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));

            match conf.substitution_model {
                // Poisson is valid for protein:
                SubstitutionModel::Poisson => {
                    run_nj::<Protein, Poisson>(msa, conf.n_bootstrap_samples)
                }
                SubstitutionModel::PDiff => run_nj::<Protein, PDiff>(msa, conf.n_bootstrap_samples),
                // DNA-only models should be rejected for proteins:
                SubstitutionModel::JukesCantor | SubstitutionModel::Kimura2P => {
                    Err("Selected model is for DNA; cannot use with Protein".into())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::SubstitutionModel;

    #[test]
    fn test_nj_wrapper_simple_tree() {
        let sequences = vec![
            SequenceObject {
                identifier: "A".into(),
                sequence: "ACGTCG".into(),
            },
            SequenceObject {
                identifier: "B".into(),
                sequence: "ACG-GC".into(),
            },
        ];
        let conf = NJConfig {
            msa: sequences,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
        };
        let newick = nj(conf).expect("NJ failed");
        assert_eq!(newick, "(A:0.167,B:0.167);");
    }

    #[test]
    fn test_nj_wrapper_adds_semicolon() {
        let sequences = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "A".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "A".into(),
            },
        ];
        let conf = NJConfig {
            msa: sequences,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
        };
        let out = nj(conf).unwrap();
        assert!(out.ends_with(';'));
    }

    #[test]
    fn test_nj_deterministic_order() {
        let sequences = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "ACGTCG".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "ACG-GC".into(),
            },
            SequenceObject {
                identifier: "Seq2".into(),
                sequence: "ACGCGT".into(),
            },
        ];
        let conf = NJConfig {
            msa: sequences,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
        };

        let t1 = nj(conf.clone()).unwrap();
        let t2 = nj(conf).unwrap();
        assert_eq!(t1, t2);
    }

    #[test]
    fn test_nj_wrapper_empty_msa() {
        let conf = NJConfig {
            msa: vec![],
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
        };
        let result = nj(conf);
        assert!(result.is_err());
    }

    #[test]
    fn test_nj_wrapper_incorrect_model_for_alphabet() {
        let sequences = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "ACGTCG".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "ACG-GC".into(),
            },
        ];
        let conf = NJConfig {
            msa: sequences,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::Poisson, // protein model for DNA MSA
        };
        let result = nj(conf);
        assert!(result.is_err());
    }

    #[test]
    fn test_nj_wrapper_incorrect_model_for_protein() {
        let sequences = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "ACDEFGH".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "ACD-FGH".into(),
            },
        ];
        let conf = NJConfig {
            msa: sequences,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::JukesCantor, // DNA model for protein MSA
        };
        let result = nj(conf);
        assert!(result.is_err());
    }

    #[test]
    fn test_detect_alphabet_dna() {
        let msa = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "ACGTACGT".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "ACG-ACGT".into(),
            },
        ];
        let alphabet = detect_alphabet(&msa).expect("detection failed");
        assert_eq!(alphabet, Alphabet::DNA);
    }

    #[test]
    fn test_detect_alphabet_protein() {
        let msa = vec![
            SequenceObject {
                identifier: "Seq0".into(),
                sequence: "ACDEFGHIK".into(),
            },
            SequenceObject {
                identifier: "Seq1".into(),
                sequence: "ACD-FGHIK".into(),
            },
        ];
        let alphabet = detect_alphabet(&msa).expect("detection failed");
        assert_eq!(alphabet, Alphabet::Protein);
    }
}
