//! Neighbor-Joining tree builder and triangular distance matrix utilities.
//!
//! This module provides a compact representation for pairwise distances (TriMat),
//! a simple tree node structure used to represent Neighbor-Joining (NJ) trees
//! (NJNode), a robust, borrow-safe implementation of the Neighbor-Joining
//! algorithm (neighbor_joining), and a helper to build a TriMat from a multiple
//! sequence alignment (tri_from_msa).
//!
//! # Types
//!
//! - TriMat
//!   - Compact triangular (upper-triangular stored as a flat Vec) distance
//!     matrix with associated node names.
//!   - Access via `get(i, j)` and `set(i, j, value)`. Diagonal distances are
//!     always treated as 0. The internal index mapping follows a strict
//!     triangular packing so only i != j is stored.
//!   - Construct with `TriMat::with_names(names)`; `dim()` returns the number
//!     of taxa.
//!
//! - TreeNode
//!   - Binary tree node representing either a leaf (no children) or an
//!     internal node with exactly two children.
//!   - Stores child branch lengths on the parent as `left_len` and `right_len`.
//!   - `to_newick(hide_internal)` emits a Newick string with branch lengths
//!     formatted to 3 decimal places. When `hide_internal` is true the internal
//!     node labels are omitted.
//!
//! # Functions
//!
//! - neighbor_joining(dist: TriMat) -> Result<NJNode, String>
//!   - Performs the Neighbor-Joining algorithm on the provided triangular
//!     distance matrix and returns the inferred tree as an `NJNode`.
//!   - The implementation is iterator-oriented and avoids borrowing issues by
//!     moving and reorganizing node ownership via a HashMap of indices -> nodes.
//!   - Uses a BitVec to track active taxa and maintains a `row_sums` vector
//!     for efficient Q-matrix computation.
//!   - Branch lengths computed for joins are clamped to >= 0.0 to avoid
//!     negative lengths caused by rounding or non-ultrametric inputs.
//!   - Returns an `Err` if the input matrix is empty or if an internal
//!     consistency check fails (e.g., no pair to join or unexpected number of
//!     remaining active nodes).
//!   - Complexity: O(n^2) memory for the distance representation and O(n^3)
//!     worst-case time if implemented naively; this implementation reduces
//!     repeated work via row sums but still iterates over active pairs while
//!     choosing minimal Q entries.
//!
//! - tri_from_msa(msa: &[String], names: Option<&[String]>) -> TriMat
//!   - Builds a `TriMat` from a multiple sequence alignment (MSA).
//!   - Pairwise distance = (number of differing positions) / (number of
//!     compared positions), ignoring any column where either sequence has a
//!     gap (`'-'`). If two sequences share no compared positions the distance
//!     is set to 0.0.
//!   - If `names` is omitted, default names are generated as `Seq0`, `Seq1`, ....
//!
//! # Usage notes and invariants
//!
//! - TriMat stores only i != j entries; `get(i,i)` is defined to be 0.0.
//! - Indices used with TriMat must be valid (0..dim()). `set` silently ignores
//!   attempts to set diagonal entries.
//! - neighbor_joining consumes the distance matrix (it is taken by value) and
//!   returns a tree whose leaves preserve the original node name strings.
//! - The NJ implementation expects a symmetric distance matrix; asymmetries
//!   will lead to undefined/incorrect results.
//! - Branch lengths are clamped to non-negative values to keep the tree valid
//!   in downstream tools; if you require different handling of negative lengths
//!   adjust the implementation accordingly.
//!
//! # Examples
//!
//! Construct a TriMat from an MSA and produce a Newick tree:
//!
//! ```ignore
//! let msa = vec!["ACG".into(), "ATG".into(), "A-G".into()];
//! let tm = tri_from_msa(&msa, None);
//! let tree = neighbor_joining(tm).expect("NJ failed");
//! let newick = tree.to_newick();
//! ```
//!
//! # Errors & Panics
//!
//! - Functions in this module return `Result` where appropriate to indicate
//!   recoverable errors (e.g., empty input). Internal indexing errors return
//!   descriptive strings rather than panicking.
//! - The code assumes reasonably small `n` such that triangular storage in a
//!   Vec is practical; very large numbers of taxa may cause memory pressure.
//!
//! # Testing
//!
//! The accompanying unit tests exercise TriMat indexing, distance-building
//! behavior with gaps, two- and three-taxon NJ behavior, and Newick output
//! formatting. They also validate that branch lengths are non-negative after
//! clamping.
mod alphabet;
mod config;
mod dist;
pub mod models;
mod msa;
mod nj;
mod tree;
use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::HashMap;

use crate::alphabet::{Alphabet, AlphabetEncoding, DNA, Protein};
use crate::config::SubstitutionModel;
pub use crate::config::{MSA, NJConfig, SequenceObject};
use crate::dist::DistMat;
use crate::models::{JukesCantor, Kimura2P, ModelCalculation, PDiff, Poisson};
use crate::tree::{NameOrSupport, TreeNode};

/// Collect leaves into a BitVec according to index map.
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

pub fn count_clades(
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
pub fn bootstrap_clade_counts<A: AlphabetEncoding, M: ModelCalculation<A>>(
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

/// Recursively adds bootstrap support values to internal nodes of the tree.
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

pub fn detect_alphabet(msa: &[SequenceObject]) -> Result<Alphabet, String> {
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

/// Runs Neighbor-Joining on the given MSA with specified model and number of bootstrap samples.
/// Returns Newick string of the resulting tree.
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

/// Main Neighbor-Joining wrapper function.
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
