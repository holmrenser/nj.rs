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
pub mod error;
pub mod event;
pub mod fasta;
pub mod models;
pub mod msa;
pub mod nj;
pub mod tree;

use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::HashMap;

use crate::alphabet::{Alphabet, AlphabetEncoding, DNA, Protein};
use crate::config::SubstitutionModel;
pub use crate::config::{DistConfig, MSA, NJConfig, SequenceObject};
use crate::distance_matrix::DistMat;
pub use crate::distance_matrix::DistanceResult;
pub use crate::error::NJError;
pub use crate::event::{LogLevel, NJEvent};
pub use crate::fasta::parse_fasta;
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
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(tree, idx, &mut bv)?;

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            counter
                .entry(bv.as_raw_slice().to_vec())
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }

        count_clades(l, idx, n_taxa, counter)?;
        count_clades(r, idx, n_taxa, counter)?;
    }
    Ok(())
}

/// Builds a Rayon thread pool with `num_threads` workers.
///
/// When `num_threads` is `None`, Rayon uses its default (one thread per logical CPU).
/// Returns `Err` if the pool cannot be constructed.
#[cfg(feature = "parallel")]
pub(crate) fn build_thread_pool(num_threads: Option<usize>) -> Result<rayon::ThreadPool, String> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    if let Some(n) = num_threads {
        builder = builder.num_threads(n);
    }
    builder.build().map_err(|e| e.to_string())
}

/// Parallel bootstrap worker: runs all replicates on Rayon threads and sends
/// per-replicate clade maps over an MPSC channel. The main thread (caller)
/// merges results and fires `on_event` — keeping the callback on a single
/// thread with no `Sync` requirement.
///
/// Uses `std::thread::scope` so `msa` and `idx_map` can be borrowed into the
/// spawned thread without requiring `'static` lifetimes.
#[cfg(feature = "parallel")]
fn bootstrap_clade_counts_parallel<A, M>(
    msa: &MSA<A>,
    n_bootstrap_samples: usize,
    idx_map: &HashMap<String, usize>,
    n_taxa: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
) -> Result<HashMap<Vec<u8>, usize>, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use std::sync::mpsc;

    let pool = build_thread_pool(num_threads)?;
    let (tx, rx) = mpsc::channel::<Result<HashMap<Vec<u8>, usize>, String>>();
    let mut counter: HashMap<Vec<u8>, usize> = HashMap::new();

    std::thread::scope(|scope| -> Result<(), String> {
        // Spawn the Rayon work onto a background OS thread so the main thread
        // can consume from `rx` concurrently, giving real per-replicate progress.
        // `pool.install` installs the thread pool as the active pool for the
        // duration of the closure, bounding parallelism to `num_threads` workers.
        scope.spawn(|| {
            pool.install(|| {
                (0..n_bootstrap_samples)
                    .into_par_iter()
                    .for_each_with(tx, |sender, _| {
                        let result: Result<HashMap<Vec<u8>, usize>, String> = (|| {
                            let tree = msa
                                .bootstrap()?
                                .into_dist::<M>()
                                .neighbor_joining()
                                .expect("NJ bootstrap iteration failed");
                            let mut local = HashMap::new();
                            count_clades(&tree, idx_map, n_taxa, &mut local)?;
                            Ok(local)
                        })();
                        // Ignore send errors: only occurs if receiver was dropped,
                        // which cannot happen while we are in the recv loop below.
                        let _ = sender.send(result);
                    });
            });
        });

        // Main thread: receive results as they arrive, merge, and fire progress.
        for completed in 1..=n_bootstrap_samples {
            match rx.recv() {
                Ok(Ok(local)) => {
                    for (clade, count) in local {
                        *counter.entry(clade).or_insert(0) += count;
                    }
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Err("bootstrap channel closed unexpectedly".into()),
            }
            if let Some(cb) = on_event {
                cb(NJEvent::BootstrapProgress {
                    completed,
                    total: n_bootstrap_samples,
                });
            }
        }
        Ok(())
    })?;

    Ok(counter)
}

/// Performs bootstrap sampling and counts clades across bootstrap trees.
///
/// When the `parallel` feature is enabled, replicates run on the Rayon thread
/// pool via [`bootstrap_clade_counts_parallel`]; the `on_event` callback is
/// still called after each replicate from the main thread. Without the feature
/// the loop is sequential and `on_event` fires in order.
fn bootstrap_clade_counts<A, M>(
    msa: &MSA<A>,
    n_bootstrap_samples: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
) -> Result<Option<HashMap<Vec<u8>, usize>>, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    if n_bootstrap_samples == 0 {
        return Ok(None);
    }
    if let Some(cb) = on_event {
        cb(NJEvent::BootstrapStarted {
            total: n_bootstrap_samples,
        });
    }
    let idx_map: HashMap<String, usize> = msa.to_index_map();
    let n_taxa = msa.n_sequences;

    #[cfg(feature = "parallel")]
    let counter = bootstrap_clade_counts_parallel::<A, M>(
        msa, n_bootstrap_samples, &idx_map, n_taxa, on_event, num_threads,
    )?;

    #[cfg(not(feature = "parallel"))]
    let counter = {
        let _ = num_threads;
        let mut c = HashMap::new();
        for i in 0..n_bootstrap_samples {
            let tree = msa
                .bootstrap()?
                .into_dist::<M>()
                .neighbor_joining()
                .expect("NJ bootstrap iteration failed");
            count_clades(&tree, &idx_map, n_taxa, &mut c)?;
            if let Some(cb) = on_event {
                cb(NJEvent::BootstrapProgress {
                    completed: i + 1,
                    total: n_bootstrap_samples,
                });
            }
        }
        c
    };

    Ok(Some(counter))
}

/// Annotates internal nodes with bootstrap support values from `counts`.
///
/// For each internal node, computes its clade `BitVec`, looks up the count in
/// `counts`, normalises it to a percentage (`count * 100 / n_bootstrap_samples`),
/// and assigns a [`NameOrSupport::Support`] label if a matching entry is found.
/// Nodes whose clade was never observed in bootstrap replicates receive no label.
fn add_bootstrap_to_tree(
    node: &mut TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counts: &HashMap<Vec<u8>, usize>,
    n_bootstrap_samples: usize,
) -> Result<(), String> {
    if node.children.is_some() {
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(node, idx, &mut bv)?;

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            if let Some(c) = counts.get(&bv.as_raw_slice().to_vec()) {
                let pct = c * 100 / n_bootstrap_samples;
                node.label = Some(NameOrSupport::Support(pct));
            }
        }

        if let Some([l, r]) = &mut node.children {
            add_bootstrap_to_tree(l, idx, n_taxa, counts, n_bootstrap_samples)?;
            add_bootstrap_to_tree(r, idx, n_taxa, counts, n_bootstrap_samples)?;
        }
    }
    Ok(())
}

/// Validates that the MSA is non-empty and all sequences have equal length.
fn validate_msa(msa: &[SequenceObject]) -> Result<(), NJError> {
    if msa.is_empty() {
        return Err(NJError::EmptyMsa);
    }
    let expected_len = msa[0].sequence.len();
    if expected_len == 0 {
        return Err(NJError::EmptySequence);
    }
    for s in msa {
        if s.sequence.len() != expected_len {
            return Err(NJError::SequenceLengthMismatch {
                expected: expected_len,
                got: s.sequence.len(),
                identifier: s.identifier.clone(),
            });
        }
    }
    Ok(())
}

/// Heuristically detects whether the MSA contains DNA or protein sequences.
///
/// Returns [`Alphabet::DNA`] unless any sequence contains a byte that is not
/// in the DNA character set (case-insensitive), in which case
/// [`Alphabet::Protein`] is returned. The DNA set includes standard bases
/// `{A, C, G, T}`, uridine `U` (RNA), `N` (unknown), gap `-`, and all 11
/// IUPAC ambiguity codes `{R, Y, S, W, K, M, B, D, H, V}`.
fn detect_alphabet(msa: &[SequenceObject]) -> Alphabet {
    let mut is_protein = false;

    'outer: for seq in msa {
        for c in seq.sequence.bytes() {
            match c.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'-'
                | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'B' | b'D' | b'H' | b'V' => { /* still possible DNA */ }
                _ => {
                    is_protein = true;
                    break 'outer;
                }
            }
        }
    }

    if is_protein { Alphabet::Protein } else { Alphabet::DNA }
}

/// Runs distance matrix computation with model `M` on alphabet `A`.
fn run_distance_matrix<A, M>(msa: MSA<A>, num_threads: Option<usize>) -> Result<DistanceResult, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist::<M>()).into_result())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist::<M>().into_result())
    }
}

/// Runs average distance computation with model `M` on alphabet `A`.
fn run_average_distance<A, M>(msa: MSA<A>, num_threads: Option<usize>) -> Result<f64, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist::<M>()).average())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist::<M>().average())
    }
}

/// Runs NJ with model `M` on alphabet `A` and returns a Newick string.
///
/// If `n_bootstrap_samples > 0`, generates that many bootstrap replicates,
/// collects clade counts via [`bootstrap_clade_counts`], runs NJ on the
/// original distances, and annotates the tree before serialising to Newick.
fn run_nj<A, M>(
    msa: MSA<A>,
    n_bootstrap_samples: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
) -> Result<String, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    let clade_counts =
        bootstrap_clade_counts::<A, M>(&msa, n_bootstrap_samples, on_event, num_threads)?;

    if let Some(cb) = on_event {
        cb(NJEvent::ComputingDistances);
    }

    #[cfg(feature = "parallel")]
    let dist = {
        let pool = build_thread_pool(num_threads)?;
        pool.install(|| msa.into_dist::<M>())
    };
    #[cfg(not(feature = "parallel"))]
    let dist = msa.into_dist::<M>();

    if let Some(cb) = on_event {
        cb(NJEvent::RunningNJ);
    }

    let mut main_tree = dist.neighbor_joining()?;

    let newick = match clade_counts {
        Some(counts) => {
            if let Some(cb) = on_event {
                cb(NJEvent::AnnotatingBootstrap);
            }
            let main_idx_map: HashMap<String, usize> = msa.to_index_map();
            add_bootstrap_to_tree(
                &mut main_tree,
                &main_idx_map,
                msa.n_sequences,
                &counts,
                n_bootstrap_samples,
            )?;
            main_tree.to_newick()
        }
        None => main_tree.to_newick(),
    };
    Ok(newick)
}

/// Infers a phylogenetic tree from an aligned MSA and returns a Newick string.
///
/// This is the single public entry point for the library. The alphabet is
/// auto-detected from the sequences unless `conf.alphabet` is set; `conf.substitution_model`
/// must be compatible with the alphabet (see the module-level compatibility
/// table). Returns `Err` for an empty MSA, an incompatible model, or any
/// internal NJ failure.
///
/// `on_event` is called with an [`NJEvent`] at each stage of the algorithm.
/// Bootstrap progress is reported via [`NJEvent::BootstrapProgress`] after
/// each replicate. Pass `None` if event reporting is not needed.
pub fn nj(
    conf: NJConfig,
    on_event: Option<Box<dyn Fn(NJEvent)>>,
) -> Result<String, NJError> {
    let cb = on_event.as_deref();
    let num_threads = conf.num_threads;
    validate_msa(&conf.msa)?;
    let n_sites = conf.msa[0].sequence.len();
    if let Some(cb) = cb {
        cb(NJEvent::MsaValidated {
            n_sequences: conf.msa.len(),
            n_sites,
        });
    }
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    if let Some(cb) = cb {
        cb(NJEvent::AlphabetDetected {
            alphabet: alphabet.clone(),
        });
    }
    let model = conf.substitution_model;
    match alphabet {
        Alphabet::DNA => {
            let msa =
                MSA::<DNA>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::PDiff => {
                    run_nj::<DNA, PDiff>(msa, conf.n_bootstrap_samples, cb, num_threads)
                        .map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor => {
                    run_nj::<DNA, JukesCantor>(msa, conf.n_bootstrap_samples, cb, num_threads)
                        .map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Kimura2P => {
                    run_nj::<DNA, Kimura2P>(msa, conf.n_bootstrap_samples, cb, num_threads)
                        .map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Poisson => Err(NJError::IncompatibleModel {
                    model,
                    alphabet: Alphabet::DNA,
                }),
            }
        }
        Alphabet::Protein => {
            let msa =
                MSA::<Protein>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::Poisson => {
                    run_nj::<Protein, Poisson>(msa, conf.n_bootstrap_samples, cb, num_threads)
                        .map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::PDiff => {
                    run_nj::<Protein, PDiff>(msa, conf.n_bootstrap_samples, cb, num_threads)
                        .map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor | SubstitutionModel::Kimura2P => {
                    Err(NJError::IncompatibleModel {
                        model,
                        alphabet: Alphabet::Protein,
                    })
                }
            }
        }
    }
}

/// Computes pairwise distances from an aligned MSA and returns a [`DistanceResult`].
///
/// The alphabet is auto-detected from the sequences unless `conf.alphabet` is set;
/// `conf.substitution_model` must be compatible with the alphabet (see the
/// module-level compatibility table). Returns `Err` for an empty MSA, incompatible
/// model, or mismatched sequence lengths. Does not run Neighbor-Joining or bootstrapping.
pub fn distance_matrix(conf: DistConfig) -> Result<DistanceResult, NJError> {
    let num_threads = conf.num_threads;
    validate_msa(&conf.msa)?;
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    let model = conf.substitution_model;
    match alphabet {
        Alphabet::DNA => {
            let msa =
                MSA::<DNA>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::PDiff => {
                    run_distance_matrix::<DNA, PDiff>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor => {
                    run_distance_matrix::<DNA, JukesCantor>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Kimura2P => {
                    run_distance_matrix::<DNA, Kimura2P>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Poisson => Err(NJError::IncompatibleModel {
                    model,
                    alphabet: Alphabet::DNA,
                }),
            }
        }
        Alphabet::Protein => {
            let msa =
                MSA::<Protein>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::Poisson => {
                    run_distance_matrix::<Protein, Poisson>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::PDiff => {
                    run_distance_matrix::<Protein, PDiff>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor | SubstitutionModel::Kimura2P => {
                    Err(NJError::IncompatibleModel {
                        model,
                        alphabet: Alphabet::Protein,
                    })
                }
            }
        }
    }
}

/// Computes the mean of all `n*(n-1)/2` unique pairwise distances.
///
/// Same alphabet auto-detection and model–alphabet compatibility as [`nj`].
/// Returns `0.0` for fewer than 2 taxa. Returns `Err` for an empty MSA,
/// incompatible model, or mismatched sequence lengths.
pub fn average_distance(conf: DistConfig) -> Result<f64, NJError> {
    let num_threads = conf.num_threads;
    validate_msa(&conf.msa)?;
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    let model = conf.substitution_model;
    match alphabet {
        Alphabet::DNA => {
            let msa =
                MSA::<DNA>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::PDiff => {
                    run_average_distance::<DNA, PDiff>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor => {
                    run_average_distance::<DNA, JukesCantor>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Kimura2P => {
                    run_average_distance::<DNA, Kimura2P>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::Poisson => Err(NJError::IncompatibleModel {
                    model,
                    alphabet: Alphabet::DNA,
                }),
            }
        }
        Alphabet::Protein => {
            let msa =
                MSA::<Protein>::from_iter(conf.msa.into_iter().map(|s| (s.identifier, s.sequence)));
            match model {
                SubstitutionModel::Poisson => {
                    run_average_distance::<Protein, Poisson>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::PDiff => {
                    run_average_distance::<Protein, PDiff>(msa, num_threads).map_err(NJError::AlgorithmFailure)
                }
                SubstitutionModel::JukesCantor | SubstitutionModel::Kimura2P => {
                    Err(NJError::IncompatibleModel {
                        model,
                        alphabet: Alphabet::Protein,
                    })
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::DistConfig;
    use crate::models::SubstitutionModel;

    #[test]
    fn test_nj_wrapper_simple_tree() {
        // ACGTCG vs ACG-GC: pos 3 is gapped (excluded), 2 diffs out of 5 comparable
        // → distance = 0.4; two taxa → each branch = 0.2.
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
            alphabet: None,
            num_threads: None,
        };
        let newick = nj(conf, None).expect("NJ failed");
        assert_eq!(newick, "(A:0.200,B:0.200);");
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
            alphabet: None,
            num_threads: None,
        };
        let out = nj(conf, None).unwrap();
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
            alphabet: None,
            num_threads: None,
        };

        let t1 = nj(conf.clone(), None).unwrap();
        let t2 = nj(conf, None).unwrap();
        assert_eq!(t1, t2);
    }

    #[test]
    fn test_nj_wrapper_empty_msa() {
        let conf = NJConfig {
            msa: vec![],
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
        };
        let result = nj(conf, None);
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
            alphabet: None,
            num_threads: None,
        };
        let result = nj(conf, None);
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
            alphabet: None,
            num_threads: None,
        };
        let result = nj(conf, None);
        assert!(result.is_err());
    }

    // --- distance_matrix ---

    fn dist_conf(pairs: &[(&str, &str)], model: SubstitutionModel) -> DistConfig {
        DistConfig {
            msa: pairs
                .iter()
                .map(|(id, seq)| SequenceObject {
                    identifier: id.to_string(),
                    sequence: seq.to_string(),
                })
                .collect(),
            substitution_model: model,
            alphabet: None,
            num_threads: None,
        }
    }

    #[test]
    fn test_distance_matrix_names_and_shape() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        assert_eq!(result.names, vec!["A", "B"]);
        assert_eq!(result.matrix.len(), 2);
        assert_eq!(result.matrix[0].len(), 2);
        assert_eq!(result.matrix[1].len(), 2);
    }

    #[test]
    fn test_distance_matrix_diagonal_zero() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        for i in 0..3 {
            assert_eq!(result.matrix[i][i], 0.0);
        }
    }

    #[test]
    fn test_distance_matrix_symmetric() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(result.matrix[i][j], result.matrix[j][i]);
            }
        }
    }

    #[test]
    fn test_distance_matrix_pdiff_known_value() {
        // one difference at position 3 (T vs A) out of 4 → 0.25
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        assert!((result.matrix[0][1] - 0.25).abs() < 1e-12);
        assert!((result.matrix[1][0] - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_distance_matrix_identical_sequences_zero() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGT")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        assert_eq!(result.matrix[0][1], 0.0);
    }

    #[test]
    fn test_distance_matrix_jukes_cantor_dna() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")], SubstitutionModel::JukesCantor);
        let result = distance_matrix(conf).unwrap();
        // JC distance for p=0.25: -0.75 * ln(1 - 4/3 * 0.25)
        let expected = -0.75_f64 * (1.0_f64 - (4.0_f64 / 3.0) * 0.25).ln();
        assert!((result.matrix[0][1] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_kimura2p_dna() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")], SubstitutionModel::Kimura2P);
        let result = distance_matrix(conf).unwrap();
        assert_eq!(result.names, vec!["A", "B", "C"]);
        assert!(result.matrix[0][0] == 0.0);
    }

    #[test]
    fn test_distance_matrix_poisson_protein() {
        let conf = dist_conf(&[("A", "ACDEFGH"), ("B", "ACDEFGK")], SubstitutionModel::Poisson);
        let result = distance_matrix(conf).unwrap();
        // 1 diff (H vs K) out of 7: p=1/7, d=-ln(1-1/7)
        let expected = -(1.0_f64 - 1.0 / 7.0).ln();
        assert!((result.matrix[0][1] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_pdiff_protein() {
        let conf = dist_conf(&[("A", "ACDEFGH"), ("B", "ACDEFGK")], SubstitutionModel::PDiff);
        let result = distance_matrix(conf).unwrap();
        assert!((result.matrix[0][1] - 1.0 / 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_distance_matrix_empty_msa_errors() {
        let conf = DistConfig { msa: vec![], substitution_model: SubstitutionModel::PDiff, alphabet: None, num_threads: None };
        assert!(distance_matrix(conf).is_err());
    }

    #[test]
    fn test_distance_matrix_incompatible_model_errors() {
        // Poisson on DNA
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::Poisson);
        assert!(distance_matrix(conf).is_err());
        // JukesCantor on Protein
        let conf = dist_conf(&[("A", "ACDEFGH"), ("B", "ACDEFGK")], SubstitutionModel::JukesCantor);
        assert!(distance_matrix(conf).is_err());
    }

    // --- average_distance ---

    #[test]
    fn test_average_distance_identical_sequences_zero() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGT")], SubstitutionModel::PDiff);
        let avg = average_distance(conf).unwrap();
        assert_eq!(avg, 0.0);
    }

    #[test]
    fn test_average_distance_two_taxa_equals_pairwise() {
        // one difference out of 4 → 0.25
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::PDiff);
        let avg = average_distance(conf).unwrap();
        assert!((avg - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_average_distance_three_taxa_known_value() {
        // A↔B: 1/4=0.25 (T→A), A↔C: 1/4=0.25 (C→G), B↔C: 2/4=0.5 → avg = 1/3
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            SubstitutionModel::PDiff,
        );
        let avg = average_distance(conf).unwrap();
        assert!((avg - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_average_distance_jukes_cantor_dna() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::JukesCantor);
        let avg = average_distance(conf).unwrap();
        let expected = -0.75_f64 * (1.0_f64 - (4.0_f64 / 3.0) * 0.25).ln();
        assert!((avg - expected).abs() < 1e-10);
    }

    #[test]
    fn test_average_distance_empty_msa_errors() {
        let conf = DistConfig { msa: vec![], substitution_model: SubstitutionModel::PDiff, alphabet: None, num_threads: None };
        assert!(average_distance(conf).is_err());
    }

    #[test]
    fn test_average_distance_incompatible_model_errors() {
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::Poisson);
        assert!(average_distance(conf).is_err());
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
        assert_eq!(detect_alphabet(&msa), Alphabet::DNA);
    }

    #[test]
    fn test_detect_alphabet_dna_iupac() {
        // IUPAC ambiguity codes and RNA U should be detected as DNA, not protein.
        let msa = vec![SequenceObject {
            identifier: "Seq0".into(),
            sequence: "ACGTRYWSMKHBDVNU".into(),
        }];
        assert_eq!(detect_alphabet(&msa), Alphabet::DNA);
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
        assert_eq!(detect_alphabet(&msa), Alphabet::Protein);
    }
}

#[cfg(all(test, feature = "parallel"))]
mod parallel_tests {
    use super::*;
    use crate::models::SubstitutionModel;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn three_seq_dna() -> Vec<SequenceObject> {
        vec![
            SequenceObject { identifier: "A".into(), sequence: "ACGTACGT".into() },
            SequenceObject { identifier: "B".into(), sequence: "ACGCACGT".into() },
            SequenceObject { identifier: "C".into(), sequence: "ACGTACGC".into() },
        ]
    }

    #[test]
    fn test_parallel_bootstrap_returns_valid_newick() {
        let conf = NJConfig {
            msa: three_seq_dna(),
            n_bootstrap_samples: 20,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
        };
        let newick = nj(conf, None).expect("parallel NJ failed");
        assert!(newick.ends_with(';'));
        assert!(newick.contains(':'));
    }

    #[test]
    fn test_parallel_progress_fires_exactly_n_times() {
        let n = 10_usize;
        let count = Arc::new(AtomicUsize::new(0));
        let count2 = count.clone();
        let conf = NJConfig {
            msa: three_seq_dna(),
            n_bootstrap_samples: n,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
        };
        let cb: Option<Box<dyn Fn(usize, usize)>> =
            Some(Box::new(move |_, _| { count2.fetch_add(1, Ordering::SeqCst); }));
        nj(conf, cb).unwrap();
        assert_eq!(count.load(Ordering::SeqCst), n);
    }

    #[test]
    fn test_parallel_progress_last_call_is_total() {
        let n = 8_usize;
        let last = Arc::new(AtomicUsize::new(0));
        let last2 = last.clone();
        let conf = NJConfig {
            msa: three_seq_dna(),
            n_bootstrap_samples: n,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
        };
        let cb: Option<Box<dyn Fn(usize, usize)>> =
            Some(Box::new(move |completed, _| { last2.store(completed, Ordering::SeqCst); }));
        nj(conf, cb).unwrap();
        assert_eq!(last.load(Ordering::SeqCst), n);
    }
}
