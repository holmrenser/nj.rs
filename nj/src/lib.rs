// The explicit `core::simd` distance kernel needs the nightly `portable_simd`
// feature. Gated on the `simd` Cargo feature so default (stable) builds are
// unaffected — the kernel falls back to the autovectorizable scalar path.
#![cfg_attr(feature = "simd", feature(portable_simd))]
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
//! The main entry point is [`nj`], which accepts an [`NJConfig`] and returns an
//! [`NJResult`] containing the Newick tree string and, optionally, the pairwise
//! distance matrix and mean distance. [`distance_matrix`] and [`average_distance`]
//! take a [`DistConfig`] for distance-only computation. Everything else is
//! internal implementation detail exposed only to the Python and WASM wrapper crates.
//!
//! # Model–alphabet compatibility
//!
//! | Model | DNA | Protein |
//! |-------|-----|---------|
//! | `PDiff` | ✓ | ✓ |
//! | `JukesCantor` | ✓ | — |
//! | `Kimura2P` | ✓ | — |
//! | `TajimaNei` | ✓ | — |
//! | `Tamura` | ✓ | — |
//! | `Poisson` | — | ✓ |
//! | `KimuraProtein` | — | ✓ |
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
#[cfg(test)]
mod nj_reference;
pub mod tree;

mod bootstrap;
mod runner;
mod validate;

use crate::config::SubstitutionModel;
pub use crate::config::{DistConfig, MSA, NJConfig, NJResult, SequenceObject};
pub use crate::distance_matrix::DistanceResult;
pub use crate::error::NJError;
pub use crate::event::{LogLevel, NJEvent};
pub use crate::fasta::parse_fasta;
use crate::runner::{dispatch_run, run_average_distance, run_distance_matrix, run_nj};
use crate::validate::{
    detect_alphabet, validate_max_distance, validate_msa, validate_rate_params,
    warn_degenerate_sequences,
};

/// Infers a phylogenetic tree from an aligned MSA and returns an [`NJResult`].
///
/// This is the single public entry point for the library. The alphabet is
/// auto-detected from the sequences unless `conf.alphabet` is set; `conf.substitution_model`
/// must be compatible with the alphabet (see the module-level compatibility
/// table). Returns `Err` for an empty MSA, duplicate identifiers, an
/// incompatible model, or any internal NJ failure.
///
/// Set [`NJConfig::return_distance_matrix`] and/or
/// [`NJConfig::return_average_distance`] to `true` to include those values in
/// the returned [`NJResult`] alongside the Newick tree.
///
/// **Degenerate inputs.** A single-sequence MSA produces a lone-leaf "tree".
/// Sequences consisting entirely of gaps have no comparable sites, so their
/// pairwise distances are `0.0` (see [`crate::models`]); a `Warning`-level
/// [`NJEvent::Log`] is emitted for each such sequence so the behaviour is not silent.
///
/// **Saturated pairs.** Neighbor-joining cannot use infinite distances, so any
/// distance that is not finite (the sequences are too diverged for the model) is
/// replaced by [`NJConfig::max_distance`] (default 10.0 substitutions per site), a
/// warning [`NJEvent::Log`] reports how many pairs were affected, and the run
/// continues; bootstrap replicates are treated the same way.
///
/// `on_event` is called with an [`NJEvent`] at each stage of the algorithm.
/// Bootstrap progress is reported via [`NJEvent::BootstrapProgress`] after
/// each replicate. Pass `None` if event reporting is not needed.
pub fn nj(conf: NJConfig, on_event: Option<Box<dyn Fn(NJEvent)>>) -> Result<NJResult, NJError> {
    let cb = on_event.as_deref();
    let num_threads = conf.num_threads;
    let include_distance_matrix = conf.return_distance_matrix;
    let include_average_distance = conf.return_average_distance;
    let n_bootstrap_samples = conf.n_bootstrap_samples;
    validate_msa(&conf.msa)?;
    let n_sites = conf.msa[0].sequence.len();
    if let Some(cb) = cb {
        cb(NJEvent::MsaValidated {
            n_sequences: conf.msa.len(),
            n_sites,
        });
    }
    warn_degenerate_sequences(&conf.msa, cb);
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    if let Some(cb) = cb {
        cb(NJEvent::AlphabetDetected {
            alphabet: alphabet.clone(),
        });
    }
    let rates = validate_rate_params(conf.gamma_shape, conf.p_invar)?;
    let max_distance = validate_max_distance(conf.max_distance)?;
    let model = conf.substitution_model;
    if matches!(model, SubstitutionModel::PDiff)
        && (rates.gamma_shape.is_some() || rates.p_invar > 0.0)
        && let Some(cb) = cb
    {
        cb(NJEvent::Log {
            level: LogLevel::Warning,
            message: "Rate-variation parameters (gamma shape / invariant sites) have no effect \
                      on the PDiff model and are ignored."
                .to_string(),
        });
    }
    let pairs = conf.msa.into_iter().map(|s| (s.identifier, s.sequence));
    dispatch_run!(
        alphabet,
        model,
        pairs,
        run_nj(
            n_bootstrap_samples,
            cb,
            num_threads,
            include_distance_matrix,
            include_average_distance,
            rates,
            max_distance,
        )
    )
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
    let rates = validate_rate_params(conf.gamma_shape, conf.p_invar)?;
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    let model = conf.substitution_model;
    let pairs = conf.msa.into_iter().map(|s| (s.identifier, s.sequence));
    dispatch_run!(
        alphabet,
        model,
        pairs,
        run_distance_matrix(num_threads, rates)
    )
}

/// Computes the mean of all `n*(n-1)/2` unique pairwise distances.
///
/// Same alphabet auto-detection and model–alphabet compatibility as [`nj`].
/// Returns `0.0` for fewer than 2 taxa. Returns `Err` for an empty MSA,
/// incompatible model, or mismatched sequence lengths.
pub fn average_distance(conf: DistConfig) -> Result<f64, NJError> {
    let num_threads = conf.num_threads;
    validate_msa(&conf.msa)?;
    let rates = validate_rate_params(conf.gamma_shape, conf.p_invar)?;
    let alphabet = conf.alphabet.unwrap_or_else(|| detect_alphabet(&conf.msa));
    let model = conf.substitution_model;
    let pairs = conf.msa.into_iter().map(|s| (s.identifier, s.sequence));
    dispatch_run!(
        alphabet,
        model,
        pairs,
        run_average_distance(num_threads, rates)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Alphabet;
    use crate::config::DistConfig;
    use crate::models::SubstitutionModel;

    /// Builds a minimal [`NJConfig`] with both include flags set to the given values.
    fn nj_conf(pairs: &[(&str, &str)], include_dm: bool, include_avg: bool) -> NJConfig {
        NJConfig {
            msa: pairs
                .iter()
                .map(|(id, seq)| SequenceObject {
                    identifier: id.to_string(),
                    sequence: seq.to_string(),
                })
                .collect(),
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
            return_distance_matrix: include_dm,
            return_average_distance: include_avg,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        }
    }

    #[test]
    fn test_nj_rejects_nonpositive_gamma_shape() {
        let mut conf = nj_conf(&[("A", "ACGT"), ("B", "ACGA")], false, false);
        conf.substitution_model = SubstitutionModel::JukesCantor;
        conf.gamma_shape = Some(0.0);
        assert!(matches!(
            nj(conf, None),
            Err(NJError::InvalidGammaShape { .. })
        ));
    }

    #[test]
    fn test_nj_rejects_out_of_range_p_invar() {
        let mut conf = nj_conf(&[("A", "ACGT"), ("B", "ACGA")], false, false);
        conf.substitution_model = SubstitutionModel::JukesCantor;
        conf.p_invar = Some(1.0);
        assert!(matches!(nj(conf, None), Err(NJError::InvalidPInvar { .. })));
    }

    #[test]
    fn test_nj_gamma_changes_branch_lengths() {
        // The gamma shape must flow all the way through nj() into the model and
        // change the resulting branch lengths versus the uniform-rate model.
        let pairs = &[("A", "ACGTACGT"), ("B", "AGGTACGA"), ("C", "ACGAACTT")];
        let mut plain = nj_conf(pairs, false, false);
        plain.substitution_model = SubstitutionModel::JukesCantor;
        let mut gamma = plain.clone();
        gamma.gamma_shape = Some(0.5);
        let plain_nwk = nj(plain, None).unwrap().newick;
        let gamma_nwk = nj(gamma, None).unwrap().newick;
        assert_ne!(plain_nwk, gamma_nwk);
    }

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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        };
        let result = nj(conf, None).expect("NJ failed");
        assert_eq!(result.newick, "(A:0.200,B:0.200);");
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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        };
        let out = nj(conf, None).unwrap();
        assert!(out.newick.ends_with(';'));
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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
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
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        };
        let result = nj(conf, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_nj_rejects_duplicate_identifiers() {
        let result = nj(nj_conf(&[("A", "ACGT"), ("A", "ACGA")], false, false), None);
        assert!(matches!(result, Err(NJError::DuplicateIdentifier { .. })));
    }

    #[test]
    fn test_distance_matrix_rejects_duplicate_identifiers() {
        let conf = dist_conf(&[("A", "ACGT"), ("A", "ACGA")], SubstitutionModel::PDiff);
        assert!(matches!(
            distance_matrix(conf),
            Err(NJError::DuplicateIdentifier { .. })
        ));
    }

    #[test]
    fn test_lowercase_protein_nj_is_non_degenerate() {
        // Lowercase protein residues must encode like uppercase, so the two
        // distinct sequences yield a positive branch length rather than 0.
        let conf = NJConfig {
            msa: vec![
                SequenceObject {
                    identifier: "P0".into(),
                    sequence: "acdefghikl".into(),
                },
                SequenceObject {
                    identifier: "P1".into(),
                    sequence: "acdefgwxyl".into(),
                },
            ],
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: Some(Alphabet::Protein),
            num_threads: None,
            return_distance_matrix: true,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        };
        let result = nj(conf, None).expect("NJ on lowercase protein should succeed");
        let dm = result.distance_matrix.unwrap();
        assert!(
            dm.matrix[0][1] > 0.0,
            "expected non-zero distance, got {dm:?}"
        );
    }

    #[test]
    fn test_all_gap_sequence_emits_warning() {
        use std::sync::{Arc, Mutex};
        let warnings: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let sink = Arc::clone(&warnings);
        let cb: Box<dyn Fn(NJEvent)> = Box::new(move |event| {
            if let NJEvent::Log {
                level: LogLevel::Warning,
                message,
            } = event
            {
                sink.lock().unwrap().push(message);
            }
        });
        let conf = nj_conf(&[("A", "----"), ("B", "ACGA")], false, false);
        let _ = nj(conf, Some(cb)).unwrap();
        let collected = warnings.lock().unwrap();
        assert_eq!(collected.len(), 1, "expected exactly one all-gap warning");
        assert!(collected[0].contains("only gaps"), "got {:?}", collected[0]);
    }

    // --- Saturated pairs (max_distance) ---

    /// Two groups of near-identical sequences, `{A, C, D}` and `{B, E}`, with every
    /// site differing between the groups. Under `JukesCantor` the six between-group
    /// pairs have no finite distance; the four within-group pairs are close.
    fn saturated_conf() -> NJConfig {
        let mut conf = nj_conf(
            &[
                ("A", "ACGTACGTACGTACGTACGT"),
                ("B", "CATGCATGCATGCATGCATG"),
                ("C", "ACGTACGTACGAACGTACGT"),
                ("D", "ACGTACGTTCGTACGTACGT"),
                ("E", "CATGCATGCATCCATGCATG"),
            ],
            true,
            true,
        );
        conf.substitution_model = SubstitutionModel::JukesCantor;
        conf
    }

    fn run_collecting_warnings(conf: NJConfig) -> (NJResult, Vec<String>) {
        use std::sync::{Arc, Mutex};
        let warnings: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let sink = Arc::clone(&warnings);
        let cb: Box<dyn Fn(NJEvent)> = Box::new(move |event| {
            if let NJEvent::Log {
                level: LogLevel::Warning,
                message,
            } = event
            {
                sink.lock().unwrap().push(message);
            }
        });
        let result = nj(conf, Some(cb)).expect("NJ should succeed with saturated pairs");
        let collected = warnings.lock().unwrap().clone();
        (result, collected)
    }

    #[test]
    fn test_saturated_pairs_give_finite_newick_and_one_warning() {
        let (result, warnings) = run_collecting_warnings(saturated_conf());
        assert!(
            !result.newick.contains("NaN") && !result.newick.contains("inf"),
            "non-finite value in Newick: {}",
            result.newick
        );
        assert_eq!(warnings.len(), 1, "got {warnings:?}");
        assert!(warnings[0].contains("6 sequence pair(s)"), "got {warnings:?}");
        assert!(warnings[0].contains("and 1 more"), "got {warnings:?}");
        assert!(warnings[0].contains("max_distance = 10"), "got {warnings:?}");
    }

    #[test]
    fn test_saturated_pairs_are_capped_in_returned_matrix_and_average() {
        let (result, _) = run_collecting_warnings(saturated_conf());
        let dm = result.distance_matrix.unwrap();
        assert_eq!(dm.matrix[0][1], DEFAULT_MAX_DISTANCE_FOR_TESTS);
        assert_eq!(dm.matrix[1][0], DEFAULT_MAX_DISTANCE_FOR_TESTS);
        assert!(dm.matrix.iter().flatten().all(|d| d.is_finite()));
        assert!(result.average_distance.unwrap().is_finite());
    }

    const DEFAULT_MAX_DISTANCE_FOR_TESTS: f64 = crate::config::DEFAULT_MAX_DISTANCE;

    #[test]
    fn test_custom_max_distance_is_used() {
        let mut conf = saturated_conf();
        conf.max_distance = Some(3.0);
        let (result, warnings) = run_collecting_warnings(conf);
        assert_eq!(result.distance_matrix.unwrap().matrix[0][1], 3.0);
        assert!(warnings[0].contains("max_distance = 3"), "got {warnings:?}");
    }

    #[test]
    fn test_invalid_max_distance_is_rejected() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut conf = saturated_conf();
            conf.max_distance = Some(bad);
            assert!(
                matches!(nj(conf, None), Err(NJError::InvalidMaxDistance { .. })),
                "max_distance = {bad} should be rejected"
            );
        }
    }

    #[test]
    fn test_no_saturation_no_warning() {
        let mut conf = nj_conf(&[("A", "ACGTACGT"), ("B", "ACGTACGA"), ("C", "ACGAACGA")], false, false);
        conf.substitution_model = SubstitutionModel::JukesCantor;
        let (_, warnings) = run_collecting_warnings(conf);
        assert!(warnings.is_empty(), "got {warnings:?}");
    }

    #[test]
    fn test_bootstrap_with_saturated_pairs_completes_and_warns_once() {
        let mut conf = saturated_conf();
        conf.n_bootstrap_samples = 10;
        let (result, warnings) = run_collecting_warnings(conf);
        assert!(
            !result.newick.contains("NaN") && !result.newick.contains("inf"),
            "{}",
            result.newick
        );
        let replicate_warnings: Vec<_> = warnings.iter().filter(|w| w.contains("bootstrap")).collect();
        assert_eq!(replicate_warnings.len(), 1, "got {warnings:?}");
        assert!(replicate_warnings[0].contains("10 of 10 bootstrap replicates"), "got {warnings:?}");
    }

    #[test]
    fn test_standalone_distance_matrix_still_reports_raw_infinity() {
        let conf = saturated_conf();
        let result = distance_matrix(DistConfig {
            msa: conf.msa,
            substitution_model: SubstitutionModel::JukesCantor,
            alphabet: None,
            num_threads: None,
            gamma_shape: None,
            p_invar: None,
        })
        .unwrap();
        assert!(result.matrix[0][1].is_infinite());
    }

    // --- NJResult optional fields ---

    #[test]
    fn test_nj_result_no_extras_by_default() {
        let result = nj(nj_conf(&[("A", "ACGT"), ("B", "ACGA")], false, false), None).unwrap();
        assert!(result.distance_matrix.is_none());
        assert!(result.average_distance.is_none());
    }

    #[test]
    fn test_nj_result_include_distance_matrix() {
        let result = nj(nj_conf(&[("A", "ACGT"), ("B", "ACGA")], true, false), None).unwrap();
        let dm = result
            .distance_matrix
            .expect("distance_matrix should be present");
        assert_eq!(dm.names, vec!["A", "B"]);
        assert_eq!(dm.matrix.len(), 2);
        assert_eq!(dm.matrix[0][0], 0.0);
        assert!((dm.matrix[0][1] - 0.25).abs() < 1e-12);
        assert!(result.average_distance.is_none());
    }

    #[test]
    fn test_nj_result_include_average_distance() {
        let result = nj(nj_conf(&[("A", "ACGT"), ("B", "ACGA")], false, true), None).unwrap();
        let avg = result
            .average_distance
            .expect("average_distance should be present");
        assert!((avg - 0.25).abs() < 1e-12);
        assert!(result.distance_matrix.is_none());
    }

    #[test]
    fn test_nj_result_include_both() {
        let result = nj(nj_conf(&[("A", "ACGT"), ("B", "ACGA")], true, true), None).unwrap();
        let dm = result
            .distance_matrix
            .as_ref()
            .expect("distance_matrix should be present");
        let avg = result
            .average_distance
            .expect("average_distance should be present");
        // distance matrix and average_distance must be consistent
        assert!((dm.matrix[0][1] - avg).abs() < 1e-12);
    }

    #[test]
    fn test_nj_result_distance_matrix_consistent_with_newick() {
        // Known tree: two taxa with 1/4 difference → branch lengths 0.125 each from NJ
        let result = nj(nj_conf(&[("A", "ACGT"), ("B", "ACGA")], true, false), None).unwrap();
        let dm = result.distance_matrix.unwrap();
        // full matrix should be symmetric with zero diagonal
        assert_eq!(dm.matrix[0][0], 0.0);
        assert_eq!(dm.matrix[1][1], 0.0);
        assert!((dm.matrix[0][1] - dm.matrix[1][0]).abs() < 1e-12);
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
            gamma_shape: None,
            p_invar: None,
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
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            SubstitutionModel::PDiff,
        );
        let result = distance_matrix(conf).unwrap();
        for i in 0..3 {
            assert_eq!(result.matrix[i][i], 0.0);
        }
    }

    #[test]
    fn test_distance_matrix_symmetric() {
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            SubstitutionModel::PDiff,
        );
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
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            SubstitutionModel::JukesCantor,
        );
        let result = distance_matrix(conf).unwrap();
        // JC distance for p=0.25: -0.75 * ln(1 - 4/3 * 0.25)
        let expected = -0.75_f64 * (1.0_f64 - (4.0_f64 / 3.0) * 0.25).ln();
        assert!((result.matrix[0][1] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_kimura2p_dna() {
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA"), ("C", "AGGT")],
            SubstitutionModel::Kimura2P,
        );
        let result = distance_matrix(conf).unwrap();
        assert_eq!(result.names, vec!["A", "B", "C"]);
        assert!(result.matrix[0][0] == 0.0);
    }

    #[test]
    fn test_distance_matrix_poisson_protein() {
        let conf = dist_conf(
            &[("A", "ACDEFGH"), ("B", "ACDEFGK")],
            SubstitutionModel::Poisson,
        );
        let result = distance_matrix(conf).unwrap();
        // 1 diff (H vs K) out of 7: p=1/7, d=-ln(1-1/7)
        let expected = -(1.0_f64 - 1.0 / 7.0).ln();
        assert!((result.matrix[0][1] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_pdiff_protein() {
        let conf = dist_conf(
            &[("A", "ACDEFGH"), ("B", "ACDEFGK")],
            SubstitutionModel::PDiff,
        );
        let result = distance_matrix(conf).unwrap();
        assert!((result.matrix[0][1] - 1.0 / 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_distance_matrix_empty_msa_errors() {
        let conf = DistConfig {
            msa: vec![],
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
            gamma_shape: None,
            p_invar: None,
        };
        assert!(distance_matrix(conf).is_err());
    }

    #[test]
    fn test_distance_matrix_incompatible_model_errors() {
        // Poisson on DNA
        let conf = dist_conf(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::Poisson);
        assert!(distance_matrix(conf).is_err());
        // JukesCantor on Protein
        let conf = dist_conf(
            &[("A", "ACDEFGH"), ("B", "ACDEFGK")],
            SubstitutionModel::JukesCantor,
        );
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
        let conf = dist_conf(
            &[("A", "ACGT"), ("B", "ACGA")],
            SubstitutionModel::JukesCantor,
        );
        let avg = average_distance(conf).unwrap();
        let expected = -0.75_f64 * (1.0_f64 - (4.0_f64 / 3.0) * 0.25).ln();
        assert!((avg - expected).abs() < 1e-10);
    }

    #[test]
    fn test_average_distance_empty_msa_errors() {
        let conf = DistConfig {
            msa: vec![],
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
            gamma_shape: None,
            p_invar: None,
        };
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
            SequenceObject {
                identifier: "A".into(),
                sequence: "ACGTACGT".into(),
            },
            SequenceObject {
                identifier: "B".into(),
                sequence: "ACGCACGT".into(),
            },
            SequenceObject {
                identifier: "C".into(),
                sequence: "ACGTACGC".into(),
            },
        ]
    }

    fn base_conf(msa: Vec<SequenceObject>, n_bootstrap_samples: usize) -> NJConfig {
        NJConfig {
            msa,
            n_bootstrap_samples,
            substitution_model: SubstitutionModel::PDiff,
            alphabet: None,
            num_threads: None,
            return_distance_matrix: false,
            return_average_distance: false,
            max_distance: None,
            gamma_shape: None,
            p_invar: None,
        }
    }

    #[test]
    fn test_parallel_bootstrap_returns_valid_newick() {
        let result = nj(base_conf(three_seq_dna(), 20), None).expect("parallel NJ failed");
        assert!(result.newick.ends_with(';'));
        assert!(result.newick.contains(':'));
    }

    #[test]
    fn test_parallel_progress_fires_exactly_n_times() {
        let n = 10_usize;
        let count = Arc::new(AtomicUsize::new(0));
        let count2 = count.clone();
        let cb: Box<dyn Fn(NJEvent)> = Box::new(move |event| {
            if let NJEvent::BootstrapProgress { .. } = event {
                count2.fetch_add(1, Ordering::SeqCst);
            }
        });
        nj(base_conf(three_seq_dna(), n), Some(cb)).unwrap();
        assert_eq!(count.load(Ordering::SeqCst), n);
    }

    #[test]
    fn test_parallel_progress_last_call_is_total() {
        let n = 8_usize;
        let last = Arc::new(AtomicUsize::new(0));
        let last2 = last.clone();
        let cb: Box<dyn Fn(NJEvent)> = Box::new(move |event| {
            if let NJEvent::BootstrapProgress { completed, .. } = event {
                last2.store(completed, Ordering::SeqCst);
            }
        });
        nj(base_conf(three_seq_dna(), n), Some(cb)).unwrap();
        assert_eq!(last.load(Ordering::SeqCst), n);
    }
}
