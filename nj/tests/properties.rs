//! Property-based tests for the public distance/NJ API.
//!
//! These assert invariants that must hold for *any* valid alignment, catching
//! classes of bugs that hand-written examples miss: distance-matrix symmetry, a
//! zero diagonal, non-negative (or saturated) entries, and an NJ tree that
//! preserves the full input leaf set.

use nj::models::SubstitutionModel;
use nj::{DistConfig, NJConfig, SequenceObject, distance_matrix, nj as run_nj};
use proptest::prelude::*;

/// Strategy producing a rectangular DNA alignment of `3..=8` sequences, each of
/// the same random length `8..=40`, with unique identifiers.
fn dna_alignment() -> impl Strategy<Value = Vec<SequenceObject>> {
    (3usize..=8, 8usize..=40).prop_flat_map(|(n_seqs, len)| {
        // Each column is one of A/C/G/T/- (gaps allowed, but never all-gap rows
        // since at least some non-gap bases appear with overwhelming probability).
        let row = proptest::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'-']), len);
        proptest::collection::vec(row, n_seqs).prop_map(|rows| {
            rows.into_iter()
                .enumerate()
                .map(|(i, bytes)| SequenceObject {
                    identifier: format!("Seq{i}"),
                    sequence: String::from_utf8(bytes).unwrap(),
                })
                .collect()
        })
    })
}

fn dist_config(msa: Vec<SequenceObject>) -> DistConfig {
    DistConfig {
        msa,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
    }
}

fn nj_config(msa: Vec<SequenceObject>) -> NJConfig {
    NJConfig {
        msa,
        n_bootstrap_samples: 0,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: false,
        return_average_distance: false,
        gamma_shape: None,
        p_invar: None,
    }
}

proptest! {
    #[test]
    fn distance_matrix_is_symmetric_with_zero_diagonal(msa in dna_alignment()) {
        let n = msa.len();
        let result = distance_matrix(dist_config(msa)).unwrap();
        prop_assert_eq!(result.matrix.len(), n);
        for i in 0..n {
            prop_assert_eq!(result.matrix[i].len(), n);
            prop_assert_eq!(result.matrix[i][i], 0.0);
            for j in 0..n {
                // Symmetric.
                prop_assert_eq!(result.matrix[i][j], result.matrix[j][i]);
                // Non-negative, or +inf at saturation (never negative / NaN).
                let d = result.matrix[i][j];
                prop_assert!(d >= 0.0, "negative distance {d}");
                prop_assert!(!d.is_nan(), "NaN distance");
            }
        }
    }

    #[test]
    fn nj_tree_preserves_all_leaves(msa in dna_alignment()) {
        let names: Vec<String> = msa.iter().map(|s| s.identifier.clone()).collect();
        let result = run_nj(nj_config(msa), None).unwrap();
        let newick = result.newick;
        prop_assert!(newick.ends_with(';'));
        prop_assert_eq!(newick.matches('(').count(), newick.matches(')').count());
        for name in &names {
            prop_assert!(newick.contains(name.as_str()), "missing {name} in {newick}");
        }
    }
}
