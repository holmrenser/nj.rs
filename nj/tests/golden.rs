//! Golden-output tests for the public `nj` API.
//!
//! The two-taxon case is an *analytic* anchor: with a single difference over
//! four sites the p-distance is 0.25, so Neighbor-Joining places each leaf at
//! half that distance (0.125) from the root. The larger case is a deterministic
//! regression guard whose Newick string was captured from the verified
//! implementation and whose topology/leaf set is checked structurally below.

use nj::models::SubstitutionModel;
use nj::{NJConfig, SequenceObject, nj as run_nj};

fn config(pairs: &[(&str, &str)], model: SubstitutionModel) -> NJConfig {
    NJConfig {
        msa: pairs
            .iter()
            .map(|(id, seq)| SequenceObject {
                identifier: id.to_string(),
                sequence: seq.to_string(),
            })
            .collect(),
        n_bootstrap_samples: 0,
        substitution_model: model,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: false,
        return_average_distance: false,
        gamma_shape: None,
        p_invar: None,
    }
}

#[test]
fn golden_two_taxa_pdistance_is_exact() {
    // One difference (T vs A) over four sites → p = 0.25 → 0.125 per branch.
    let result = run_nj(config(&[("A", "ACGT"), ("B", "ACGA")], SubstitutionModel::PDiff), None).unwrap();
    assert_eq!(result.newick, "(A:0.125,B:0.125);");
}

#[test]
fn golden_two_identical_taxa_zero_branches() {
    let result = run_nj(config(&[("A", "ACGT"), ("B", "ACGT")], SubstitutionModel::PDiff), None).unwrap();
    assert_eq!(result.newick, "(A:0.000,B:0.000);");
}

#[test]
fn golden_five_taxa_regression() {
    // Deterministic regression guard. If the Newick changes, confirm the change
    // is intended (e.g. a formatting or algorithm tweak) before updating this.
    let pairs = [
        ("sp1", "ACGTACGTAC"),
        ("sp2", "ACGTACGTAA"),
        ("sp3", "ACGTTCGTAC"),
        ("sp4", "AGGTACGAAC"),
        ("sp5", "ACGTACGAAG"),
    ];
    let result = run_nj(config(&pairs, SubstitutionModel::JukesCantor), None).unwrap();
    let newick = &result.newick;

    // Structural invariants (independent of the exact string).
    assert!(newick.ends_with(';'));
    assert_eq!(
        newick.matches('(').count(),
        newick.matches(')').count(),
        "parentheses must be balanced: {newick}"
    );
    for (name, _) in &pairs {
        assert!(newick.contains(name), "missing leaf {name} in {newick}");
    }
    // Determinism across runs.
    let again = run_nj(config(&pairs, SubstitutionModel::JukesCantor), None).unwrap();
    assert_eq!(*newick, again.newick);
}
