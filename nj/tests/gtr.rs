//! The `Gtr` model: ML distances under explicit exchangeabilities, public API.

use nj::models::SubstitutionModel;
use nj::{
    DistConfig, NJConfig, NJError, SequenceObject, average_distance, distance_matrix, nj as run_nj,
};

fn seq(id: &str, s: &str) -> SequenceObject {
    SequenceObject {
        identifier: id.to_string(),
        sequence: s.to_string(),
    }
}

/// `n` DNA sequences of 300 sites derived from one base sequence; sequence `i` has about `4·i`%
/// of its sites substituted, and the composition is skewed towards A and T.
fn family(n: usize) -> Vec<SequenceObject> {
    let letters = b"AAAACCGGTT";
    let base: Vec<u8> = (0..300).map(|i| letters[(i * 7 + i / 3) % 10]).collect();
    (0..n)
        .map(|i| {
            let s: Vec<u8> = base
                .iter()
                .enumerate()
                .map(|(c, &b)| {
                    if (c * 31 + i * 17) % 100 < 4 * i {
                        letters[(c + i + 3) % 10]
                    } else {
                        b
                    }
                })
                .collect();
            seq(&format!("s{i}"), &String::from_utf8(s).unwrap())
        })
        .collect()
}

fn dist_conf(
    msa: Vec<SequenceObject>,
    model: SubstitutionModel,
    rates: Option<&[f64]>,
    freqs: Option<&[f64]>,
) -> DistConfig {
    DistConfig {
        msa,
        substitution_model: model,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
        matrix: None,
        gtr_rates: rates.map(<[f64]>::to_vec),
        base_frequencies: freqs.map(<[f64]>::to_vec),
    }
}

fn nj_conf(msa: Vec<SequenceObject>, rates: Option<&[f64]>, bootstrap: usize) -> NJConfig {
    NJConfig {
        msa,
        n_bootstrap_samples: bootstrap,
        substitution_model: SubstitutionModel::Gtr,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: true,
        return_average_distance: true,
        gamma_shape: None,
        p_invar: None,
        matrix: None,
        gtr_rates: rates.map(<[f64]>::to_vec),
        base_frequencies: None,
        max_distance: None,
    }
}

const HKY: [f64; 6] = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0];

#[test]
fn equal_rates_and_frequencies_reproduce_jukes_cantor() {
    let msa = family(5);
    let gtr = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&[1.0; 6]),
        Some(&[0.25; 4]),
    ))
    .unwrap();
    let jc = distance_matrix(dist_conf(msa, SubstitutionModel::JukesCantor, None, None)).unwrap();
    for i in 0..5 {
        for j in 0..5 {
            assert!(
                (gtr.matrix[i][j] - jc.matrix[i][j]).abs() < 1e-8,
                "({i},{j})"
            );
        }
    }
}

#[test]
fn pooled_frequencies_are_the_default_and_match_the_same_values_supplied() {
    let msa = family(6);
    // Pooled over every base of every sequence.
    let mut counts = [0f64; 4];
    for s in &msa {
        for b in s.sequence.bytes() {
            counts["ACGT".find(b as char).unwrap()] += 1.0;
        }
    }
    let total: f64 = counts.iter().sum();
    let pooled = counts.map(|c| c / total);
    let default = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        None,
    ))
    .unwrap();
    let explicit = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        Some(&pooled),
    ))
    .unwrap();
    assert_eq!(default.matrix, explicit.matrix);
    // Other frequencies give other distances.
    let other = distance_matrix(dist_conf(
        msa,
        SubstitutionModel::Gtr,
        Some(&HKY),
        Some(&[0.25; 4]),
    ))
    .unwrap();
    assert_ne!(default.matrix, other.matrix);
}

#[test]
fn exchangeabilities_are_relative_and_symmetric_in_their_scale() {
    let msa = family(4);
    let a = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        None,
    ))
    .unwrap();
    let scaled: Vec<f64> = HKY.iter().map(|r| r * 17.5).collect();
    let b = distance_matrix(dist_conf(msa, SubstitutionModel::Gtr, Some(&scaled), None)).unwrap();
    for i in 0..4 {
        for j in 0..4 {
            assert!((a.matrix[i][j] - b.matrix[i][j]).abs() < 1e-9);
        }
    }
}

#[test]
fn average_distance_agrees_with_the_matrix() {
    let conf = dist_conf(family(6), SubstitutionModel::Gtr, Some(&HKY), None);
    let m = distance_matrix(conf.clone()).unwrap().matrix;
    let mean = average_distance(conf).unwrap();
    let sum: f64 = (0..6)
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .map(|(i, j)| m[i][j])
        .sum();
    assert!((mean - sum / 15.0).abs() < 1e-12 && mean > 0.0);
}

#[test]
fn full_nj_run_with_bootstrap_and_rate_variation() {
    let mut conf = nj_conf(family(8), Some(&[0.9, 3.2, 1.1, 0.7, 4.4, 1.0]), 25);
    conf.gamma_shape = Some(0.8);
    conf.p_invar = Some(0.15);
    let result = run_nj(conf, None).unwrap();
    for i in 0..8 {
        assert!(
            result.newick.contains(&format!("s{i}:")),
            "{}",
            result.newick
        );
    }
    assert!(!result.newick.contains("NaN") && !result.newick.contains("inf"));
    assert!(result.average_distance.unwrap() > 0.0);
    assert_eq!(result.distance_matrix.unwrap().matrix.len(), 8);
}

#[test]
fn parameters_are_required_for_gtr_and_only_then() {
    let msa = family(4);
    let missing = run_nj(nj_conf(msa.clone(), None, 0), None).unwrap_err();
    assert!(
        matches!(&missing, NJError::InvalidMatrix { reason } if reason.contains("gtr_rates")),
        "{missing:?}"
    );

    for model in [SubstitutionModel::Kimura2P, SubstitutionModel::PDiff] {
        let stray = distance_matrix(dist_conf(msa.clone(), model, Some(&HKY), None)).unwrap_err();
        assert!(matches!(stray, NJError::InvalidMatrix { .. }), "{stray:?}");
        let stray = distance_matrix(dist_conf(
            msa.clone(),
            SubstitutionModel::JukesCantor,
            None,
            Some(&[0.25; 4]),
        ))
        .unwrap_err();
        assert!(matches!(stray, NJError::InvalidMatrix { .. }), "{stray:?}");
    }
}

#[test]
fn malformed_parameters_are_rejected_with_a_reason() {
    let msa = family(4);
    let cases: [(&[f64], Option<&[f64]>, &str); 6] = [
        (&[1.0, 2.0, 3.0], None, "6 exchangeabilities"),
        (&HKY, Some(&[0.5, 0.5]), "4 base frequencies"),
        (&HKY, Some(&[0.5, 0.5, 0.0, 0.0]), "frequencies"),
        (&[1.0, -2.0, 1.0, 1.0, 1.0, 1.0], None, "non-negative"),
        (&[0.0; 6], None, "no substitutions"),
        (&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], None, "irreducible"),
    ];
    for (rates, freqs, expect) in cases {
        let err = distance_matrix(dist_conf(
            msa.clone(),
            SubstitutionModel::Gtr,
            Some(rates),
            freqs,
        ))
        .unwrap_err();
        assert!(
            matches!(&err, NJError::InvalidMatrix { reason } if reason.contains(expect)),
            "{expect}: {err:?}"
        );
    }
}

#[test]
fn a_base_that_never_occurs_needs_explicit_frequencies() {
    let msa = vec![
        seq("a", "AACCTTAACCTT"),
        seq("b", "AACCTTAACCTA"),
        seq("c", "AACTTTAACCTT"),
    ];
    let err = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        None,
    ))
    .unwrap_err();
    assert!(
        matches!(&err, NJError::InvalidMatrix { reason } if reason.contains("'G'")),
        "{err:?}"
    );
    let ok = distance_matrix(dist_conf(
        msa,
        SubstitutionModel::Gtr,
        Some(&HKY),
        Some(&[0.3, 0.2, 0.2, 0.3]),
    ))
    .unwrap();
    assert!(ok.matrix[0][1] > 0.0);
}

#[test]
fn gtr_is_dna_only() {
    let protein = vec![
        seq("a", "ARNDCQEGHILKMFPSTWYV"),
        seq("b", "ARNDCQEGHILKMFPSTWYA"),
    ];
    let err =
        distance_matrix(dist_conf(protein, SubstitutionModel::Gtr, Some(&HKY), None)).unwrap_err();
    assert!(matches!(err, NJError::IncompatibleModel { .. }), "{err:?}");
}

#[test]
fn saturated_pairs_are_infinite_raw_and_capped_in_nj() {
    // Two sequences that differ at every site and share nothing in composition.
    let msa = vec![
        seq("a", &"ACGT".repeat(15)),
        seq("b", &"CATG".repeat(15)),
        seq("c", &"ACGT".repeat(15)),
    ];
    let raw = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        Some(&[0.25; 4]),
    ))
    .unwrap();
    assert!(raw.matrix[0][1].is_infinite());
    assert_eq!(raw.matrix[0][2], 0.0);
    let mut conf = nj_conf(msa, Some(&HKY), 0);
    conf.base_frequencies = Some(vec![0.25; 4]);
    let result = run_nj(conf, None).unwrap();
    assert!(
        !result.newick.contains("inf") && !result.newick.contains("NaN"),
        "{}",
        result.newick
    );
    assert_eq!(result.distance_matrix.unwrap().matrix[0][1], 10.0);
}
