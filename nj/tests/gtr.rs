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
        kappa: None,
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
        kappa: None,
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

fn kappa_conf(
    msa: Vec<SequenceObject>,
    model: SubstitutionModel,
    kappa: Option<f64>,
) -> DistConfig {
    let mut conf = dist_conf(msa, model, None, None);
    conf.kappa = kappa;
    conf
}

#[test]
fn f84_with_kappa_zero_is_f81_and_kappa_changes_the_distances() {
    let msa = family(6);
    let f84 = distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::F84, Some(0.0))).unwrap();
    let f81 = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&[1.0; 6]),
        None,
    ))
    .unwrap();
    assert_eq!(f84.matrix, f81.matrix);
    let bigger = distance_matrix(kappa_conf(msa, SubstitutionModel::F84, Some(3.0))).unwrap();
    assert_ne!(bigger.matrix, f84.matrix);
    assert!(bigger.matrix[5][0].is_finite() && bigger.matrix[5][0] > 0.0);
}

#[test]
fn f84_takes_optional_frequencies_and_runs_nj_with_rate_variation() {
    let mut conf = kappa_conf(family(5), SubstitutionModel::F84, Some(1.0));
    conf.base_frequencies = Some(vec![0.3, 0.2, 0.2, 0.3]);
    assert!(distance_matrix(conf).unwrap().matrix[1][0] > 0.0);

    let mut nj = nj_conf(family(7), None, 10);
    nj.substitution_model = SubstitutionModel::F84;
    nj.kappa = Some(2.0);
    nj.gamma_shape = Some(0.7);
    nj.p_invar = Some(0.1);
    let result = run_nj(nj, None).unwrap();
    for i in 0..7 {
        assert!(
            result.newick.contains(&format!("s{i}:")),
            "{}",
            result.newick
        );
    }
    assert!(!result.newick.contains("NaN") && !result.newick.contains("inf"));
}

#[test]
fn f84_needs_kappa_and_kappa_is_only_for_f84() {
    let msa = family(4);
    let missing =
        distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::F84, None)).unwrap_err();
    assert!(
        matches!(&missing, NJError::InvalidMatrix { reason } if reason.contains("kappa")),
        "{missing:?}"
    );
    for model in [
        SubstitutionModel::Gtr,
        SubstitutionModel::Kimura2P,
        SubstitutionModel::EmpiricalGtr,
    ] {
        let stray = distance_matrix(kappa_conf(msa.clone(), model, Some(1.0))).unwrap_err();
        assert!(matches!(&stray, NJError::InvalidMatrix { .. }), "{stray:?}");
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let err = distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::F84, Some(bad)))
            .unwrap_err();
        assert!(
            matches!(&err, NJError::InvalidMatrix { reason } if reason.contains("kappa")),
            "{bad}: {err:?}"
        );
    }
    // F84 does not take the GTR rates.
    let mut conf = kappa_conf(msa, SubstitutionModel::F84, Some(1.0));
    conf.gtr_rates = Some(HKY.to_vec());
    assert!(matches!(
        distance_matrix(conf).unwrap_err(),
        NJError::InvalidMatrix { .. }
    ));
}

#[test]
fn f84_and_empirical_gtr_are_dna_only() {
    let protein = vec![
        seq("a", "ARNDCQEGHILKMFPSTWYV"),
        seq("b", "ARNDCQEGHILKMFPSTWYA"),
    ];
    let err = distance_matrix(kappa_conf(
        protein.clone(),
        SubstitutionModel::F84,
        Some(1.0),
    ))
    .unwrap_err();
    assert!(matches!(err, NJError::IncompatibleModel { .. }), "{err:?}");
    let err = distance_matrix(dist_conf(
        protein,
        SubstitutionModel::EmpiricalGtr,
        None,
        None,
    ))
    .unwrap_err();
    assert!(matches!(err, NJError::IncompatibleModel { .. }), "{err:?}");
}

/// Two sequences of `len` sites that differ in the given proportion of sites, with the changed bases
/// chosen so that the pair table is rich (all 12 off-diagonal cells occur).
fn diverged_pair(len: usize, proportion: f64) -> Vec<SequenceObject> {
    let bases = b"ACGT";
    let a: Vec<u8> = (0..len).map(|i| bases[(i * 5 + i / 7) % 4]).collect();
    let b: Vec<u8> = a
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            if ((i * 37) % 1000) as f64 / 1000.0 < proportion {
                let k = bases.iter().position(|&y| y == x).unwrap();
                bases[(k + 1 + (i / 3) % 3) % 4]
            } else {
                x
            }
        })
        .collect();
    vec![
        seq("a", &String::from_utf8(a).unwrap()),
        seq("b", &String::from_utf8(b).unwrap()),
    ]
}

#[test]
fn empirical_gtr_needs_no_options_and_agrees_with_ml_on_a_long_low_divergence_pair() {
    let msa = diverged_pair(20_000, 0.15);
    let empirical = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::EmpiricalGtr,
        None,
        None,
    ))
    .unwrap()
    .matrix[0][1];
    let ml = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&[1.0; 6]),
        None,
    ))
    .unwrap()
    .matrix[0][1];
    assert!(
        empirical.is_finite() && empirical > 0.1 && empirical < 0.4,
        "{empirical}"
    );
    assert!(
        (empirical - ml).abs() < 0.05 * ml,
        "empirical {empirical} vs ML {ml}"
    );

    for opts in [(Some(&HKY[..]), None), (None, Some(&[0.25; 4][..]))] {
        let err = distance_matrix(dist_conf(
            msa.clone(),
            SubstitutionModel::EmpiricalGtr,
            opts.0,
            opts.1,
        ))
        .unwrap_err();
        assert!(matches!(err, NJError::InvalidMatrix { .. }), "{err:?}");
    }
}

#[test]
fn empirical_gtr_rate_variation_and_nj_with_a_saturated_pair() {
    // Gamma and invariant sites are accepted and change the distance.
    let msa = diverged_pair(20_000, 0.15);
    let plain = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::EmpiricalGtr,
        None,
        None,
    ))
    .unwrap()
    .matrix[0][1];
    let mut gamma = dist_conf(msa, SubstitutionModel::EmpiricalGtr, None, None);
    gamma.gamma_shape = Some(0.5);
    let with_gamma = distance_matrix(gamma).unwrap().matrix[0][1];
    assert!(
        with_gamma > plain,
        "gamma rates imply more hidden changes: {with_gamma} vs {plain}"
    );

    // A pair with no closed form is infinite in the raw matrix and capped, with a warning, in nj().
    let swapped = vec![
        seq("a", &"AC".repeat(30)),
        seq("b", &"CA".repeat(30)),
        seq("c", &"AC".repeat(30)),
    ];
    let raw = distance_matrix(dist_conf(
        swapped.clone(),
        SubstitutionModel::EmpiricalGtr,
        None,
        None,
    ))
    .unwrap();
    assert!(raw.matrix[0][1].is_infinite() && raw.matrix[0][2] == 0.0);
    let mut conf = nj_conf(swapped, None, 5);
    conf.substitution_model = SubstitutionModel::EmpiricalGtr;
    let events = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = events.clone();
    let result = run_nj(
        conf,
        Some(Box::new(move |e| {
            sink.lock().unwrap().push(format!("{e:?}"))
        })),
    )
    .unwrap();
    assert!(
        !result.newick.contains("inf") && !result.newick.contains("NaN"),
        "{}",
        result.newick
    );
    assert_eq!(result.distance_matrix.unwrap().matrix[0][1], 10.0);
    assert!(
        events
            .lock()
            .unwrap()
            .iter()
            .any(|e| e.contains("saturated"))
    );
}

#[test]
fn hky_equals_gtr_with_the_hky_rates_and_kappa_one_is_f81() {
    let msa = family(6);
    let hky = distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::Hky, Some(4.0))).unwrap();
    let gtr = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&HKY),
        None,
    ))
    .unwrap();
    assert_eq!(hky.matrix, gtr.matrix);
    let one = distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::Hky, Some(1.0))).unwrap();
    let f81 = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Gtr,
        Some(&[1.0; 6]),
        None,
    ))
    .unwrap();
    assert_eq!(one.matrix, f81.matrix);
    // Not the same model as F84 with the same number.
    let f84 = distance_matrix(kappa_conf(msa, SubstitutionModel::F84, Some(4.0))).unwrap();
    assert_ne!(hky.matrix, f84.matrix);
}

#[test]
fn hky_takes_optional_frequencies_runs_nj_and_validates_its_options() {
    let mut conf = kappa_conf(family(5), SubstitutionModel::Hky, Some(2.0));
    conf.base_frequencies = Some(vec![0.3, 0.2, 0.2, 0.3]);
    assert!(distance_matrix(conf).unwrap().matrix[1][0] > 0.0);

    let mut nj = nj_conf(family(7), None, 10);
    nj.substitution_model = SubstitutionModel::Hky;
    nj.kappa = Some(3.0);
    nj.gamma_shape = Some(0.7);
    nj.p_invar = Some(0.1);
    let result = run_nj(nj, None).unwrap();
    for i in 0..7 {
        assert!(
            result.newick.contains(&format!("s{i}:")),
            "{}",
            result.newick
        );
    }

    let msa = family(4);
    let missing =
        distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::Hky, None)).unwrap_err();
    assert!(
        matches!(&missing, NJError::InvalidMatrix { reason } if reason.contains("kappa")),
        "{missing:?}"
    );
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let err = distance_matrix(kappa_conf(msa.clone(), SubstitutionModel::Hky, Some(bad)))
            .unwrap_err();
        assert!(
            matches!(&err, NJError::InvalidMatrix { reason } if reason.contains("kappa")),
            "{bad}: {err:?}"
        );
    }
    let mut stray = kappa_conf(msa.clone(), SubstitutionModel::Hky, Some(1.0));
    stray.gtr_rates = Some(HKY.to_vec());
    assert!(matches!(
        distance_matrix(stray).unwrap_err(),
        NJError::InvalidMatrix { .. }
    ));
    let protein = vec![
        seq("a", "ARNDCQEGHILKMFPSTWYV"),
        seq("b", "ARNDCQEGHILKMFPSTWYA"),
    ];
    let err = distance_matrix(kappa_conf(protein, SubstitutionModel::Hky, Some(2.0))).unwrap_err();
    assert!(matches!(err, NJError::IncompatibleModel { .. }), "{err:?}");
}
