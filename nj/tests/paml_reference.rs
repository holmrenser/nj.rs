//! Cross-check against an independent implementation: PAML 4.9j `baseml` (Debian package
//! 4.9j+dfsg-6), HKY85 with a fixed ts/tv rate ratio `kappa`.
//!
//! `baseml` was run once per pair on a two-sequence alignment (`runmode = 0`, tree `(a,b);`,
//! `model = 4`, `fix_kappa = 1`, no rate variation), so its base frequencies are those of the two
//! sequences; the golden value is the "tree length" it prints (the sum of the two branches, six
//! decimals). The same two-sequence alignments are given to nj.rs, whose default frequencies are pooled
//! over the sequences present. nj.rs expresses HKY85 as GTR with exchangeabilities
//! `[1, kappa, 1, 1, kappa, 1]`.
//!
//! The same pairs also give golden values for F84 with a fixed `kappa` (`model = 3`, same definition of
//! `kappa` as nj.rs) and for the empirical GTR distance, which equals the maximum-likelihood distance of
//! `baseml`'s REV model (`model = 7`, all rates free, observed frequencies) when the closed-form estimate
//! implies a valid rate matrix (non-negative off-diagonal rates). Where it does not (one of the ten pairs here), the
//! closed form and the constrained ML estimate differ, so that pair is left out.
//!
//! The five sequences have 200 sites, a composition skewed towards A and T, a 70:30 transition
//! bias in the substitutions and 5% to 29% of their sites substituted from a common base sequence.

use nj::models::SubstitutionModel;
use nj::{DistConfig, SequenceObject, distance_matrix};

const SEQUENCES: [&str; 5] = [
    "TGAAGAATACAGCGCTCAAGAGGTCATATAGTAGACCTACGAATTCCATTCTTGCTTCAGGACGCAGTATATAAAAGTAAAGAACTAAGGCACAGGTATAAAACGTAAGGAGGTATAATCATCGTGTTTCAATTTAATCACACCCCATTAATGTATATGATGATGTTTAAACAATGAGCAGGGCGGGCAGTGTTTTTTGT",
    "TGAAGAATGCAGCGCACAGAGGGTTATATAGTAGACATACGAACTCCATTCATGCTTCAGGACGTAGTATATAAAAGTAAAGAACTAAGGCACATTTATAAAACGTAAGGAGGTATAATCATCGCGTTTCAATTTAGTCACCAACCATTAATGTATATGATGATGTACAAACAATGAGCAGGGCGGGCAGTGTATTTTGT",
    "TCGAGAATATGGCGCTTAAAAAGTTATATAGTAGACCTACGAATCTCATTCGTGCTTCAGGCTGTAGCATATATAGGTGAATAACTAGGGCGCATGTATAAAACGTGAAGAGTTATAATTATCGTCTTTCAATGTAATCACAAACCAGTGATGTACGTGAGGATGTTTAAACAATGAACAGGGCGGGCAGCGTTTTATGT",
    "CGAAGAACACAGCGCTCAATGGGTTGTGTAGTAAACCCACCAAGTCCACTCATGTTTCAGAACATAATATTTAAAAGATTAGCACAAGGGCACAGGTATAAAGAGTAAGGAGGTATGATCGTCGTAATTCAGTTCAGTCACAAACCAGCAGTGGATGTGATAGTGTTTAAACAATGAGCAGGGCGAACAGTGTTTTTAGT",
    "TCAATCCTATGACGCTCATGAAATTATACGGAAGACCTAGGACCTTCATTCACACGTAAAGACGTTGCGCGTAAAAGTAAAGAGCCAAGGCTCATGTATGTAACATAAACGGCTATAATCATCGTATTTCAATGTGATTACAAACTATTAATGTATATAGCGGCGTTTGGACAGTAACTAGGACGCGCGGTGTCTCTTCT",
];

/// Pairs (i, j), i > j, in the order of the golden values: (1,0), (2,0), (2,1), (3,0), ...
const N_PAIRS: usize = 10;

const KAPPA_3_0: [f64; 10] = [
    0.108660, 0.222090, 0.270650, 0.276940, 0.309600, 0.437350, 0.431190, 0.468790, 0.480370,
    0.723430,
];

const KAPPA_0_7: [f64; 10] = [
    0.108420, 0.225830, 0.277210, 0.282590, 0.312470, 0.457570, 0.450880, 0.489950, 0.512880,
    0.789450,
];

fn nj_distance(i: usize, j: usize, kappa: f64) -> f64 {
    let msa = vec![
        SequenceObject {
            identifier: "a".into(),
            sequence: SEQUENCES[i].into(),
        },
        SequenceObject {
            identifier: "b".into(),
            sequence: SEQUENCES[j].into(),
        },
    ];
    let result = distance_matrix(DistConfig {
        msa,
        substitution_model: SubstitutionModel::Gtr,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
        matrix: None,
        gtr_rates: Some(vec![1.0, kappa, 1.0, 1.0, kappa, 1.0]),
        base_frequencies: None,
        kappa: None,
    })
    .unwrap();
    result.matrix[0][1]
}

fn check(kappa: f64, golden: &[f64]) {
    let pairs: Vec<(usize, usize)> = (0..5).flat_map(|i| (0..i).map(move |j| (i, j))).collect();
    assert_eq!(pairs.len(), N_PAIRS);
    for (&(i, j), &want) in pairs.iter().zip(golden) {
        let got = nj_distance(i, j, kappa);
        // baseml prints six decimals for each of the two branches.
        assert!(
            (got - want).abs() < 1e-5,
            "kappa {kappa}, pair ({i},{j}): nj {got} vs baseml {want}"
        );
    }
}

#[test]
fn hky_with_kappa_3_matches_baseml() {
    check(3.0, &KAPPA_3_0);
}

#[test]
fn the_hky_preset_matches_baseml() {
    // `Hky` with kappa is `Gtr` with [1, kappa, 1, 1, kappa, 1]; check the preset itself against baseml.
    let pairs: Vec<(usize, usize)> = (0..5).flat_map(|i| (0..i).map(move |j| (i, j))).collect();
    for (kappa, golden) in [(3.0, &KAPPA_3_0), (0.7, &KAPPA_0_7)] {
        for (&(i, j), &want) in pairs.iter().zip(golden.iter()) {
            let mut conf = pair_config(i, j, SubstitutionModel::Hky);
            conf.kappa = Some(kappa);
            let got = distance_matrix(conf).unwrap().matrix[0][1];
            assert!(
                (got - want).abs() < 1e-5,
                "Hky kappa {kappa}, pair ({i},{j}): nj {got} vs baseml {want}"
            );
        }
    }
}

#[test]
fn hky_with_kappa_below_one_matches_baseml() {
    check(0.7, &KAPPA_0_7);
}

#[test]
fn the_reference_distances_are_not_degenerate() {
    // From close pairs to well over half a substitution per site, and kappa changes the distances.
    let (lo, hi) = KAPPA_3_0
        .iter()
        .fold((f64::MAX, 0.0f64), |(l, h), &d| (l.min(d), h.max(d)));
    assert!(lo < 0.15 && hi > 0.5, "{lo} .. {hi}");
    let largest_change = KAPPA_3_0
        .iter()
        .zip(&KAPPA_0_7)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        largest_change > 0.03,
        "kappa should matter for the distant pairs: {largest_change}"
    );
}

const F84_KAPPA_2_0: [f64; 10] = [
    0.110430, 0.228620, 0.279960, 0.286990, 0.322960, 0.461320, 0.456950, 0.500020, 0.510380,
    0.796110,
];

const F84_KAPPA_0_5: [f64; 10] = [
    0.107630, 0.218820, 0.266070, 0.272260, 0.302800, 0.427030, 0.421730, 0.457540, 0.470710,
    0.702630,
];

/// `baseml` REV (all rates free, observed frequencies) on the pairs where the closed-form estimate implies a valid
/// rate matrix; `(i, j, distance)`.
const REV_EMBEDDABLE: [(usize, usize, f64); 9] = [
    (2, 0, 0.227170),
    (2, 1, 0.277960),
    (3, 0, 0.285680),
    (3, 1, 0.314500),
    (3, 2, 0.466680),
    (4, 0, 0.464270),
    (4, 1, 0.516020),
    (4, 2, 0.539300),
    (4, 3, 0.933200),
];

fn pair_config(i: usize, j: usize, model: SubstitutionModel) -> DistConfig {
    DistConfig {
        msa: vec![
            SequenceObject {
                identifier: "a".into(),
                sequence: SEQUENCES[i].into(),
            },
            SequenceObject {
                identifier: "b".into(),
                sequence: SEQUENCES[j].into(),
            },
        ],
        substitution_model: model,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
        matrix: None,
        gtr_rates: None,
        base_frequencies: None,
        kappa: None,
    }
}

fn check_f84(kappa: f64, golden: &[f64]) {
    let pairs: Vec<(usize, usize)> = (0..5).flat_map(|i| (0..i).map(move |j| (i, j))).collect();
    for (&(i, j), &want) in pairs.iter().zip(golden) {
        let mut conf = pair_config(i, j, SubstitutionModel::F84);
        conf.kappa = Some(kappa);
        let got = distance_matrix(conf).unwrap().matrix[0][1];
        assert!(
            (got - want).abs() < 1e-5,
            "F84 kappa {kappa}, pair ({i},{j}): nj {got} vs baseml {want}"
        );
    }
}

#[test]
fn f84_matches_baseml() {
    check_f84(2.0, &F84_KAPPA_2_0);
    check_f84(0.5, &F84_KAPPA_0_5);
}

#[test]
fn f84_with_kappa_zero_is_f81_and_differs_from_hky_with_the_same_number() {
    // kappa = 0 makes every exchangeability 1: F81, i.e. Gtr with equal rates and pooled frequencies.
    let mut f84 = pair_config(4, 1, SubstitutionModel::F84);
    f84.kappa = Some(0.0);
    let mut f81 = pair_config(4, 1, SubstitutionModel::Gtr);
    f81.gtr_rates = Some(vec![1.0; 6]);
    let (a, b) = (
        distance_matrix(f84).unwrap().matrix[0][1],
        distance_matrix(f81).unwrap().matrix[0][1],
    );
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    // F84's kappa is not HKY85's: the same number gives a different distance.
    let mut hky = pair_config(4, 1, SubstitutionModel::Gtr);
    hky.gtr_rates = Some(vec![1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    let mut f84 = pair_config(4, 1, SubstitutionModel::F84);
    f84.kappa = Some(2.0);
    let (a, b) = (
        distance_matrix(f84).unwrap().matrix[0][1],
        distance_matrix(hky).unwrap().matrix[0][1],
    );
    assert!((a - b).abs() > 1e-3, "{a} vs {b}");
}

#[test]
fn empirical_gtr_equals_baseml_rev_where_the_estimate_is_a_valid_rate_matrix() {
    for &(i, j, want) in &REV_EMBEDDABLE {
        let got = distance_matrix(pair_config(i, j, SubstitutionModel::EmpiricalGtr))
            .unwrap()
            .matrix[0][1];
        assert!(
            (got - want).abs() < 1e-5,
            "pair ({i},{j}): nj {got} vs baseml REV {want}"
        );
    }
}
