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
