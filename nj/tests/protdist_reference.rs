//! Cross-check against an independent implementation: PHYLIP 3.697 `protdist` (Debian package
//! 1:3.697+dfsg-4), which fits the same quantity, the maximum-likelihood distance under a fixed
//! rate matrix, with its own eigendecomposition and Newton iteration.
//!
//! The alignment is five 150-residue sequences with 5% to 41% of the columns resampled from a base
//! sequence. `protdist` was run with its JTT and PAM (Dayhoff) models, uniform rates and gamma rates with
//! coefficient of variation 1 (shape 1), and the distances printed to six decimals are the golden
//! values below (lower triangle by rows: (1,0), (2,0), (2,1), (3,0), ...). They agree with the
//! recomputed "DCMut" matrices of Kosiol and Goldman (2005), `JttDcmut` and `DayhoffDcmut`; the original
//! `Jtt` and `Dayhoff` files differ from `protdist` by a few tenths of a percent, so they are not used here.

use nj::models::SubstitutionModel;
use nj::{DistConfig, SequenceObject, distance_matrix};

const SEQUENCES: [&str; 5] = [
    "FCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYKLKVSYPNNHSNRIYPIMKAPKQVDSREIDGMMSNQPMWHCFWHFKMGCNQCGPASYQHI",
    "RCMRNWEKYRTERNFFNGNWFRYDGYRYYDRGRWCICCWDYPWQDYYEKDWNYRVWSWFLPYGKHGQGNYITSLPIVNDTFQLNSGDNWYLLTVSYPNNHSNRIYPIAKAPKQFDSREICGMMSKQPHWMCFWHFKMGCNQCGGASMHHI",
    "LHMYAWDTYVTERVFTNGNTFRYYGYRYRKRGRACIFNWDYISQHYYEKDWIYRVCSWFLSHPEIGQGNWNTPTPIVENTFQLCSKRNWALLKVSYLLNESNRIYPIMKADIQHDSREIMGMMSFQPMWWSFWHFKMGCNQTGGLSYEHW",
    "GCMANWDKYRHERNFFNGNWFAFDGYRMYMRGRWCIFCWDYIWCRYCEKDWIYRVVWPLLPYGFIAFGNYRTSLPIVNEEGIVCGFRVMALRKPSYPNNHTNMIYPNKKAPKQVFSWEKCGMMRNQEMWHVLWVFKMGPNQFGSISYQHK",
    "LCMQNWDKFRTNFNFFFGNWFRIYHYGGYNGGDWCSFPWIRYNQDAVERCHEYRKIRWFDPYWQHIQGYYIKSMAIDNKTACLCVFRKTNSLKVRYPMVHSNRIYSETKGPKWVDSMEICFGMSNQVMWHCFWHNFPGCCLCGTASYNHI",
];

const JTT: [f64; 10] = [
    0.201314, 0.498424, 0.657195, 0.578199, 0.809723, 1.159288, 0.996436, 1.202677, 1.867978,
    1.713499,
];

const PAM: [f64; 10] = [
    0.262607, 0.515714, 0.729021, 0.594067, 0.895696, 1.149418, 1.054007, 1.348652, 2.034552,
    1.923031,
];

const JTT_G1: [f64; 10] = [
    0.252119, 0.716824, 0.984005, 0.866259, 1.276535, 2.005691, 1.631199, 2.054554, 3.571949,
    3.200473,
];

const PAM_G1: [f64; 10] = [
    0.324403, 0.747143, 1.142943, 0.887496, 1.424290, 2.052796, 1.771552, 2.391612, 4.010615,
    3.583829,
];

fn distances(model: SubstitutionModel, gamma_shape: Option<f64>) -> Vec<f64> {
    let msa = SEQUENCES
        .iter()
        .enumerate()
        .map(|(i, s)| SequenceObject {
            identifier: format!("T{i}"),
            sequence: s.to_string(),
        })
        .collect();
    let result = distance_matrix(DistConfig {
        msa,
        substitution_model: model,
        alphabet: None,
        num_threads: None,
        gamma_shape,
        p_invar: None,
        matrix: None,
        gtr_rates: None,
        base_frequencies: None,
        kappa: None,
    })
    .unwrap();
    (0..SEQUENCES.len())
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .map(|(i, j)| result.matrix[i][j])
        .collect()
}

fn assert_matches(label: &str, got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len());
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        // protdist prints six decimals and stops its Newton iteration at a coarser tolerance than
        // ours (1e-8 relative), so compare at 1e-5 relative.
        let tolerance = 1e-5 * w.max(1.0);
        assert!(
            (g - w).abs() < tolerance,
            "{label}, pair {k}: nj {g} vs protdist {w}"
        );
    }
}

#[test]
fn jtt_dcmut_matches_protdist_jtt() {
    assert_matches(
        "JTT-DCMut, uniform",
        &distances(SubstitutionModel::JttDcmut, None),
        &JTT,
    );
}

#[test]
fn dayhoff_dcmut_matches_protdist_pam() {
    assert_matches(
        "DCMut-Dayhoff, uniform",
        &distances(SubstitutionModel::DayhoffDcmut, None),
        &PAM,
    );
}

#[test]
fn gamma_rates_match_protdist() {
    assert_matches(
        "JTT-DCMut, gamma 1",
        &distances(SubstitutionModel::JttDcmut, Some(1.0)),
        &JTT_G1,
    );
    assert_matches(
        "DCMut-Dayhoff, gamma 1",
        &distances(SubstitutionModel::DayhoffDcmut, Some(1.0)),
        &PAM_G1,
    );
}

#[test]
fn the_reference_distances_span_a_useful_range() {
    // Guard against a degenerate alignment: from close pairs to well beyond one substitution per site.
    let (lo, hi) = JTT
        .iter()
        .fold((f64::MAX, 0.0f64), |(l, h), &d| (l.min(d), h.max(d)));
    assert!(lo < 0.25 && hi > 0.8, "{lo} .. {hi}");
}
