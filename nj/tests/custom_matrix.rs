//! The `Custom` model: ML distances under a user-supplied protein rate matrix.

use nj::alphabet::Alphabet;
use nj::models::SubstitutionModel;
use nj::{
    DistConfig, NJConfig, NJError, SequenceObject, average_distance, distance_matrix, nj as run_nj,
};

/// 190 equal exchangeabilities and 20 equal frequencies: the 20-state equal-rates model.
const EQUAL_RATES: &str = include_str!("data/equal_rates_protein.dat");

fn seq(id: &str, s: &str) -> SequenceObject {
    SequenceObject {
        identifier: id.to_string(),
        sequence: s.to_string(),
    }
}

/// `n` protein sequences of 120 columns derived from one base sequence; sequence `i` has about
/// `4·i`% of its columns substituted.
fn family(n: usize) -> Vec<SequenceObject> {
    let letters = b"ARNDCQEGHILKMFPSTWYV";
    let base: Vec<u8> = (0..120).map(|i| letters[(i * 7 + i / 3) % 20]).collect();
    (0..n)
        .map(|i| {
            let s: Vec<u8> = base
                .iter()
                .enumerate()
                .map(|(c, &b)| {
                    if (c * 31 + i * 17) % 100 < 4 * i {
                        letters[(c + i + 5) % 20]
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
    matrix: Option<&str>,
) -> DistConfig {
    DistConfig {
        msa,
        substitution_model: model,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
        matrix: matrix.map(str::to_string),
        gtr_rates: None,
        base_frequencies: None,
    }
}

fn nj_conf(msa: Vec<SequenceObject>, matrix: Option<&str>, bootstrap: usize) -> NJConfig {
    NJConfig {
        msa,
        n_bootstrap_samples: bootstrap,
        substitution_model: SubstitutionModel::Custom,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: true,
        return_average_distance: true,
        gamma_shape: None,
        p_invar: None,
        matrix: matrix.map(str::to_string),
        gtr_rates: None,
        base_frequencies: None,
        max_distance: None,
    }
}

#[test]
fn equal_rates_matrix_reproduces_the_closed_form() {
    // 100 columns, 15 of them different: p = 0.15, d = −(19/20)·ln(1 − p/(19/20)).
    let a = "A".repeat(100);
    let b = format!("{}{}", "W".repeat(15), "A".repeat(85));
    // Mostly-`A` sequences would be auto-detected as DNA, so name the alphabet.
    let mut conf = dist_conf(
        vec![seq("a", &a), seq("b", &b)],
        SubstitutionModel::Custom,
        Some(EQUAL_RATES),
    );
    conf.alphabet = Some(Alphabet::Protein);
    let d = distance_matrix(conf).unwrap();
    let want = -(19.0_f64 / 20.0) * (1.0_f64 - 0.15 / (19.0_f64 / 20.0)).ln();
    assert!(
        (d.matrix[0][1] - want).abs() < 1e-8,
        "{} vs {want}",
        d.matrix[0][1]
    );
    assert_eq!(d.matrix[0][1], d.matrix[1][0]);
    assert_eq!(d.matrix[0][0], 0.0);
}

#[test]
fn average_distance_and_distance_matrix_agree() {
    let msa = family(6);
    let conf = dist_conf(msa, SubstitutionModel::Custom, Some(EQUAL_RATES));
    let m = distance_matrix(conf.clone()).unwrap().matrix;
    let mean = average_distance(conf).unwrap();
    let mut sum = 0.0;
    for i in 0..6 {
        for j in 0..i {
            sum += m[i][j];
        }
    }
    assert!((mean - sum / 15.0).abs() < 1e-12);
    assert!(mean > 0.0);
}

#[test]
fn full_nj_run_with_bootstrap_rate_variation_and_outputs() {
    let mut conf = nj_conf(family(8), Some(EQUAL_RATES), 25);
    conf.gamma_shape = Some(0.8);
    conf.p_invar = Some(0.1);
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
fn the_matrix_must_be_given_for_custom_and_only_then() {
    let msa = family(4);
    let missing = run_nj(nj_conf(msa.clone(), None, 0), None).unwrap_err();
    assert!(
        matches!(missing, NJError::InvalidMatrix { .. }),
        "{missing:?}"
    );
    assert_eq!(missing.code(), "InvalidMatrix");

    let mut conf = nj_conf(msa.clone(), Some(EQUAL_RATES), 0);
    conf.substitution_model = SubstitutionModel::Poisson;
    assert!(matches!(
        run_nj(conf, None).unwrap_err(),
        NJError::InvalidMatrix { .. }
    ));

    let conf = dist_conf(msa, SubstitutionModel::KimuraProtein, Some(EQUAL_RATES));
    assert!(matches!(
        distance_matrix(conf).unwrap_err(),
        NJError::InvalidMatrix { .. }
    ));
}

#[test]
fn malformed_matrices_are_rejected_with_a_reason() {
    for (text, why) in [
        ("1 2 3", "too few numbers"),
        ("not a matrix", "not numbers"),
        (&EQUAL_RATES.replace("0.05", "0"), "zero frequencies"),
    ] {
        let err = distance_matrix(dist_conf(family(3), SubstitutionModel::Custom, Some(text)))
            .unwrap_err();
        assert!(
            matches!(&err, NJError::InvalidMatrix { reason } if !reason.is_empty()),
            "{why}: {err:?}"
        );
        assert!(
            err.to_string().contains("Invalid substitution matrix"),
            "{why}"
        );
    }
}

#[test]
fn custom_is_protein_only_for_now() {
    let dna = vec![seq("a", "ACGTACGTAC"), seq("b", "ACGTACGTAA")];
    let err =
        distance_matrix(dist_conf(dna, SubstitutionModel::Custom, Some(EQUAL_RATES))).unwrap_err();
    assert!(matches!(err, NJError::IncompatibleModel { .. }), "{err:?}");
}

#[test]
fn saturated_pairs_are_capped_in_nj_and_infinite_in_the_raw_matrix() {
    // Two sequences that disagree in every column: beyond what the model can resolve.
    let letters = b"ARNDCQEGHILKMFPSTWYV";
    let a: String = (0..60).map(|i| letters[i % 20] as char).collect();
    let b: String = (0..60).map(|i| letters[(i + 7) % 20] as char).collect();
    let msa = vec![seq("a", &a), seq("b", &b), seq("c", &a)];
    let raw = distance_matrix(dist_conf(
        msa.clone(),
        SubstitutionModel::Custom,
        Some(EQUAL_RATES),
    ))
    .unwrap();
    assert!(raw.matrix[0][1].is_infinite());
    assert_eq!(raw.matrix[0][2], 0.0);

    let events = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = events.clone();
    let result = run_nj(
        nj_conf(msa, Some(EQUAL_RATES), 0),
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

const BUILTINS: [(SubstitutionModel, &str); 6] = [
    (
        SubstitutionModel::Lg,
        include_str!("../data/matrices/LG.dat"),
    ),
    (
        SubstitutionModel::Wag,
        include_str!("../data/matrices/WAG.dat"),
    ),
    (
        SubstitutionModel::Jtt,
        include_str!("../data/matrices/JTT.dat"),
    ),
    (
        SubstitutionModel::JttDcmut,
        include_str!("../data/matrices/DCMut-JTT.dat"),
    ),
    (
        SubstitutionModel::Dayhoff,
        include_str!("../data/matrices/Dayhoff.dat"),
    ),
    (
        SubstitutionModel::DayhoffDcmut,
        include_str!("../data/matrices/DCMut-Dayhoff.dat"),
    ),
];

#[test]
fn builtin_models_equal_custom_with_the_same_file() {
    // The embedded matrix and the same file passed at run time must give identical distances.
    for (model, text) in BUILTINS {
        let builtin = distance_matrix(dist_conf(family(6), model.clone(), None)).unwrap();
        let custom =
            distance_matrix(dist_conf(family(6), SubstitutionModel::Custom, Some(text))).unwrap();
        assert_eq!(builtin.matrix, custom.matrix, "{model:?}");
        assert!(
            builtin.matrix[0][5] > 0.0 && builtin.matrix[0][5].is_finite(),
            "{model:?}"
        );
    }
}

#[test]
fn builtin_models_run_nj_with_bootstrap() {
    for (model, _) in BUILTINS {
        let mut conf = nj_conf(family(7), None, 10);
        conf.substitution_model = model.clone();
        conf.gamma_shape = Some(1.0);
        let result = run_nj(conf, None).unwrap();
        for i in 0..7 {
            assert!(
                result.newick.contains(&format!("s{i}:")),
                "{model:?}: {}",
                result.newick
            );
        }
        assert!(
            !result.newick.contains("NaN") && !result.newick.contains("inf"),
            "{model:?}"
        );
    }
}

#[test]
fn builtin_models_are_protein_only_and_take_no_matrix() {
    let dna = vec![seq("a", "ACGTACGTAC"), seq("b", "ACGTACGTAA")];
    for (model, text) in BUILTINS {
        let err = distance_matrix(dist_conf(dna.clone(), model.clone(), None)).unwrap_err();
        assert!(
            matches!(err, NJError::IncompatibleModel { .. }),
            "{model:?}: {err:?}"
        );
        let err = distance_matrix(dist_conf(family(3), model.clone(), Some(text))).unwrap_err();
        assert!(
            matches!(err, NJError::InvalidMatrix { .. }),
            "{model:?}: {err:?}"
        );
    }
}
