//! Empirical protein models: maximum-likelihood distances under a fixed rate matrix.
//!
//! A matrix such as WAG or LG is a table of exchangeabilities and frequencies, fixed once
//! and for all; only the distance `t` of each pair is fitted ([`crate::ml`]). A matrix is a
//! zero-sized marker type implementing [`ProteinMatrix`], whose [`ReversibleModel`] (with its
//! eigensystem) is built on first use and shared by every thread. [`ProteinMl`] turns a matrix
//! into a [`ModelCalculation`], so it plugs into the distance-matrix, NJ and bootstrap code
//! like every other model.
//!
//! Adding a matrix is one data file in PAML layout (strict lower triangle of the
//! exchangeabilities, then the 20 frequencies, in the order `ARNDCQEGHILKMFPSTWYV`) and one
//! `protein_matrix!` line.

use std::marker::PhantomData;

use crate::alphabet::Protein;
use crate::counts::ProteinCounts;
use crate::ml::ml_distance;
use crate::models::{DistanceModel, ModelCalculation, RateHet};
use crate::prepared::Prepared;
use crate::reversible::ReversibleModel;

/// A fixed protein rate matrix.
pub trait ProteinMatrix: 'static {
    /// The model, built once.
    fn model() -> &'static ReversibleModel;
}

/// The ML distance under matrix `M`, as a [`ModelCalculation`] for protein alignments.
///
/// Returns `0.0` for identical (or non-overlapping) sequences and `f64::INFINITY` for a
/// saturated pair. Gamma rate variation and invariant sites enter the likelihood through
/// `rates` (see [`ReversibleModel::transition_matrix`]).
pub struct ProteinMl<M>(PhantomData<fn() -> M>);

impl<M: ProteinMatrix> ModelCalculation<Protein> for ProteinMl<M> {
    fn distance(
        s1: &[<Protein as crate::alphabet::AlphabetEncoding>::Symbol],
        s2: &[<Protein as crate::alphabet::AlphabetEncoding>::Symbol],
        rates: RateHet,
    ) -> f64 {
        let counts = ProteinCounts::from_symbols(s1, s2);
        ml_distance(M::model(), &counts.to_f64(), rates)
    }

    fn distance_prepared(data: &Prepared<'_, Protein>, i: usize, j: usize, rates: RateHet) -> f64 {
        ml_distance(M::model(), &data.protein_counts(i, j).to_f64(), rates)
    }
}

/// A protein rate matrix supplied at run time (the `Custom` model): the same ML distance as
/// [`ProteinMl`], with the model held as a value instead of a type.
pub struct CustomProtein {
    model: ReversibleModel,
}

impl CustomProtein {
    /// Parses a matrix in PAML layout (see [`ReversibleModel::from_paml`]) for the 20 amino acids.
    pub fn from_paml(text: &str) -> Result<Self, String> {
        Ok(Self {
            model: ReversibleModel::from_paml(text, 20)?,
        })
    }

    /// Wraps an existing 20-state model.
    pub fn new(model: ReversibleModel) -> Result<Self, String> {
        if model.n_states() != 20 {
            return Err(format!(
                "a protein model needs 20 states, got {}",
                model.n_states()
            ));
        }
        Ok(Self { model })
    }
}

impl DistanceModel<Protein> for CustomProtein {
    fn distance(&self, data: &Prepared<'_, Protein>, i: usize, j: usize, rates: RateHet) -> f64 {
        ml_distance(&self.model, &data.protein_counts(i, j).to_f64(), rates)
    }
}

/// Defines a zero-sized [`ProteinMatrix`] from PAML-layout text, parsed on first use.
macro_rules! protein_matrix {
    ($(#[$meta:meta])* $name:ident, $text:expr) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug)]
        pub struct $name;

        impl $crate::empirical::ProteinMatrix for $name {
            fn model() -> &'static $crate::reversible::ReversibleModel {
                static MODEL: ::std::sync::OnceLock<$crate::reversible::ReversibleModel> = ::std::sync::OnceLock::new();
                MODEL.get_or_init(|| {
                    $crate::reversible::ReversibleModel::from_paml($text, 20)
                        .expect(concat!("embedded matrix ", stringify!($name), " is invalid"))
                })
            }
        }
    };
}

// The built-in matrices. Provenance, citations and checksums: `nj/data/matrices/README.md`.
protein_matrix!(
    /// LG: Le SQ, Gascuel O (2008) *Mol Biol Evol* 25:1307–1320, doi:10.1093/molbev/msn067.
    Lg,
    include_str!("../data/matrices/LG.dat")
);
protein_matrix!(
    /// WAG: Whelan S, Goldman N (2001) *Mol Biol Evol* 18:691–699, doi:10.1093/oxfordjournals.molbev.a003851.
    Wag,
    include_str!("../data/matrices/WAG.dat")
);
protein_matrix!(
    /// JTT: Jones DT, Taylor WR, Thornton JM (1992) *Comput Appl Biosci* 8:275–282, doi:10.1093/bioinformatics/8.3.275.
    Jtt,
    include_str!("../data/matrices/JTT.dat")
);
protein_matrix!(
    /// JTT recomputed by Kosiol C, Goldman N (2005) *Mol Biol Evol* 22:193–199, doi:10.1093/molbev/msi005.
    JttDcmut,
    include_str!("../data/matrices/DCMut-JTT.dat")
);
protein_matrix!(
    /// Dayhoff MO, Schwartz RM, Orcutt BC (1978) *Atlas of Protein Sequence and Structure* 5:345–352.
    Dayhoff,
    include_str!("../data/matrices/Dayhoff.dat")
);
protein_matrix!(
    /// Dayhoff recomputed by Kosiol C, Goldman N (2005) *Mol Biol Evol* 22:193–199, doi:10.1093/molbev/msi005.
    DayhoffDcmut,
    include_str!("../data/matrices/DCMut-Dayhoff.dat")
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Protein;
    use crate::distance_matrix::DistMat;
    use crate::msa::MSA;
    use crate::reversible::testing::{Rng, random_model};
    use std::sync::OnceLock;

    /// A random 20-state model standing in for a published matrix.
    struct Synthetic;
    impl ProteinMatrix for Synthetic {
        fn model() -> &'static ReversibleModel {
            static MODEL: OnceLock<ReversibleModel> = OnceLock::new();
            MODEL.get_or_init(|| random_model(&mut Rng(0x5EED_0020), 20))
        }
    }
    type Ml = ProteinMl<Synthetic>;

    /// `n` sequences of length `len`: a random base sequence with progressively more
    /// substitutions, a gap and an `X` here and there.
    fn alignment(n: usize, len: usize, seed: u64) -> MSA<Protein> {
        let mut rng = Rng(seed);
        let letters = b"ARNDCQEGHILKMFPSTWYV";
        let base: Vec<u8> = (0..len)
            .map(|_| letters[(rng.next() % 20) as usize])
            .collect();
        let pairs = (0..n).map(|i| {
            let mut s = base.clone();
            for c in s.iter_mut() {
                if rng.unit() < 0.04 * i as f64 {
                    *c = letters[(rng.next() % 20) as usize];
                }
            }
            if i % 3 == 1 {
                s[len / 2] = b'-';
                s[len / 3] = b'X';
            }
            (format!("s{i}"), String::from_utf8(s).unwrap())
        });
        MSA::<Protein>::try_from_pairs(pairs).unwrap()
    }

    #[test]
    fn test_prepared_path_matches_the_symbol_path_and_the_solver() {
        let msa = alignment(9, 400, 1);
        let prepared = Prepared::new(&msa);
        for i in 0..9 {
            for j in 0..9 {
                let by_symbols = Ml::distance(&msa.sequences[i], &msa.sequences[j], RateHet::NONE);
                let by_prepared = Ml::distance_prepared(&prepared, i, j, RateHet::NONE);
                assert_eq!(
                    by_symbols.to_bits(),
                    by_prepared.to_bits(),
                    "pair ({i},{j})"
                );
                let counts = ProteinCounts::from_symbols(&msa.sequences[i], &msa.sequences[j]);
                let direct = ml_distance(Synthetic::model(), &counts.to_f64(), RateHet::NONE);
                assert_eq!(by_symbols.to_bits(), direct.to_bits(), "pair ({i},{j})");
            }
        }
    }

    #[test]
    fn test_distance_matrix_is_symmetric_finite_and_orders_by_divergence() {
        let msa = alignment(8, 500, 2);
        let d: DistMat = msa.into_dist::<Ml>(RateHet::NONE);
        let n = 8;
        for i in 0..n {
            assert_eq!(d.get(i, i), 0.0);
            for j in 0..n {
                assert_eq!(d.get(i, j).to_bits(), d.get(j, i).to_bits());
                assert!(d.get(i, j).is_finite() && d.get(i, j) >= 0.0);
            }
        }
        // Sequence 0 is the ancestor of the others; later ones diverge further from it.
        for i in 2..n {
            assert!(
                d.get(0, i) > d.get(0, i - 1) * 0.8,
                "divergence should grow with i"
            );
        }
        assert!(d.get(0, n - 1) > d.get(0, 1));
    }

    #[test]
    fn test_gaps_and_x_in_either_sequence_change_nothing() {
        let a = "ARNDCQEGHILKMFPSTWYVARNDCQEGHIL";
        let b = "ARNDCQEGHILKMFPSTWYAARNDCQEGHIK";
        let msa = |x: &str, y: &str| {
            MSA::<Protein>::try_from_pairs(
                [
                    ("a".to_string(), x.to_string()),
                    ("b".to_string(), y.to_string()),
                ]
                .into_iter(),
            )
            .unwrap()
        };
        let base = msa(a, b);
        let want = Ml::distance(&base.sequences[0], &base.sequences[1], RateHet::NONE);
        assert!(want > 0.0 && want.is_finite());
        for (extra_a, extra_b) in [("-", "W"), ("W", "-"), ("X", "L"), ("L", "X"), ("-", "X")] {
            let m = msa(&format!("{a}{extra_a}"), &format!("{b}{extra_b}"));
            let got = Ml::distance(&m.sequences[0], &m.sequences[1], RateHet::NONE);
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "extra column ({extra_a}, {extra_b})"
            );
        }
    }

    #[test]
    fn test_full_nj_run_with_bootstrap_and_rate_variation() {
        let msa = alignment(10, 300, 3);
        let rates = RateHet {
            gamma_shape: Some(0.7),
            p_invar: 0.1,
        };
        let result = crate::runner::run_nj(
            msa,
            &crate::models::Static::<Ml>::new(),
            20,
            None,
            None,
            true,
            true,
            rates,
            10.0,
        )
        .unwrap();
        assert!(result.newick.ends_with(';'));
        assert!(!result.newick.contains("NaN") && !result.newick.contains("inf"));
        assert!(result.average_distance.unwrap() > 0.0);
        // All ten taxa appear in the tree.
        for i in 0..10 {
            assert!(
                result.newick.contains(&format!("s{i}:")),
                "taxon s{i} missing: {}",
                result.newick
            );
        }
    }

    fn builtin_matrices() -> Vec<(&'static str, &'static ReversibleModel, &'static str)> {
        vec![
            ("LG", Lg::model(), include_str!("../data/matrices/LG.dat")),
            (
                "WAG",
                Wag::model(),
                include_str!("../data/matrices/WAG.dat"),
            ),
            (
                "JTT",
                Jtt::model(),
                include_str!("../data/matrices/JTT.dat"),
            ),
            (
                "JTT-DCMut",
                JttDcmut::model(),
                include_str!("../data/matrices/DCMut-JTT.dat"),
            ),
            (
                "Dayhoff",
                Dayhoff::model(),
                include_str!("../data/matrices/Dayhoff.dat"),
            ),
            (
                "DCMut",
                DayhoffDcmut::model(),
                include_str!("../data/matrices/DCMut-Dayhoff.dat"),
            ),
        ]
    }

    #[test]
    fn test_builtin_matrices_are_valid_reversible_models() {
        for (name, model, text) in builtin_matrices() {
            assert_eq!(model.n_states(), 20, "{name}");
            let pi = model.frequencies();
            assert!(
                (pi.iter().sum::<f64>() - 1.0).abs() < 1e-12,
                "{name}: frequencies sum to one"
            );
            assert!(
                pi.iter().all(|&f| f > 0.005 && f < 0.2),
                "{name}: frequencies are plausible"
            );

            // The frequencies are the 20 numbers after the 190 exchangeabilities in the file
            // (renormalised, since the files give six digits).
            let numbers: Vec<f64> = text
                .split_whitespace()
                .take(210)
                .map(|t| t.parse().unwrap())
                .collect();
            let file_pi = &numbers[190..];
            let total: f64 = file_pi.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-3,
                "{name}: file frequencies sum to {total}"
            );
            for a in 0..20 {
                assert!(
                    (pi[a] - file_pi[a] / total).abs() < 1e-12,
                    "{name}: frequency {a}"
                );
            }

            for &t in &[0.02, 0.5, 3.0] {
                let p = model.transition_matrix(t, RateHet::NONE);
                for a in 0..20 {
                    let row: f64 = (0..20).map(|b| p[a * 20 + b]).sum();
                    assert!((row - 1.0).abs() < 1e-11, "{name}: row {a} at t={t}");
                    for b in 0..20 {
                        assert!(p[a * 20 + b] > -1e-13, "{name}: P non-negative");
                        let (l, r) = (pi[a] * p[a * 20 + b], pi[b] * p[b * 20 + a]);
                        assert!(
                            (l - r).abs() < 1e-12,
                            "{name}: detailed balance ({a},{b}) at t={t}"
                        );
                    }
                }
            }
            // One substitution per site per unit time.
            let t = 1e-6;
            let p = model.transition_matrix(t, RateHet::NONE);
            let moved: f64 = (0..20).map(|a| pi[a] * (1.0 - p[a * 20 + a])).sum();
            assert!(
                (moved / t - 1.0).abs() < 1e-4,
                "{name}: mean rate {}",
                moved / t
            );
            // All eigenvalues are negative apart from the stationary one.
            assert_eq!(
                model.eigenvalues().iter().filter(|&&l| l == 0.0).count(),
                1,
                "{name}"
            );
        }
    }

    #[test]
    fn test_builtin_matrices_use_the_amino_acid_order_of_the_alphabet() {
        // Exchangeabilities are proportional to P_ab(t) / π_b for small t. Conservative
        // replacements (similar size and chemistry) are the most exchangeable pairs in every
        // empirical matrix and clearly dissimilar pairs are among the least. A matrix in another
        // residue order than `ARNDCQEGHILKMFPSTWYV` would fail this with probability ≈ 1 − 10⁻⁷.
        let idx = |c: char| {
            b"ARNDCQEGHILKMFPSTWYV"
                .iter()
                .position(|&x| x as char == c)
                .unwrap()
        };
        let conservative = [
            "IV", "IL", "LM", "IM", "DE", "ND", "QE", "KR", "FY", "ST", "AS", "FL",
        ];
        let dissimilar = ["WD", "CD", "PW"];
        for (name, model, _) in builtin_matrices() {
            let p = model.transition_matrix(1e-4, RateHet::NONE);
            let pi = model.frequencies();
            let s = |a: usize, b: usize| p[a * 20 + b] / pi[b];
            let mut all: Vec<f64> = (0..20)
                .flat_map(|a| ((a + 1)..20).map(move |b| (a, b)))
                .map(|(a, b)| s(a, b))
                .collect();
            all.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let (median, q75) = (all[all.len() / 2], all[all.len() * 3 / 4]);
            for pair in conservative {
                let (a, b) = (
                    idx(pair.chars().next().unwrap()),
                    idx(pair.chars().nth(1).unwrap()),
                );
                assert!(
                    s(a, b) > q75,
                    "{name}: {pair} should be among the most exchangeable pairs"
                );
            }
            for pair in dissimilar {
                let (a, b) = (
                    idx(pair.chars().next().unwrap()),
                    idx(pair.chars().nth(1).unwrap()),
                );
                assert!(
                    s(a, b) < median,
                    "{name}: {pair} should be below the median"
                );
            }
            // W and C are among the rarest residues, L and A among the commonest.
            assert!(
                pi[idx('W')] < 0.02 && pi[idx('C')] < 0.04,
                "{name}: W/C frequencies"
            );
            assert!(
                pi[idx('L')] > 0.07 && pi[idx('A')] > 0.06,
                "{name}: L/A frequencies"
            );
        }
    }

    #[test]
    fn test_builtin_matrices_differ_from_each_other() {
        // Guard against two models silently loading the same file.
        let all = builtin_matrices();
        for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                let (a, b) = (
                    all[i].1.transition_matrix(0.5, RateHet::NONE),
                    all[j].1.transition_matrix(0.5, RateHet::NONE),
                );
                let diff: f64 = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).sum();
                assert!(
                    diff > 1e-3,
                    "{} and {} are indistinguishable",
                    all[i].0,
                    all[j].0
                );
            }
        }
    }

    #[test]
    fn test_builtin_distances_are_finite_ordered_and_matrix_dependent() {
        let msa = alignment(6, 400, 4);
        let prepared = Prepared::new(&msa);
        let mut seen = Vec::new();
        for (name, _, _) in builtin_matrices() {
            let d = match name {
                "LG" => ProteinMl::<Lg>::distance_prepared(&prepared, 0, 5, RateHet::NONE),
                "WAG" => ProteinMl::<Wag>::distance_prepared(&prepared, 0, 5, RateHet::NONE),
                "JTT" => ProteinMl::<Jtt>::distance_prepared(&prepared, 0, 5, RateHet::NONE),
                "JTT-DCMut" => {
                    ProteinMl::<JttDcmut>::distance_prepared(&prepared, 0, 5, RateHet::NONE)
                }
                "Dayhoff" => {
                    ProteinMl::<Dayhoff>::distance_prepared(&prepared, 0, 5, RateHet::NONE)
                }
                _ => ProteinMl::<DayhoffDcmut>::distance_prepared(&prepared, 0, 5, RateHet::NONE),
            };
            assert!(d.is_finite() && d > 0.0, "{name}: {d}");
            seen.push(d);
        }
        // The matrices give different (but comparable) distances for the same pair.
        let (lo, hi) = seen
            .iter()
            .fold((f64::MAX, 0.0f64), |(l, h), &d| (l.min(d), h.max(d)));
        assert!(hi > lo * 1.0001 && hi < lo * 1.5, "distances {seen:?}");
    }

    #[test]
    fn test_protein_matrix_macro_parses_paml_text() {
        // 190 exchangeabilities of 1.0 and 20 equal frequencies is the equal-rates model.
        const TEXT: &str = include_str!("../tests/data/equal_rates_protein.dat");
        protein_matrix!(
            /// Test matrix defined by the macro.
            EqualRates,
            TEXT
        );
        let model = EqualRates::model();
        assert_eq!(model.n_states(), 20);
        let reference = ReversibleModel::equal_rates(20);
        let (p, q) = (
            model.transition_matrix(0.6, RateHet::NONE),
            reference.transition_matrix(0.6, RateHet::NONE),
        );
        for (x, y) in p.iter().zip(&q) {
            assert!((x - y).abs() < 1e-13);
        }
        // Same instance on every call.
        assert!(std::ptr::eq(EqualRates::model(), model));
    }
}
