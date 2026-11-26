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
#[allow(unused_macros)]
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
#[allow(unused_imports)]
pub(crate) use protein_matrix;

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
