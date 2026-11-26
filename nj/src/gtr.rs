//! The general time-reversible (GTR) model for DNA: maximum-likelihood distances with explicit
//! exchangeabilities.
//!
//! GTR has six exchangeabilities `s_ab` and four base frequencies `π`; the rate from `a` to `b` is
//! `s_ab·π_b` (see [`crate::reversible`]). As in the design (`docs/design/substitution-models.md`
//! §5.1), the rates are **explicit**, with no built-in default, because any default would be
//! arbitrary: they are meant to come from an ML program that has fitted the model on a tree. The
//! frequencies are, unless given, pooled over the whole alignment (all unambiguous bases of all
//! sequences), which is unbiased under stationarity and far less noisy than a per-pair estimate.
//! Only the distance of each pair is fitted ([`crate::ml`]).
//!
//! Special cases are GTR with particular rates: all equal is F81 (JC with equal frequencies),
//! `[1, κ, 1, 1, κ, 1]` is HKY85 with ts/tv rate ratio `κ`, and so on.
//!
//! The order of the exchangeabilities is `AC, AG, AT, CG, CT, GT` and of the frequencies
//! `A, C, G, T`.

use crate::alphabet::{AlphabetEncoding, DNA};
use crate::ml::ml_distance;
use crate::models::{DistanceModel, RateHet};
use crate::msa::MSA;
use crate::prepared::Prepared;
use crate::reversible::ReversibleModel;

/// Number of exchangeabilities of a GTR model.
pub const N_RATES: usize = 6;

/// A GTR model with fixed exchangeabilities and frequencies, used as a [`DistanceModel`].
#[derive(Debug)]
pub struct Gtr {
    model: ReversibleModel,
}

impl Gtr {
    /// Builds the model from the six exchangeabilities (`AC, AG, AT, CG, CT, GT`, relative
    /// values, only their ratios matter) and the four frequencies (`A, C, G, T`, renormalised
    /// to sum to one).
    pub fn new(rates: &[f64], frequencies: &[f64]) -> Result<Self, String> {
        if rates.len() != N_RATES {
            return Err(format!(
                "GTR needs 6 exchangeabilities in the order AC, AG, AT, CG, CT, GT, got {}",
                rates.len()
            ));
        }
        if frequencies.len() != 4 {
            return Err(format!(
                "GTR needs 4 base frequencies in the order A, C, G, T, got {}",
                frequencies.len()
            ));
        }
        let mut s = vec![0.0; 16];
        let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for (&(a, b), &r) in pairs.iter().zip(rates) {
            s[a * 4 + b] = r;
            s[b * 4 + a] = r;
        }
        let model = ReversibleModel::new(&s, frequencies).map_err(|e| format!("GTR: {e}"))?;
        Ok(Self { model })
    }

    /// Like [`new`](Gtr::new), with the frequencies pooled from `msa` when `frequencies` is `None`.
    pub fn from_alignment(
        rates: &[f64],
        frequencies: Option<&[f64]>,
        msa: &MSA<DNA>,
    ) -> Result<Self, String> {
        match frequencies {
            Some(f) => Self::new(rates, f),
            None => Self::new(rates, &pooled_base_frequencies(msa)?),
        }
    }

    /// The underlying reversible model.
    pub fn model(&self) -> &ReversibleModel {
        &self.model
    }
}

impl DistanceModel<DNA> for Gtr {
    fn distance(&self, data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        ml_distance(&self.model, &data.dna_counts(i, j).to_f64(), rates)
    }
}

/// Base frequencies `A, C, G, T` pooled over every unambiguous base of every sequence (gaps and
/// `N` are not counted). Errors if a base never occurs, since a zero frequency makes the model
/// degenerate; supply the frequencies explicitly in that case.
pub fn pooled_base_frequencies(msa: &MSA<DNA>) -> Result<[f64; 4], String> {
    let mut counts = [0u64; 4];
    for seq in &msa.sequences {
        for &b in DNA::as_bytes(seq) {
            if b < 4 {
                counts[b as usize] += 1;
            }
        }
    }
    if let Some(missing) = counts.iter().position(|&c| c == 0) {
        return Err(format!(
            "GTR: base '{}' does not occur in the alignment, so its pooled frequency is zero; \
             supply base_frequencies",
            b"ACGT"[missing] as char
        ));
    }
    let total: u64 = counts.iter().sum();
    Ok(counts.map(|c| c as f64 / total as f64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counts::DnaCounts;
    use crate::models::{JukesCantor, pairwise_distance};
    use crate::reversible::testing::Rng;

    fn msa(seqs: &[&str]) -> MSA<DNA> {
        MSA::<DNA>::try_from_pairs(
            seqs.iter()
                .enumerate()
                .map(|(i, s)| (format!("s{i}"), s.to_string())),
        )
        .unwrap()
    }

    /// `n` sequences of `len` sites: a base sequence with composition skewed towards A and T,
    /// and sequence `i` substituted at about `5·i`% of its sites.
    fn skewed_alignment(n: usize, len: usize, seed: u64) -> MSA<DNA> {
        let mut rng = Rng(seed);
        let draw = |rng: &mut Rng| match rng.unit() {
            u if u < 0.4 => 'A',
            u if u < 0.55 => 'C',
            u if u < 0.7 => 'G',
            _ => 'T',
        };
        let base: Vec<char> = (0..len).map(|_| draw(&mut rng)).collect();
        let pairs = (0..n).map(|i| {
            let s: String = base
                .iter()
                .map(|&c| {
                    if rng.unit() < 0.05 * i as f64 {
                        draw(&mut rng)
                    } else {
                        c
                    }
                })
                .collect();
            (format!("s{i}"), s)
        });
        MSA::<DNA>::try_from_pairs(pairs).unwrap()
    }

    #[test]
    fn test_pooled_frequencies_count_only_unambiguous_bases() {
        let m = msa(&["AACG-N", "ACGTNN", "A-GTAC"]);
        let f = pooled_base_frequencies(&m).unwrap();
        // A: 2+1+2 = 5, C: 1+1+1 = 3, G: 1+1+1 = 3, T: 0+1+1 = 2  (13 bases)
        assert_eq!(f, [5.0 / 13.0, 3.0 / 13.0, 3.0 / 13.0, 2.0 / 13.0]);
        let missing = pooled_base_frequencies(&msa(&["AAAC", "AACC"])).unwrap_err();
        assert!(
            missing.contains("'G'") && missing.contains("base_frequencies"),
            "{missing}"
        );
    }

    #[test]
    fn test_equal_exchangeabilities_with_pooled_frequencies_are_f81() {
        // With all exchangeabilities equal the model is F81, whose transition probabilities are
        // explicit: with B = 1 − Σπ² and e = exp(−t/B),
        //   P_aa = e + π_a (1 − e),   P_ab = π_b (1 − e).
        // The ML distance is not the moment estimator −B·ln(1 − p/B) (it uses the per-base counts),
        // so check it against a direct maximisation of that likelihood by golden-section search.
        // This pins the use of the pooled frequencies, which the closed form of a pair does not.
        let (n, len) = (6, 600);
        let m = skewed_alignment(n, len, 11);
        let gtr = Gtr::from_alignment(&[1.0; 6], None, &m).unwrap();
        let pi = pooled_base_frequencies(&m).unwrap();
        assert!(pi[0] > 0.3 && pi[1] < 0.2, "composition is skewed: {pi:?}");
        let b = 1.0 - pi.iter().map(|p| p * p).sum::<f64>();
        let loglik = |f: &[f64; 16], t: f64| -> f64 {
            let e = (-t / b).exp();
            (0..16)
                .filter(|&k| f[k] > 0.0)
                .map(|k| {
                    let (a, c) = (k / 4, k % 4);
                    let p = if a == c {
                        e + pi[a] * (1.0 - e)
                    } else {
                        pi[c] * (1.0 - e)
                    };
                    f[k] * p.ln()
                })
                .sum()
        };
        let d = m.into_dist_with(&gtr, RateHet::NONE);
        let m = skewed_alignment(n, len, 11);
        for i in 0..n {
            for j in 0..i {
                let f = DnaCounts::from_symbols(&m.sequences[i], &m.sequences[j]).to_f64();
                let (mut lo, mut hi) = (1e-9, 5.0);
                let phi = (5f64.sqrt() - 1.0) / 2.0;
                for _ in 0..200 {
                    let (x1, x2) = (hi - phi * (hi - lo), lo + phi * (hi - lo));
                    if loglik(&f, x1) < loglik(&f, x2) {
                        lo = x1;
                    } else {
                        hi = x2;
                    }
                }
                let want = 0.5 * (lo + hi);
                assert!(
                    (d.get(i, j) - want).abs() < 1e-6,
                    "pair ({i},{j}): {} vs {want}",
                    d.get(i, j)
                );
            }
        }
    }

    #[test]
    fn test_equal_rates_and_frequencies_match_jukes_cantor() {
        let m = skewed_alignment(5, 400, 12);
        let gtr = Gtr::new(&[1.0; 6], &[0.25; 4]).unwrap();
        let d = m.into_dist_with(&gtr, RateHet::NONE);
        let m = skewed_alignment(5, 400, 12);
        for i in 0..5 {
            for j in 0..i {
                let want = pairwise_distance::<JukesCantor, DNA>(&m.sequences[i], &m.sequences[j]);
                assert!((d.get(i, j) - want).abs() < 1e-8, "pair ({i},{j})");
            }
        }
    }

    #[test]
    fn test_hky_at_the_pairs_own_kappa_reproduces_the_kimura_two_parameter_distance() {
        // Kimura's closed form estimates the ts/tv ratio of each pair as well as its distance.
        // GTR with equal frequencies and exchangeabilities [1, κ, 1, 1, κ, 1] at that same κ
        // must therefore return the same distance. This pins the order of the six rates:
        // AG and CT (positions 2 and 5) are the transitions.
        let cases: [(f64, f64, f64); 4] = [
            (700.0, 60.0, 40.0),
            (850.0, 90.0, 20.0),
            (500.0, 100.0, 150.0),
            (950.0, 10.0, 40.0),
        ];
        for (same, ts, tv) in cases {
            let n = same + ts + tv;
            let (p, q) = (ts / n, tv / n);
            let k2p = -0.5 * (1.0 - 2.0 * p - q).ln() - 0.25 * (1.0 - 2.0 * q).ln();
            let kappa = 2.0 * (1.0 - 2.0 * p - q).ln() / (1.0 - 2.0 * q).ln() - 1.0;
            let gtr = Gtr::new(&[1.0, kappa, 1.0, 1.0, kappa, 1.0], &[0.25; 4]).unwrap();
            // Counts: transitions split over A<->G and C<->T, transversions over the other pairs.
            let mut f = [0.0; 16];
            f[0] = same * 0.5;
            f[5] = same * 0.5;
            f[2] = ts * 0.25;
            f[8] = ts * 0.25;
            f[7] = ts * 0.25;
            f[13] = ts * 0.25;
            f[1] = tv * 0.25;
            f[4] = tv * 0.25;
            f[3] = tv * 0.25;
            f[12] = tv * 0.25;
            let got = ml_distance(gtr.model(), &f, RateHet::NONE);
            assert!(
                (got - k2p).abs() < 1e-8,
                "same {same}, ts {ts}, tv {tv}: GTR {got} vs K2P {k2p} (kappa {kappa})"
            );
        }
    }

    /// `exp(Q·t)` for a 4×4 matrix by scaling and squaring a Taylor series: deliberately unrelated
    /// to the eigendecomposition used by the model.
    fn expm(q: &[[f64; 4]; 4], t: f64) -> [[f64; 4]; 4] {
        let squarings = 12;
        let scale = t / (1u64 << squarings) as f64;
        let mut a = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                a[i][j] = q[i][j] * scale;
            }
        }
        let mul = |x: &[[f64; 4]; 4], y: &[[f64; 4]; 4]| {
            let mut r = [[0.0; 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    r[i][j] = (0..4).map(|k| x[i][k] * y[k][j]).sum();
                }
            }
            r
        };
        let mut result = [[0.0; 4]; 4];
        let mut term = [[0.0; 4]; 4];
        for i in 0..4 {
            result[i][i] = 1.0;
            term[i][i] = 1.0;
        }
        for k in 1..=18 {
            term = mul(&term, &a);
            for i in 0..4 {
                for j in 0..4 {
                    term[i][j] /= k as f64;
                    result[i][j] += term[i][j];
                }
            }
        }
        for _ in 0..squarings {
            result = mul(&result, &result);
        }
        result
    }

    #[test]
    fn test_asymmetric_rates_match_an_independent_matrix_exponential() {
        // Every exchangeability different, so a wrong pairing of rates and base pairs changes the
        // answer. Q is built from *named* pairs here, then exponentiated by an unrelated routine,
        // and the likelihood is maximised by golden-section search.
        let named = [
            ("AC", 0.5),
            ("AG", 3.0),
            ("AT", 0.8),
            ("CG", 1.3),
            ("CT", 6.0),
            ("GT", 0.2),
        ];
        let rates: Vec<f64> = named.iter().map(|&(_, r)| r).collect();
        let pi = [0.35, 0.15, 0.2, 0.3];
        let base = |c: char| "ACGT".find(c).unwrap();
        let mut q = [[0.0; 4]; 4];
        for &(pair, r) in &named {
            let (a, b) = (
                base(pair.chars().next().unwrap()),
                base(pair.chars().nth(1).unwrap()),
            );
            q[a][b] = r * pi[b];
            q[b][a] = r * pi[a];
        }
        let mut mean_rate = 0.0;
        for a in 0..4 {
            q[a][a] = -(0..4).filter(|&b| b != a).map(|b| q[a][b]).sum::<f64>();
            mean_rate -= pi[a] * q[a][a];
        }
        for row in q.iter_mut() {
            for v in row.iter_mut() {
                *v /= mean_rate;
            }
        }
        // An asymmetric table of aligned pairs (row = base in the first sequence).
        #[rustfmt::skip]
        let f: [f64; 16] = [
            80.0,  5.0, 12.0,  4.0,
             3.0, 30.0,  2.0, 15.0,
             7.0,  3.0, 40.0,  4.0,
             6.0, 10.0,  5.0, 60.0,
        ];
        let loglik = |t: f64| {
            let p = expm(&q, t);
            (0..16).map(|k| f[k] * p[k / 4][k % 4].ln()).sum::<f64>()
        };
        let (mut lo, mut hi) = (1e-6, 6.0);
        let phi = (5f64.sqrt() - 1.0) / 2.0;
        for _ in 0..200 {
            let (x1, x2) = (hi - phi * (hi - lo), lo + phi * (hi - lo));
            if loglik(x1) < loglik(x2) {
                lo = x1;
            } else {
                hi = x2;
            }
        }
        let want = 0.5 * (lo + hi);
        let gtr = Gtr::new(&rates, &pi).unwrap();
        let got = ml_distance(gtr.model(), &f, RateHet::NONE);
        assert!(
            (got - want).abs() < 1e-6,
            "model {got} vs independent {want}"
        );
        // Swapping two exchangeabilities changes the distance noticeably: the check can fail.
        let mut swapped = rates.clone();
        swapped.swap(1, 4);
        let other = ml_distance(Gtr::new(&swapped, &pi).unwrap().model(), &f, RateHet::NONE);
        assert!(
            (other - want).abs() > 1e-3,
            "swapping AG and CT must matter ({other} vs {want})"
        );
    }

    #[test]
    fn test_gaps_and_n_change_nothing_and_identical_sequences_are_zero() {
        let a = "ACGTACGTACGTACGGACGT";
        let b = "ACGTTCGTACGAACGGACGT";
        let gtr = Gtr::new(&[1.0, 4.0, 0.8, 1.2, 3.5, 1.0], &[0.3, 0.2, 0.2, 0.3]).unwrap();
        let dist = |x: &str, y: &str| {
            let m = msa(&[x, y]);
            m.into_dist_with(&gtr, RateHet::NONE).get(1, 0)
        };
        let base = dist(a, b);
        assert!(base > 0.0 && base.is_finite());
        for (xa, xb) in [("-", "A"), ("N", "C"), ("G", "N"), ("-", "-")] {
            assert_eq!(
                dist(&format!("{a}{xa}"), &format!("{b}{xb}")).to_bits(),
                base.to_bits(),
                "({xa},{xb})"
            );
        }
        assert_eq!(dist(a, a), 0.0);
    }

    #[test]
    fn test_invalid_parameters_are_rejected() {
        let ok = [1.0, 2.0, 1.0, 1.0, 2.0, 1.0];
        assert!(
            Gtr::new(&ok[..5], &[0.25; 4])
                .unwrap_err()
                .contains("6 exchangeabilities")
        );
        assert!(
            Gtr::new(&ok, &[0.5, 0.5, 0.0])
                .unwrap_err()
                .contains("4 base frequencies")
        );
        assert!(
            Gtr::new(&ok, &[0.5, 0.5, 0.0, 0.0]).is_err(),
            "zero frequency"
        );
        assert!(
            Gtr::new(&[1.0, -1.0, 1.0, 1.0, 1.0, 1.0], &[0.25; 4]).is_err(),
            "negative rate"
        );
        assert!(Gtr::new(&[0.0; 6], &[0.25; 4]).is_err(), "no substitutions");
        assert!(
            Gtr::new(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], &[0.25; 4]).is_err(),
            "only A<->C: reducible"
        );
        assert!(Gtr::new(&[1.0, f64::NAN, 1.0, 1.0, 1.0, 1.0], &[0.25; 4]).is_err());
    }
}
