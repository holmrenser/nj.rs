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
//!
//! [`EmpiricalGtr`] is the other GTR distance: it needs no rates, because it estimates the whole
//! model from each pair (see there).

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol};
use crate::counts::DnaCounts;
use crate::ml::{T_MAX, ml_distance};
use crate::models::{DistanceModel, ModelCalculation, RateHet};
use crate::msa::MSA;
use crate::prepared::Prepared;
use crate::reversible::{ReversibleModel, jacobi_eigen};

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

    /// The F84 model of Felsenstein (1984), with the `κ` of PAML's F84 (`model = 3` in `baseml`; checked
    /// against it in `tests/paml_reference.rs`).
    ///
    /// Transversions have exchangeability 1, and the two kinds of transition get an extra rate
    /// within their class: `1 + κ/π_R` for `A↔G` and `1 + κ/π_Y` for `C↔T`, where `π_R = π_A + π_G` and
    /// `π_Y = π_C + π_T`. With `κ = 0` this is F81. Unlike HKY85, F84 has two different transition
    /// rates unless `π_R = π_Y`; `κ` is not the transition/transversion rate ratio.
    pub fn f84(kappa: f64, frequencies: &[f64]) -> Result<Self, String> {
        if !(kappa.is_finite() && kappa >= 0.0) {
            return Err(format!(
                "F84: kappa must be finite and non-negative, got {kappa}"
            ));
        }
        if frequencies.len() != 4 {
            return Err(format!(
                "F84 needs 4 base frequencies in the order A, C, G, T, got {}",
                frequencies.len()
            ));
        }
        let total: f64 = frequencies.iter().sum();
        let pi: Vec<f64> = frequencies.iter().map(|f| f / total).collect();
        let (pi_r, pi_y) = (pi[0] + pi[2], pi[1] + pi[3]);
        if !(pi_r > 0.0 && pi_y > 0.0) {
            return Err("F84: base frequencies must be positive".into());
        }
        let rates = [1.0, 1.0 + kappa / pi_r, 1.0, 1.0, 1.0 + kappa / pi_y, 1.0];
        Self::new(&rates, &pi).map_err(|e| e.replacen("GTR", "F84", 1))
    }

    /// Like [`f84`](Gtr::f84), with the frequencies pooled from `msa` when `frequencies` is `None`.
    pub fn f84_from_alignment(
        kappa: f64,
        frequencies: Option<&[f64]>,
        msa: &MSA<DNA>,
    ) -> Result<Self, String> {
        match frequencies {
            Some(f) => Self::f84(kappa, f),
            None => Self::f84(kappa, &pooled_base_frequencies(msa)?),
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

/// The empirical GTR distance of Lanave et al. (1984), Tavaré (1986) and Rodríguez et al. (1990):
/// the closed-form distance of a general time-reversible model whose rates and frequencies are all
/// estimated from the pair itself.
///
/// With `J` the table of aligned base pairs as joint frequencies, symmetrised (`(J + Jᵀ)/2`, since
/// a reversible model has `J = Π P(t)` symmetric) and `π` its marginal, a GTR model predicts
/// `Π⁻¹J = exp(Qt)`, so `Qt = ln(Π⁻¹J)` and the distance (substitutions per site) is
///
/// ```text
/// d = −tr(Π ln(Π⁻¹J)) = −Σ_k w_k ln M_k ,   S = Π^-½ J Π^-½ = V M Vᵀ ,   w_k = Σ_i π_i V_ik² .
/// ```
///
/// `S` is symmetric, so this is one 4×4 eigenproblem. It is the maximum-likelihood estimate for the
/// model with a free symmetric `J` (invariance of ML estimates), so it needs no rates and no
/// numerical search; the price is that all parameters are estimated from one pair, which is
/// noisy for short sequences. If an eigenvalue of `S` is not positive, `ln` does not exist: the
/// pair is saturated (distance infinity), which for short alignments or large distances happens
/// to a sizeable share of pairs (about 38% at 100 sites and `t = 1.2`).
///
/// Among-site rate variation enters through the eigenvalues, as for the closed-form models:
/// gamma rates change `−ln M` to `α(M^{−1/α} − 1)`, and invariant sites change `M` to
/// `(M − p_inv)/(1 − p_inv)` and scale the result by `1 − p_inv`. A base that does not occur in the
/// pair has no frequency and is left out of the model. As in [`crate::ml`], a distance beyond
/// [`T_MAX`] substitutions per site counts as saturated.
#[derive(Clone, Copy, Debug, Default)]
pub struct EmpiricalGtr;

impl ModelCalculation<DNA> for EmpiricalGtr {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], rates: RateHet) -> f64 {
        empirical_gtr_distance(&DnaCounts::from_symbols(s1, s2).as_table(), rates)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        empirical_gtr_distance(&data.dna_counts(i, j).as_table(), rates)
    }
}

/// The empirical GTR distance of a table of aligned base pairs (`table[a][b]`: bases `a` and `b` in
/// the two sequences; counts or frequencies, only the proportions matter). See [`EmpiricalGtr`].
pub(crate) fn empirical_gtr_distance(table: &[[f64; 4]; 4], het: RateHet) -> f64 {
    let total: f64 = table.iter().flatten().sum();
    let identical: f64 = (0..4).map(|a| table[a][a]).sum();
    if total <= 0.0 || identical >= total {
        return 0.0;
    }

    // Symmetrised joint frequencies and their marginal.
    let mut joint = [[0.0; 4]; 4];
    let mut pi = [0.0; 4];
    for a in 0..4 {
        for b in 0..4 {
            joint[a][b] = (table[a][b] + table[b][a]) / (2.0 * total);
            pi[a] += joint[a][b];
        }
    }
    let kept: Vec<usize> = (0..4).filter(|&a| pi[a] > 0.0).collect();
    let k = kept.len();
    let mut s = vec![0.0; k * k];
    for (x, &a) in kept.iter().enumerate() {
        for (y, &b) in kept.iter().enumerate() {
            s[x * k + y] = joint[a][b] / (pi[a] * pi[b]).sqrt();
        }
    }
    let (eigenvalues, v) = jacobi_eigen(&mut s, k);

    let p_inv = het.p_invar;
    let mut d = 0.0;
    for (m, &lambda) in eigenvalues.iter().enumerate() {
        let lambda = (lambda - p_inv) / (1.0 - p_inv);
        if !(lambda > 1e-9) {
            return f64::INFINITY;
        }
        let weight: f64 = kept
            .iter()
            .enumerate()
            .map(|(x, &a)| pi[a] * v[x * k + m] * v[x * k + m])
            .sum();
        d += match het.gamma_shape {
            None => -weight * lambda.ln(),
            Some(alpha) => weight * alpha * (lambda.powf(-1.0 / alpha) - 1.0),
        };
    }
    if d > T_MAX {
        return f64::INFINITY;
    }
    (1.0 - p_inv) * d.max(0.0)
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

    /// The rate matrix (one substitution per unit time) of a GTR model given by *named* pairs, built
    /// by hand so that it does not depend on the model code under test.
    fn named_q(named: &[(&str, f64); 6], pi: &[f64; 4]) -> [[f64; 4]; 4] {
        let base = |c: char| "ACGT".find(c).unwrap();
        let mut q = [[0.0; 4]; 4];
        for &(pair, r) in named {
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
        q
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
        let q = named_q(&named, &pi);
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

    /// Expected pair table `π_a P_ab(t)` of a GTR model (no rate variation), by matrix exponential.
    fn expected_table(q: &[[f64; 4]; 4], pi: &[f64; 4], t: f64) -> [[f64; 4]; 4] {
        let p = expm(q, t);
        let mut f = [[0.0; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                f[a][b] = pi[a] * p[a][b];
            }
        }
        f
    }

    #[test]
    fn test_empirical_gtr_recovers_the_distance_of_a_known_model_from_its_expected_table() {
        // With the exact expected table of a GTR model (symmetric, J = Π P(t)) the closed form must
        // return t exactly, whatever the rates and frequencies; Q is scaled to one substitution per
        // unit time, so d = t.
        let models = [
            (
                [
                    ("AC", 0.5),
                    ("AG", 3.0),
                    ("AT", 0.8),
                    ("CG", 1.3),
                    ("CT", 6.0),
                    ("GT", 0.2),
                ],
                [0.35, 0.15, 0.2, 0.3],
            ),
            (
                [
                    ("AC", 1.0),
                    ("AG", 1.0),
                    ("AT", 1.0),
                    ("CG", 1.0),
                    ("CT", 1.0),
                    ("GT", 1.0),
                ],
                [0.25; 4],
            ),
            (
                [
                    ("AC", 1.0),
                    ("AG", 8.0),
                    ("AT", 1.0),
                    ("CG", 1.0),
                    ("CT", 8.0),
                    ("GT", 1.0),
                ],
                [0.4, 0.1, 0.1, 0.4],
            ),
            (
                [
                    ("AC", 2.2),
                    ("AG", 0.4),
                    ("AT", 5.0),
                    ("CG", 0.9),
                    ("CT", 1.7),
                    ("GT", 3.1),
                ],
                [0.1, 0.4, 0.3, 0.2],
            ),
        ];
        for (named, pi) in models {
            let q = named_q(&named, &pi);
            for &t in &[0.02, 0.2, 0.7, 1.5, 2.5] {
                let f = expected_table(&q, &pi, t);
                let d = empirical_gtr_distance(&f, RateHet::NONE);
                assert!(
                    (d - t).abs() < 1e-8 * t.max(1.0),
                    "{named:?}, t {t}: got {d}"
                );
            }
        }
    }

    #[test]
    fn test_empirical_gtr_with_rate_variation_matches_the_ml_solver_on_expected_tables() {
        // Gamma and invariant sites through the eigenvalues must agree with the ML solver, which
        // handles them through the transition probabilities, on the exact expected table.
        let mut rng = Rng(0x5EED_0040);
        for trial in 0..12 {
            let rates: Vec<f64> = (0..6).map(|_| 0.2 + 4.0 * rng.unit()).collect();
            let pi_raw: Vec<f64> = (0..4).map(|_| 0.15 + rng.unit()).collect();
            let gtr = Gtr::new(&rates, &pi_raw).unwrap();
            let pi = gtr.model().frequencies().to_vec();
            let het = [
                RateHet {
                    gamma_shape: Some(0.5),
                    p_invar: 0.0,
                },
                RateHet {
                    gamma_shape: Some(2.0),
                    p_invar: 0.0,
                },
                RateHet {
                    gamma_shape: None,
                    p_invar: 0.25,
                },
                RateHet {
                    gamma_shape: Some(0.8),
                    p_invar: 0.3,
                },
            ][trial % 4];
            let t = [0.1, 0.4, 0.9][trial % 3];
            let p = gtr.model().transition_matrix(t, het);
            let mut table = [[0.0; 4]; 4];
            let mut flat = [0.0; 16];
            for a in 0..4 {
                for b in 0..4 {
                    table[a][b] = pi[a] * p[a * 4 + b];
                    flat[a * 4 + b] = table[a][b] * 1e6;
                }
            }
            let empirical = empirical_gtr_distance(&table, het);
            let ml = ml_distance(gtr.model(), &flat, het);
            let want = (1.0 - het.p_invar) * t;
            assert!(
                (empirical - want).abs() < 1e-7,
                "trial {trial}, {het:?}: empirical {empirical} vs {want}"
            );
            assert!(
                (ml - want).abs() < 1e-6,
                "trial {trial}, {het:?}: ml {ml} vs {want}"
            );
        }
    }

    #[test]
    fn test_empirical_gtr_two_states_is_the_binary_symmetric_closed_form() {
        // Only A and C present, counts [[a, b], [b, a]]: π = ½, S has eigenvalues 1 and 1 − 2p, and
        // d = −½·ln(1 − 2p) with p = b/(a + b).
        let (a, b) = (170.0, 30.0);
        let mut f = [[0.0; 4]; 4];
        f[0][0] = a;
        f[1][1] = a;
        f[0][1] = b;
        f[1][0] = b;
        let p = 2.0 * b / (2.0 * a + 2.0 * b);
        let want = -0.5 * (1.0_f64 - 2.0 * p).ln();
        assert!((empirical_gtr_distance(&f, RateHet::NONE) - want).abs() < 1e-12);
    }

    #[test]
    fn test_empirical_gtr_saturated_tables_are_infinite() {
        // Two bases swapped in every column: eigenvalue −1.
        let mut swapped = [[0.0; 4]; 4];
        swapped[0][1] = 50.0;
        swapped[1][0] = 50.0;
        assert!(empirical_gtr_distance(&swapped, RateHet::NONE).is_infinite());
        // Independent sequences: J = π πᵀ, so every eigenvalue but one is 0.
        let pi = [0.4, 0.3, 0.2, 0.1];
        let mut independent = [[0.0; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                independent[a][b] = 1000.0 * pi[a] * pi[b];
            }
        }
        assert!(empirical_gtr_distance(&independent, RateHet::NONE).is_infinite());
        // A far but resolvable pair becomes saturated once half the sites are invariant, because the
        // variable sites then have to be explained as far more diverged.
        let q = named_q(
            &[
                ("AC", 1.0),
                ("AG", 1.0),
                ("AT", 1.0),
                ("CG", 1.0),
                ("CT", 1.0),
                ("GT", 1.0),
            ],
            &[0.25; 4],
        );
        let far = expected_table(&q, &[0.25; 4], 3.0);
        assert!(empirical_gtr_distance(&far, RateHet::NONE).is_finite());
        assert!(
            empirical_gtr_distance(
                &far,
                RateHet {
                    gamma_shape: None,
                    p_invar: 0.5
                }
            )
            .is_infinite()
        );
    }

    #[test]
    fn test_empirical_gtr_identical_and_empty_pairs_are_zero_and_scale_does_not_matter() {
        let mut same = [[0.0; 4]; 4];
        same[0][0] = 10.0;
        same[2][2] = 30.0;
        assert_eq!(empirical_gtr_distance(&same, RateHet::NONE), 0.0);
        assert_eq!(empirical_gtr_distance(&[[0.0; 4]; 4], RateHet::NONE), 0.0);
        let pi = [0.3, 0.2, 0.2, 0.3];
        let q = named_q(
            &[
                ("AC", 1.0),
                ("AG", 4.0),
                ("AT", 1.0),
                ("CG", 1.0),
                ("CT", 4.0),
                ("GT", 1.0),
            ],
            &pi,
        );
        let f = expected_table(&q, &pi, 0.6);
        let mut scaled = f;
        for row in scaled.iter_mut() {
            for v in row.iter_mut() {
                *v *= 12345.0;
            }
        }
        let (a, b) = (
            empirical_gtr_distance(&f, RateHet::NONE),
            empirical_gtr_distance(&scaled, RateHet::NONE),
        );
        assert!((a - b).abs() < 1e-12 && (a - 0.6).abs() < 1e-9);
    }

    #[test]
    fn test_empirical_gtr_through_the_distance_matrix_ignores_gaps_and_n_and_agrees_with_ml_on_long_data()
     {
        // Long, low-divergence alignment: the 9 parameters estimated per pair are well determined,
        // so the empirical distance is close to the ML distance under the generating model (F81 with
        // the skewed composition, i.e. equal exchangeabilities and pooled frequencies).
        let (n, len) = (5, 40_000);
        let m = skewed_alignment(n, len, 21);
        let ml_model = Gtr::from_alignment(&[1.0; 6], None, &m).unwrap();
        let ml = m.into_dist_with(&ml_model, RateHet::NONE);
        let m = skewed_alignment(n, len, 21);
        let emp = m.into_dist::<EmpiricalGtr>(RateHet::NONE);
        for i in 1..n {
            for j in 0..i {
                let (e, l) = (emp.get(i, j), ml.get(i, j));
                assert!(e.is_finite() && l > 0.0, "pair ({i},{j})");
                assert!(
                    (e - l).abs() < 0.03 * l + 2e-3,
                    "pair ({i},{j}): empirical {e} vs ML {l}"
                );
            }
        }
        // Gaps and N change nothing.
        let (a, b) = (
            "ACGTACGTAAGTACCTTCGATTAGCAAGT",
            "ACGTTCGTACGTACCTTAGATTAGCTAGT",
        );
        let dist = |x: &str, y: &str| {
            msa(&[x, y])
                .into_dist::<EmpiricalGtr>(RateHet::NONE)
                .get(1, 0)
        };
        let base = dist(a, b);
        assert!(base.is_finite() && base > 0.0, "{base}");
        for (xa, xb) in [("-", "A"), ("N", "C"), ("G", "N")] {
            assert_eq!(
                dist(&format!("{a}{xa}"), &format!("{b}{xb}")).to_bits(),
                base.to_bits()
            );
        }
    }

    /// Share of simulated pairs that the empirical GTR distance declares saturated, by length and
    /// distance. Run with `cargo test --release -p nj --lib empirical_gtr_failure -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn empirical_gtr_failure_rate_by_length_and_distance() {
        let mut rng = Rng(0x5EED_0041);
        let named = [
            ("AC", 1.0),
            ("AG", 4.0),
            ("AT", 1.0),
            ("CG", 1.0),
            ("CT", 4.0),
            ("GT", 1.0),
        ];
        let pi = [0.3, 0.2, 0.2, 0.3];
        let q = named_q(&named, &pi);
        println!("sites   t    saturated   median d    (ML solver saturated)");
        let gtr = Gtr::new(&named.map(|(_, r)| r), &pi).unwrap();
        for &len in &[100usize, 300, 1000, 5000] {
            for &t in &[0.1, 0.3, 0.6, 1.0, 1.5, 2.5] {
                let p = expm(&q, t);
                let trials = 400;
                let (mut sat, mut ml_sat) = (0, 0);
                let mut ds = Vec::new();
                for _ in 0..trials {
                    let mut f = [[0.0; 4]; 4];
                    for _ in 0..len {
                        let a = pick(&mut rng, &pi);
                        let b = pick(&mut rng, &p[a]);
                        f[a][b] += 1.0;
                    }
                    let d = empirical_gtr_distance(&f, RateHet::NONE);
                    if d.is_infinite() {
                        sat += 1;
                    } else {
                        ds.push(d);
                    }
                    let flat: Vec<f64> = f.iter().flatten().copied().collect();
                    ml_sat += ml_distance(gtr.model(), &flat, RateHet::NONE).is_infinite() as usize;
                }
                ds.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let median = ds.get(ds.len() / 2).copied().unwrap_or(f64::NAN);
                println!(
                    "{len:5} {t:4.1}   {:5.1}%      {median:6.3}     ({:.1}%)",
                    100.0 * sat as f64 / trials as f64,
                    100.0 * ml_sat as f64 / trials as f64
                );
            }
        }
    }

    fn pick(rng: &mut Rng, weights: &[f64]) -> usize {
        let u = rng.unit();
        let mut acc = 0.0;
        for (i, w) in weights.iter().enumerate() {
            acc += w;
            if u < acc {
                return i;
            }
        }
        weights.len() - 1
    }

    #[test]
    fn test_f84_rates_follow_felsenstein_with_the_kishino_hasegawa_parametrisation() {
        let (kappa, pi) = (1.7, [0.1, 0.4, 0.3, 0.2]);
        let (pi_r, pi_y) = (pi[0] + pi[2], pi[1] + pi[3]);
        let f84 = Gtr::f84(kappa, &pi).unwrap();
        let by_hand = Gtr::new(
            &[1.0, 1.0 + kappa / pi_r, 1.0, 1.0, 1.0 + kappa / pi_y, 1.0],
            &pi,
        )
        .unwrap();
        for &t in &[0.05, 0.6, 2.0] {
            let (a, b) = (
                f84.model().transition_matrix(t, RateHet::NONE),
                by_hand.model().transition_matrix(t, RateHet::NONE),
            );
            assert_eq!(a, b);
        }
        // With equal frequencies π_R = π_Y = 1/2, so both transitions get 1 + 2κ: HKY85 with that ratio.
        let equal = Gtr::f84(0.8, &[0.25; 4]).unwrap();
        let hky = Gtr::new(&[1.0, 2.6, 1.0, 1.0, 2.6, 1.0], &[0.25; 4]).unwrap();
        assert_eq!(
            equal.model().transition_matrix(0.9, RateHet::NONE),
            hky.model().transition_matrix(0.9, RateHet::NONE)
        );
        // kappa = 0 is F81: all exchangeabilities equal.
        let f81 = Gtr::new(&[1.0; 6], &pi).unwrap();
        assert_eq!(
            Gtr::f84(0.0, &pi)
                .unwrap()
                .model()
                .transition_matrix(0.4, RateHet::NONE),
            f81.model().transition_matrix(0.4, RateHet::NONE)
        );
        // Frequencies are renormalised before π_R and π_Y are formed.
        let scaled = Gtr::f84(kappa, &[1.0, 4.0, 3.0, 2.0]).unwrap();
        assert_eq!(
            scaled.model().transition_matrix(0.7, RateHet::NONE),
            f84.model().transition_matrix(0.7, RateHet::NONE)
        );
    }

    #[test]
    fn test_f84_invalid_parameters_are_rejected() {
        let pi = [0.25; 4];
        assert!(Gtr::f84(-0.1, &pi).unwrap_err().contains("kappa"));
        assert!(Gtr::f84(f64::NAN, &pi).is_err());
        assert!(Gtr::f84(f64::INFINITY, &pi).is_err());
        assert!(
            Gtr::f84(1.0, &[0.5, 0.5, 0.0])
                .unwrap_err()
                .contains("4 base frequencies")
        );
        assert!(
            Gtr::f84(1.0, &[0.5, 0.5, 0.0, 0.0]).is_err(),
            "zero frequencies"
        );
        assert!(Gtr::f84(1.0, &[0.0, 0.5, 0.0, 0.5]).is_err(), "no purines");
    }

    #[test]
    fn test_f84_pools_frequencies_from_the_alignment() {
        let m = skewed_alignment(4, 500, 31);
        let pooled = pooled_base_frequencies(&m).unwrap();
        let a = Gtr::f84_from_alignment(1.2, None, &m).unwrap();
        let b = Gtr::f84(1.2, &pooled).unwrap();
        assert_eq!(
            a.model().transition_matrix(0.5, RateHet::NONE),
            b.model().transition_matrix(0.5, RateHet::NONE)
        );
    }
}
