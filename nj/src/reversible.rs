//! Time-reversible substitution models: rate matrix, eigensystem and transition
//! probabilities.
//!
//! A reversible model is a stationary distribution `π` and symmetric exchangeabilities
//! `s_ab`, with rates `Q_ab = s_ab·π_b` for `a ≠ b`, `Q_aa = −Σ_b Q_ab`, scaled to one
//! expected substitution per unit time (`−Σ_a π_a Q_aa = 1`), so a time `t` is a distance
//! in substitutions per site.
//!
//! Reversibility makes `A = Π^½ Q Π^-½` symmetric (`A_ab = √(π_a π_b)·s_ab`), so one
//! symmetric eigendecomposition `A = R Λ Rᵀ` gives every transition probability:
//!
//! ```text
//! P_ab(t) = Σ_m  R_am R_bm √(π_b/π_a) · e^{λ_m t}  =  Σ_m X_ab[m] · e^{λ_m t}
//! ```
//!
//! The coefficients `X_ab[m]` are computed once per model; evaluating `P`, `P′` and `P″`
//! for one entry is then three length-`K` dot products. Among-site rate variation only
//! changes the factor `e^{λt}` (see `ReversibleModel::factors`), so it costs nothing
//! extra per entry. Background and references: `docs/theory/distance-estimation.md`.

use crate::models::RateHet;

/// A time-reversible substitution model on `K` states with a cached eigensystem.
#[derive(Clone, Debug)]
pub struct ReversibleModel {
    k: usize,
    pi: Vec<f64>,
    /// Eigenvalues of the rate matrix: one exactly `0`, the rest negative.
    lambda: Vec<f64>,
    /// `x[(a·K + b)·K + m] = R_am R_bm √(π_b/π_a)`.
    x: Vec<f64>,
}

impl ReversibleModel {
    /// Builds a model from a symmetric exchangeability matrix (row-major `K×K`, the
    /// diagonal is ignored) and stationary frequencies, which are renormalised to sum
    /// to one. Rates are scaled to a mean of one substitution per unit time.
    pub fn new(exchangeabilities: &[f64], frequencies: &[f64]) -> Result<Self, String> {
        let k = frequencies.len();
        if k < 2 {
            return Err("a substitution model needs at least two states".into());
        }
        if exchangeabilities.len() != k * k {
            return Err(format!(
                "expected a {k}×{k} exchangeability matrix ({} values), got {}",
                k * k,
                exchangeabilities.len()
            ));
        }
        if frequencies.iter().any(|&f| !(f.is_finite() && f > 0.0)) {
            return Err("state frequencies must be finite and positive".into());
        }
        let total: f64 = frequencies.iter().sum();
        let pi: Vec<f64> = frequencies.iter().map(|f| f / total).collect();

        let mut s = vec![0.0; k * k];
        for a in 0..k {
            for b in 0..k {
                if a == b {
                    continue;
                }
                let (u, v) = (exchangeabilities[a * k + b], exchangeabilities[b * k + a]);
                if !(u.is_finite() && u >= 0.0) {
                    return Err("exchangeabilities must be finite and non-negative".into());
                }
                if (u - v).abs() > 1e-9 * u.abs().max(v.abs()).max(1e-300) {
                    return Err(format!("exchangeabilities are not symmetric at ({a}, {b})"));
                }
                s[a * k + b] = 0.5 * (u + v);
            }
        }

        // Row sums of Q give the diagonal; the mean rate fixes the time scale.
        let mut diag = vec![0.0; k];
        for a in 0..k {
            diag[a] = -(0..k)
                .filter(|&b| b != a)
                .map(|b| s[a * k + b] * pi[b])
                .sum::<f64>();
        }
        let mean_rate: f64 = -(0..k).map(|a| pi[a] * diag[a]).sum::<f64>();
        if !(mean_rate.is_finite() && mean_rate > 0.0) {
            return Err(
                "the rate matrix has no substitutions (all exchangeabilities are zero)".into(),
            );
        }

        // A = Π^½ Q Π^-½ (symmetric), scaled by 1/mean_rate.
        let mut a_mat = vec![0.0; k * k];
        for a in 0..k {
            for b in 0..k {
                a_mat[a * k + b] = if a == b {
                    diag[a] / mean_rate
                } else {
                    (pi[a] * pi[b]).sqrt() * s[a * k + b] / mean_rate
                };
            }
        }
        let (mut lambda, r) = jacobi_eigen(&mut a_mat, k);
        // The stationary eigenvalue is exactly zero in theory; make it so.
        for l in lambda.iter_mut() {
            if l.abs() < 1e-12 {
                *l = 0.0;
            }
        }
        if lambda.iter().filter(|&&l| l == 0.0).count() != 1 || lambda.iter().any(|&l| l > 0.0) {
            return Err(
                "the rate matrix is not irreducible (expected one zero eigenvalue, the rest negative)".into(),
            );
        }

        let mut x = vec![0.0; k * k * k];
        for a in 0..k {
            for b in 0..k {
                let scale = (pi[b] / pi[a]).sqrt();
                for m in 0..k {
                    x[(a * k + b) * k + m] = r[a * k + m] * r[b * k + m] * scale;
                }
            }
        }
        Ok(Self { k, pi, lambda, x })
    }

    /// Parses the PAML/PHYLIP-style text layout used for empirical matrices: the strict
    /// lower triangle of the exchangeability matrix by rows (`K−1` rows of `1, 2, …, K−1`
    /// numbers), then `K` state frequencies. Anything after those numbers is ignored.
    pub fn from_paml(text: &str, k: usize) -> Result<Self, String> {
        let need = k * (k - 1) / 2 + k;
        let mut numbers = Vec::with_capacity(need);
        for token in text.split_whitespace() {
            if numbers.len() == need {
                break;
            }
            match token.parse::<f64>() {
                Ok(v) => numbers.push(v),
                Err(_) => {
                    return Err(format!(
                        "unexpected token '{token}' while reading the matrix"
                    ));
                }
            }
        }
        if numbers.len() < need {
            return Err(format!(
                "expected {need} numbers for {k} states, found {}",
                numbers.len()
            ));
        }
        let mut s = vec![0.0; k * k];
        let mut it = numbers.iter();
        for a in 1..k {
            for b in 0..a {
                let v = *it.next().expect("length checked above");
                s[a * k + b] = v;
                s[b * k + a] = v;
            }
        }
        let freqs: Vec<f64> = it.copied().collect();
        Self::new(&s, &freqs)
    }

    /// The equal-rates model on `K` states (Jukes–Cantor for 4, Poisson for 20).
    pub fn equal_rates(k: usize) -> Self {
        let mut s = vec![1.0; k * k];
        for a in 0..k {
            s[a * k + a] = 0.0;
        }
        Self::new(&s, &vec![1.0; k]).expect("equal-rates model is valid for k >= 2")
    }

    /// Number of states `K`.
    pub fn n_states(&self) -> usize {
        self.k
    }

    /// Stationary frequencies `π` (sum to one).
    pub fn frequencies(&self) -> &[f64] {
        &self.pi
    }

    /// Eigenvalues of the rate matrix (one is exactly zero, the others negative).
    pub fn eigenvalues(&self) -> &[f64] {
        &self.lambda
    }

    /// The `K` coefficients `X_ab[m]` of entry `(a, b)`.
    #[inline]
    pub(crate) fn coefficients(&self, a: usize, b: usize) -> &[f64] {
        let o = (a * self.k + b) * self.k;
        &self.x[o..o + self.k]
    }

    /// Per-eigenvalue time factors and their first two derivatives, for the variable
    /// sites: `g_m(t)`, `g′_m(t)`, `g″_m(t)`.
    ///
    /// Uniform rates: `g = e^{λt}`. Gamma-distributed rates with shape `α` (mean one):
    /// `g = (1 − λt/α)^{−α}`, the Laplace transform of the rate distribution
    /// [Gu & Li 1996; Waddell & Steel 1997].
    pub(crate) fn factors(
        &self,
        t: f64,
        gamma_shape: Option<f64>,
        g: &mut [f64],
        g1: &mut [f64],
        g2: &mut [f64],
    ) {
        match gamma_shape {
            None => {
                for (m, &l) in self.lambda.iter().enumerate() {
                    let e = (l * t).exp();
                    g[m] = e;
                    g1[m] = l * e;
                    g2[m] = l * l * e;
                }
            }
            Some(alpha) => {
                for (m, &l) in self.lambda.iter().enumerate() {
                    let base = 1.0 - l * t / alpha;
                    let v = base.powf(-alpha);
                    g[m] = v;
                    g1[m] = l * v / base;
                    g2[m] = l * l * (1.0 + 1.0 / alpha) * v / (base * base);
                }
            }
        }
    }

    /// The factors in the limit `t → ∞`: `1` for the stationary eigenvalue, `0` otherwise.
    pub(crate) fn limit_factors(&self, g: &mut [f64]) {
        for (m, &l) in self.lambda.iter().enumerate() {
            g[m] = if l == 0.0 { 1.0 } else { 0.0 };
        }
    }

    /// Transition probability matrix `P(t)` (row-major `K×K`), with among-site rate
    /// variation: invariant sites (`p_invar`) mix in the identity, gamma rates change the
    /// time factors. `t` is the expected number of substitutions per *variable* site.
    pub fn transition_matrix(&self, t: f64, het: RateHet) -> Vec<f64> {
        let k = self.k;
        let (mut g, mut g1, mut g2) = (vec![0.0; k], vec![0.0; k], vec![0.0; k]);
        self.factors(t, het.gamma_shape, &mut g, &mut g1, &mut g2);
        let mut p = vec![0.0; k * k];
        for a in 0..k {
            for b in 0..k {
                let v: f64 = self
                    .coefficients(a, b)
                    .iter()
                    .zip(&g)
                    .map(|(x, g)| x * g)
                    .sum();
                p[a * k + b] = (1.0 - het.p_invar) * v + if a == b { het.p_invar } else { 0.0 };
            }
        }
        p
    }
}

/// Cyclic Jacobi eigenvalue algorithm for a real symmetric `n×n` matrix (row-major).
///
/// Returns the eigenvalues and the matrix whose *columns* are the corresponding
/// orthonormal eigenvectors. `a` is destroyed. Converges quadratically; a 20×20 matrix
/// needs about eight sweeps [Press et al. 2007, §11.1].
pub(crate) fn jacobi_eigen(a: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut eigenvalues = vec![0.0; n];
    let mut v = vec![0.0; n * n];
    jacobi_eigen_into(a, n, &mut eigenvalues, &mut v);
    (eigenvalues, v)
}

/// [`jacobi_eigen`] writing into caller-provided buffers (`eigenvalues`: `n` values, `v`: `n×n`), so
/// that a caller with a small fixed size, such as the 4×4 problem of every DNA pair, needs no heap.
pub(crate) fn jacobi_eigen_into(a: &mut [f64], n: usize, eigenvalues: &mut [f64], v: &mut [f64]) {
    v[..n * n].fill(0.0);
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    let scale: f64 = a
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
        .max(f64::MIN_POSITIVE);
    for _sweep in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off.sqrt() <= 1e-15 * scale {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    for i in 0..n {
        eigenvalues[i] = a[i * n + i];
    }
}

#[cfg(test)]
pub(crate) mod testing {
    use super::ReversibleModel;

    /// Small deterministic generator (xorshift64*), so tests need no dependency.
    pub(crate) struct Rng(pub u64);
    impl Rng {
        pub fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        pub fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    pub(crate) fn random_model(rng: &mut Rng, k: usize) -> ReversibleModel {
        let mut s = vec![0.0; k * k];
        for a in 0..k {
            for b in 0..a {
                let v = 0.05 + 2.0 * rng.unit() * rng.unit();
                s[a * k + b] = v;
                s[b * k + a] = v;
            }
        }
        let pi: Vec<f64> = (0..k).map(|_| 0.2 + rng.unit()).collect();
        ReversibleModel::new(&s, &pi).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::testing::{Rng, random_model};
    use super::*;

    #[test]
    fn test_jacobi_reconstructs_and_is_orthogonal() {
        let mut rng = Rng(0xA11CE);
        for n in [2usize, 3, 4, 7, 20] {
            let mut m = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..=i {
                    let v = rng.unit() * 2.0 - 1.0;
                    m[i * n + j] = v;
                    m[j * n + i] = v;
                }
            }
            let orig = m.clone();
            let (lam, v) = jacobi_eigen(&mut m, n);
            for i in 0..n {
                for j in 0..n {
                    let rec: f64 = (0..n).map(|k| v[i * n + k] * lam[k] * v[j * n + k]).sum();
                    assert!(
                        (rec - orig[i * n + j]).abs() < 1e-12,
                        "n {n}: reconstruction ({i},{j})"
                    );
                    let dot: f64 = (0..n).map(|k| v[k * n + i] * v[k * n + j]).sum();
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!((dot - want).abs() < 1e-12, "n {n}: orthogonality ({i},{j})");
                }
            }
        }
    }

    #[test]
    fn test_transition_matrix_is_a_reversible_markov_semigroup() {
        let mut rng = Rng(0xBEEF);
        for k in [4usize, 20] {
            let model = random_model(&mut rng, k);
            let pi = model.frequencies();
            let p0 = model.transition_matrix(0.0, RateHet::NONE);
            for a in 0..k {
                for b in 0..k {
                    let want = if a == b { 1.0 } else { 0.0 };
                    assert!((p0[a * k + b] - want).abs() < 1e-12, "P(0) is the identity");
                }
            }
            for &t in &[0.01, 0.3, 1.0, 4.0] {
                let p = model.transition_matrix(t, RateHet::NONE);
                for a in 0..k {
                    let row: f64 = (0..k).map(|b| p[a * k + b]).sum();
                    assert!((row - 1.0).abs() < 1e-12, "rows sum to one (k {k}, t {t})");
                    for b in 0..k {
                        assert!(p[a * k + b] > -1e-14, "non-negative");
                        let lhs = pi[a] * p[a * k + b];
                        let rhs = pi[b] * p[b * k + a];
                        assert!((lhs - rhs).abs() < 1e-13, "detailed balance (k {k}, t {t})");
                    }
                }
                // Chapman–Kolmogorov: P(t) P(t) = P(2t).
                let p2 = model.transition_matrix(2.0 * t, RateHet::NONE);
                for a in 0..k {
                    for b in 0..k {
                        let prod: f64 = (0..k).map(|c| p[a * k + c] * p[c * k + b]).sum();
                        assert!(
                            (prod - p2[a * k + b]).abs() < 1e-12,
                            "semigroup (k {k}, t {t})"
                        );
                    }
                }
            }
            let far = model.transition_matrix(1e4, RateHet::NONE);
            for a in 0..k {
                for b in 0..k {
                    assert!((far[a * k + b] - pi[b]).abs() < 1e-10, "P(t) tends to π");
                }
            }
        }
    }

    #[test]
    fn test_time_scale_is_one_substitution_per_unit_time() {
        let mut rng = Rng(0xC0FFEE);
        for k in [4usize, 20] {
            let model = random_model(&mut rng, k);
            let t = 1e-6;
            let p = model.transition_matrix(t, RateHet::NONE);
            let moved: f64 = (0..k)
                .map(|a| model.frequencies()[a] * (1.0 - p[a * k + a]))
                .sum();
            assert!(
                (moved / t - 1.0).abs() < 1e-5,
                "k {k}: mean rate {}",
                moved / t
            );
        }
    }

    #[test]
    fn test_equal_rates_match_jukes_cantor_and_poisson() {
        for k in [4usize, 20] {
            let model = ReversibleModel::equal_rates(k);
            let kf = k as f64;
            for &t in &[0.05, 0.7, 2.5] {
                let p = model.transition_matrix(t, RateHet::NONE);
                let same = 1.0 / kf + (kf - 1.0) / kf * (-kf * t / (kf - 1.0)).exp();
                let diff = (1.0 - same) / (kf - 1.0);
                for a in 0..k {
                    for b in 0..k {
                        let want = if a == b { same } else { diff };
                        assert!((p[a * k + b] - want).abs() < 1e-13, "k {k}, t {t}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_derivatives_match_finite_differences() {
        let mut rng = Rng(0xD00D);
        let model = random_model(&mut rng, 20);
        let k = 20;
        for het in [
            RateHet::NONE,
            RateHet {
                gamma_shape: Some(0.5),
                p_invar: 0.0,
            },
            RateHet {
                gamma_shape: Some(2.0),
                p_invar: 0.3,
            },
        ] {
            let (mut g, mut g1, mut g2) = (vec![0.0; k], vec![0.0; k], vec![0.0; k]);
            let t = 0.8;
            let h = 1e-5;
            let (mut gp, mut gm) = (vec![0.0; k], vec![0.0; k]);
            let (mut d1p, mut d1m) = (vec![0.0; k], vec![0.0; k]);
            let mut scratch = vec![0.0; k];
            model.factors(t, het.gamma_shape, &mut g, &mut g1, &mut g2);
            model.factors(t + h, het.gamma_shape, &mut gp, &mut d1p, &mut scratch);
            model.factors(t - h, het.gamma_shape, &mut gm, &mut d1m, &mut scratch);
            for m in 0..k {
                let fd1 = (gp[m] - gm[m]) / (2.0 * h);
                let fd2 = (d1p[m] - d1m[m]) / (2.0 * h);
                assert!(
                    (fd1 - g1[m]).abs() < 1e-7 * (1.0 + g1[m].abs()),
                    "g' at {m}, {het:?}"
                );
                assert!(
                    (fd2 - g2[m]).abs() < 1e-6 * (1.0 + g2[m].abs()),
                    "g'' at {m}, {het:?}"
                );
            }
        }
    }

    #[test]
    fn test_gamma_factor_is_the_laplace_transform_of_the_rate_distribution() {
        // E[exp(λ r t)] over r ~ Gamma(shape α, mean 1), by direct numerical integration.
        let (alpha, lambda, t) = (1.7f64, -1.3f64, 0.9f64);
        let ln_norm = alpha * alpha.ln() - ln_gamma(alpha);
        let (steps, upper) = (400_000usize, 80.0f64);
        let h = upper / steps as f64;
        let mut sum = 0.0;
        for i in 0..steps {
            let r = (i as f64 + 0.5) * h;
            let density = (ln_norm + (alpha - 1.0) * r.ln() - alpha * r).exp();
            sum += density * (lambda * r * t).exp() * h;
        }
        let closed = (1.0 - lambda * t / alpha).powf(-alpha);
        assert!(
            (sum - closed).abs() < 1e-6,
            "numerical {sum} vs closed form {closed}"
        );
    }

    /// Lanczos approximation of ln Γ(x), enough for the test above.
    fn ln_gamma(x: f64) -> f64 {
        const G: [f64; 6] = [
            76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5,
        ];
        let mut y = x;
        let tmp = x + 5.5 - (x + 0.5) * (x + 5.5).ln();
        let mut ser = 1.000000000190015;
        for g in G {
            y += 1.0;
            ser += g / y;
        }
        -tmp + (2.5066282746310005 * ser / x).ln()
    }

    #[test]
    fn test_paml_parser_round_trips_and_rejects_bad_input() {
        let k = 4;
        let text = "1\n2 3\n4 5 6\n\n0.1 0.2 0.3 0.4\n\nnotes after the numbers are ignored";
        let model = ReversibleModel::from_paml(text, k).unwrap();
        let mut s = vec![0.0; 16];
        let mut v = 1.0;
        for a in 1..k {
            for b in 0..a {
                s[a * k + b] = v;
                s[b * k + a] = v;
                v += 1.0;
            }
        }
        let direct = ReversibleModel::new(&s, &[0.1, 0.2, 0.3, 0.4]).unwrap();
        let (p, q) = (
            model.transition_matrix(0.4, RateHet::NONE),
            direct.transition_matrix(0.4, RateHet::NONE),
        );
        assert_eq!(p, q);
        assert!(ReversibleModel::from_paml("1 2 3", 4).is_err());
        assert!(ReversibleModel::from_paml("1 2 x 4 5 6 .1 .2 .3 .4", 4).is_err());
    }

    #[test]
    fn test_invalid_models_are_rejected() {
        let s = vec![0.0, 1.0, 1.0, 0.0];
        assert!(
            ReversibleModel::new(&s, &[0.5, 0.0]).is_err(),
            "zero frequency"
        );
        assert!(ReversibleModel::new(&s, &[0.5]).is_err(), "one state");
        assert!(
            ReversibleModel::new(&[0.0, 1.0, 2.0, 0.0], &[0.5, 0.5]).is_err(),
            "asymmetric"
        );
        assert!(
            ReversibleModel::new(&[0.0; 4], &[0.5, 0.5]).is_err(),
            "no substitutions"
        );
        // Two disconnected pairs of states: more than one zero eigenvalue.
        let mut blocks = vec![0.0; 16];
        blocks[1] = 1.0;
        blocks[4] = 1.0;
        blocks[11] = 1.0;
        blocks[14] = 1.0;
        assert!(
            ReversibleModel::new(&blocks, &[0.25; 4]).is_err(),
            "reducible"
        );
    }
}
