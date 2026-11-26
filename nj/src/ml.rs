//! Maximum-likelihood distance between two sequences under a reversible model.
//!
//! The sufficient statistic of a pair is its table of aligned state pairs `F` (`F_ab`
//! columns with `a` in the first and `b` in the second sequence). With `P(t)` from a
//! [`ReversibleModel`] the log-likelihood of a distance `t` is
//!
//! ```text
//! ℓ(t) = Σ_ab F_ab ln P_ab(t)        ℓ′ = Σ F_ab P′_ab / P_ab        ℓ″ = Σ F_ab (P″_ab P_ab − P′_ab²) / P_ab²
//! ```
//!
//! (the stationary term `Σ F_ab ln π_a` does not depend on `t`), and the distance is
//! `d̂ = argmax ℓ`. The model is fixed; only `t` is fitted, so this is a one-dimensional
//! problem per pair. Background: `docs/theory/distance-estimation.md`.
//!
//! # Unordered pairs
//!
//! Reversibility gives `π_a P_ab(t) = π_b P_ba(t)`, so `F_ab ln P_ab + F_ba ln P_ba =
//! (F_ab + F_ba) ln P_ab + F_ba ln(π_a/π_b)`: the second term does not depend on `t`. Only the
//! merged counts of unordered pairs `a ≤ b` are needed, which roughly halves the work per
//! evaluation. `ℓ` below is this reduced log-likelihood (it differs from the full one by a
//! constant), and `ℓ(∞)` is reduced in the same way.
//!
//! # Solver
//!
//! `ℓ′(0⁺) = +∞` whenever the pair differs, but `ℓ′` can change sign again in the tail, so the
//! likelihood is not assumed unimodal. The solver runs Newton's method on `ℓ′` from the
//! Poisson-type closed-form estimate (typically within a few evaluations of the answer) and
//! keeps the tightest bracket seen so far, a point with `ℓ′ > 0` on the left and one with
//! `ℓ′ < 0` on the right. A step that leaves the bracket, or that has the wrong curvature, is
//! replaced by bisection; while one side is still unknown the step doubles or halves `t`.
//!
//! A pair is **saturated** (`d̂ = ∞`) when `ℓ′` is still positive at [`T_MAX`], or when the
//! maximum does not beat the limit `ℓ(∞)`: the data are then no better explained by any finite
//! distance than by two independent sequences.

use crate::models::RateHet;
use crate::reversible::ReversibleModel;

/// Largest distance searched, in substitutions per variable site. Jukes–Cantor at
/// `p = 0.7499` is already 6.7, so 10 is beyond anything an estimate can resolve.
pub const T_MAX: f64 = 10.0;

/// Smallest distance considered; below this `P_ab(t)` loses its relative accuracy.
const T_MIN: f64 = 1e-9;

/// Value of `ℓ`, `ℓ′` and `ℓ″` at one `t`.
struct Eval {
    ll: f64,
    d1: f64,
    d2: f64,
}

/// Likelihood evaluation for one pair: the nonzero entries of `F` plus scratch space.
struct Pair<'a> {
    model: &'a ReversibleModel,
    het: RateHet,
    /// `(a, b, F_ab)` for every `F_ab > 0`.
    entries: Vec<(usize, usize, f64)>,
    g: Vec<f64>,
    g1: Vec<f64>,
    g2: Vec<f64>,
}

impl Pair<'_> {
    /// `ℓ(t)` and its first two derivatives, or `None` if some `P_ab(t)` is not positive
    /// (`t` so small that rounding swamps an off-diagonal entry).
    fn eval(&mut self, t: f64) -> Option<Eval> {
        let model = self.model;
        model.factors(
            t,
            self.het.gamma_shape,
            &mut self.g,
            &mut self.g1,
            &mut self.g2,
        );
        let (w, p_inv) = (1.0 - self.het.p_invar, self.het.p_invar);
        let mut e = Eval {
            ll: 0.0,
            d1: 0.0,
            d2: 0.0,
        };
        for &(a, b, n) in &self.entries {
            let x = model.coefficients(a, b);
            let (mut p, mut p1, mut p2) = (0.0, 0.0, 0.0);
            for m in 0..x.len() {
                p += x[m] * self.g[m];
                p1 += x[m] * self.g1[m];
                p2 += x[m] * self.g2[m];
            }
            let (mut p, p1, p2) = (w * p, w * p1, w * p2);
            if a == b {
                p += p_inv;
            }
            if !(p > 0.0) {
                return None;
            }
            e.ll += n * p.ln();
            e.d1 += n * p1 / p;
            e.d2 += n * (p2 * p - p1 * p1) / (p * p);
        }
        Some(e)
    }

    /// `ℓ(∞)`: every `P_ab` at its limit, `p_inv·δ_ab + (1 − p_inv)·π_b`.
    fn loglik_at_infinity(&mut self) -> f64 {
        self.model.limit_factors(&mut self.g);
        let (w, p_inv) = (1.0 - self.het.p_invar, self.het.p_invar);
        let pi = self.model.frequencies();
        self.entries
            .iter()
            .map(|&(a, b, n)| n * (w * pi[b] + if a == b { p_inv } else { 0.0 }).ln())
            .sum()
    }
}

/// Maximum-likelihood distance for the pair-count table `counts` (row-major `K×K`
/// non-negative counts, `K = model.n_states()`), in substitutions per site over all sites.
///
/// Returns `0.0` when no column differs (including when nothing is comparable) and
/// `f64::INFINITY` when the pair is saturated. Under rate variation the result is
/// `(1 − p_invar)·t̂`, where `t̂` is the fitted distance per variable site.
pub fn ml_distance(model: &ReversibleModel, counts: &[f64], het: RateHet) -> f64 {
    let k = model.n_states();
    debug_assert_eq!(counts.len(), k * k);

    // Unordered pairs a <= b (see the module docs); `total` and `differing` use the raw table.
    let mut entries = Vec::new();
    let (mut total, mut differing) = (0.0, 0.0);
    for a in 0..k {
        let n = counts[a * k + a];
        if n > 0.0 {
            entries.push((a, a, n));
            total += n;
        }
        for b in (a + 1)..k {
            let n = counts[a * k + b] + counts[b * k + a];
            if n > 0.0 {
                entries.push((a, b, n));
                total += n;
                differing += n;
            }
        }
    }
    if differing == 0.0 {
        return 0.0;
    }

    let mut pair = Pair {
        model,
        het,
        entries,
        g: vec![0.0; k],
        g1: vec![0.0; k],
        g2: vec![0.0; k],
    };
    let Some(t) = solve(&mut pair, total, differing) else {
        return f64::INFINITY;
    };
    (1.0 - het.p_invar) * t
}

/// Finds the interior maximum of `ℓ`, or `None` if the pair is saturated.
fn solve(pair: &mut Pair<'_>, total: f64, differing: f64) -> Option<f64> {
    let k = pair.model.n_states() as f64;
    let b = (k - 1.0) / k;

    // Closed-form start (equal rates, variable sites only).
    let p = differing / total / (1.0 - pair.het.p_invar);
    let start = if p < b { -b * (1.0 - p / b).ln() } else { 1.0 };

    // Bracket: ℓ′ > 0 at `lo` (0 = nothing known yet), ℓ′ < 0 at `hi` (∞ = nothing known yet).
    let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
    let mut t = start.clamp(10.0 * T_MIN, T_MAX / 2.0);
    let mut ll = f64::NEG_INFINITY;
    for _ in 0..200 {
        let mut next;
        match pair.eval(t) {
            None => {
                // t is too small to evaluate: the slope there is positive.
                lo = lo.max(t);
                next = if hi.is_finite() {
                    0.5 * (lo + hi)
                } else {
                    2.0 * t
                };
            }
            Some(e) => {
                ll = e.ll;
                if e.d1 == 0.0 {
                    next = t;
                } else {
                    if e.d1 > 0.0 {
                        lo = lo.max(t);
                    } else {
                        hi = hi.min(t);
                    }
                    let newton = t - e.d1 / e.d2;
                    let usable = e.d2 < 0.0 && newton > lo && newton < hi;
                    next = if usable {
                        newton
                    } else if hi.is_finite() && lo > 0.0 {
                        0.5 * (lo + hi)
                    } else if hi.is_finite() {
                        (0.5 * t).max(T_MIN) // only a right end known: go left
                    } else {
                        2.0 * t // only a left end known: go right
                    };
                }
            }
        }
        if next >= T_MAX {
            if lo >= T_MAX {
                return None; // still rising at t_max: no interior maximum
            }
            next = T_MAX;
        }
        if next <= T_MIN {
            if hi <= T_MIN {
                return Some(T_MIN); // falling from the smallest resolvable distance
            }
            next = T_MIN;
        }
        let converged = (next - t).abs() <= 1e-12 + 1e-8 * t;
        t = next;
        if converged {
            break;
        }
    }

    // The maximum must beat "the sequences are independent" (ℓ from the last evaluation,
    // which is at a point within the convergence tolerance of `t`).
    (ll > pair.loglik_at_infinity()).then_some(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::DNA;
    use crate::models::{JukesCantor, pairwise_distance_with};
    use crate::reversible::testing::{Rng, random_model};

    fn dense(k: usize, entries: &[(usize, usize, f64)]) -> Vec<f64> {
        let mut f = vec![0.0; k * k];
        for &(a, b, n) in entries {
            f[a * k + b] += n;
        }
        f
    }

    /// Expected pair counts `N·π_a·P_ab(t)` (real-valued), so the exact maximum is `t`.
    fn expected_counts(model: &ReversibleModel, t: f64, het: RateHet, n: f64) -> Vec<f64> {
        let k = model.n_states();
        let p = model.transition_matrix(t, het);
        let mut f = vec![0.0; k * k];
        for a in 0..k {
            for b in 0..k {
                f[a * k + b] = n * model.frequencies()[a] * p[a * k + b];
            }
        }
        f
    }

    #[test]
    fn test_recovers_the_true_distance_from_expected_counts() {
        let mut rng = Rng(0x5EED_0010);
        let hets = [
            RateHet::NONE,
            RateHet {
                gamma_shape: Some(0.5),
                p_invar: 0.0,
            },
            RateHet {
                gamma_shape: None,
                p_invar: 0.25,
            },
            RateHet {
                gamma_shape: Some(1.5),
                p_invar: 0.4,
            },
        ];
        for k in [4usize, 20] {
            for _ in 0..3 {
                let model = random_model(&mut rng, k);
                for het in hets {
                    for &t in &[0.003, 0.05, 0.4, 1.0, 2.0] {
                        let f = expected_counts(&model, t, het, 1e8);
                        let d = ml_distance(&model, &f, het);
                        let want = (1.0 - het.p_invar) * t;
                        assert!(
                            (d - want).abs() < 1e-6 * want.max(1e-3),
                            "k {k}, {het:?}, t {t}: got {d}, want {want}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_equal_rates_give_the_closed_forms() {
        // JC (K=4) and Poisson (K=20): the ML distance depends on the differing fraction only.
        for k in [4usize, 20] {
            let model = ReversibleModel::equal_rates(k);
            let b = (k as f64 - 1.0) / k as f64;
            for &(same, diff) in &[(990.0, 10.0), (700.0, 300.0), (450.0, 400.0), (5000.0, 1.0)] {
                let f = dense(k, &[(0, 0, same), (0, 1, diff * 0.3), (3, 2, diff * 0.7)]);
                let p = diff / (same + diff);
                if p >= b {
                    continue;
                }
                let want = -b * (1.0 - p / b).ln();
                let got = ml_distance(&model, &f, RateHet::NONE);
                assert!((got - want).abs() < 1e-9, "k {k}, p {p}: {got} vs {want}");
            }
        }
    }

    /// Two sequences with `same` identical columns and `diff` columns A/C, as symbols.
    fn dna_pair(
        same: usize,
        diff: usize,
    ) -> (
        Vec<crate::alphabet::DnaSymbol>,
        Vec<crate::alphabet::DnaSymbol>,
    ) {
        use crate::alphabet::{AlphabetEncoding, DNA};
        let mut x: Vec<_> = vec![DNA::encode(b'A'); same];
        let mut y = x.clone();
        x.extend(vec![DNA::encode(b'A'); diff]);
        y.extend(vec![DNA::encode(b'C'); diff]);
        (x, y)
    }

    #[test]
    fn test_matches_the_existing_jukes_cantor_formulas_with_rate_variation() {
        // The existing JC distance has closed forms for gamma rates and invariant sites;
        // the equal-rates ML solver must reproduce each of them.
        let model = ReversibleModel::equal_rates(4);
        let (x, y) = dna_pair(800, 120);
        let f = dense(4, &[(0, 0, 800.0), (0, 1, 120.0)]);
        for het in [
            RateHet::NONE,
            RateHet {
                gamma_shape: Some(0.3),
                p_invar: 0.0,
            },
            RateHet {
                gamma_shape: Some(2.0),
                p_invar: 0.0,
            },
            RateHet {
                gamma_shape: None,
                p_invar: 0.2,
            },
            RateHet {
                gamma_shape: Some(0.8),
                p_invar: 0.35,
            },
        ] {
            let want = pairwise_distance_with::<JukesCantor, DNA>(&x, &y, het);
            let got = ml_distance(&model, &f, het);
            assert!(
                (got - want).abs() < 1e-8 * want,
                "{het:?}: ML {got} vs closed form {want}"
            );
        }
    }

    /// Closed form for the equal-rates model on `K` states with rate variation:
    /// with `b = (K−1)/K` and `p` the differing fraction among variable sites,
    /// `d = −b·ln(1 − p/b)` (uniform) or `b·α·((1 − p/b)^(−1/α) − 1)` (gamma), times `1 − p_invar`.
    fn equal_rates_closed_form(k: usize, p_observed: f64, het: RateHet) -> f64 {
        let b = (k as f64 - 1.0) / k as f64;
        let p = p_observed / (1.0 - het.p_invar);
        let arg = 1.0 - p / b;
        let t = match het.gamma_shape {
            None => -b * arg.ln(),
            Some(alpha) => b * alpha * (arg.powf(-1.0 / alpha) - 1.0),
        };
        (1.0 - het.p_invar) * t
    }

    #[test]
    fn test_equal_rates_on_20_states_with_rate_variation() {
        // The library's `Poisson` is −ln(1 − p), the infinite-states limit, which is *not*
        // the 20-state equal-rates model (−(19/20)·ln(1 − 20p/19)); the solver must give the
        // latter, with the same gamma / invariant-sites generalisation as Jukes–Cantor.
        let model = ReversibleModel::equal_rates(20);
        let f = dense(20, &[(0, 0, 700.0), (0, 17, 150.0)]);
        let p = 150.0 / 850.0;
        for het in [
            RateHet::NONE,
            RateHet {
                gamma_shape: Some(0.5),
                p_invar: 0.0,
            },
            RateHet {
                gamma_shape: Some(1.0),
                p_invar: 0.3,
            },
        ] {
            let want = equal_rates_closed_form(20, p, het);
            let got = ml_distance(&model, &f, het);
            assert!(
                (got - want).abs() < 1e-8 * want,
                "{het:?}: ML {got} vs closed form {want}"
            );
        }
        // The two differ by about 0.5% here; guard against confusing them.
        let poisson = -(1.0_f64 - p).ln();
        assert!((ml_distance(&model, &f, RateHet::NONE) - poisson).abs() > 5e-4);
    }

    /// Draws `n` aligned pairs from `π_a P_ab(t)` and tabulates them (row-major `K×K`).
    fn simulate_counts(rng: &mut Rng, model: &ReversibleModel, t: f64, n: usize) -> Vec<f64> {
        let k = model.n_states();
        let p = model.transition_matrix(t, RateHet::NONE);
        let pick = |rng: &mut Rng, weights: &mut dyn Iterator<Item = f64>| {
            let u = rng.unit();
            let mut acc = 0.0;
            let mut last = 0;
            for (i, w) in weights.enumerate() {
                acc += w;
                last = i;
                if u < acc {
                    return i;
                }
            }
            last
        };
        let mut f = vec![0.0; k * k];
        for _ in 0..n {
            let a = pick(rng, &mut model.frequencies().iter().copied());
            let b = pick(rng, &mut p[a * k..(a + 1) * k].iter().copied());
            f[a * k + b] += 1.0;
        }
        f
    }

    /// The full (unmerged) log-likelihood on a grid, the independent reference for the solver.
    fn grid_maximum(
        model: &ReversibleModel,
        f: &[f64],
        step: f64,
        steps: usize,
    ) -> (f64, f64, f64) {
        let k = model.n_states();
        let ll = |p: &[f64]| -> f64 {
            (0..k * k)
                .filter(|&i| f[i] > 0.0)
                .map(|i| f[i] * p[i].ln())
                .sum()
        };
        let (mut best_t, mut best_ll) = (0.0, f64::NEG_INFINITY);
        for i in 1..=steps {
            let t = step * i as f64;
            let v = ll(&model.transition_matrix(t, RateHet::NONE));
            if v > best_ll {
                (best_t, best_ll) = (t, v);
            }
        }
        let limit: f64 = (0..k * k)
            .filter(|&i| f[i] > 0.0)
            .map(|i| f[i] * model.frequencies()[i % k].ln())
            .sum();
        (best_t, best_ll, limit)
    }

    /// Compares the solver with the grid maximum of the full likelihood on simulated pairs.
    /// Short sequences at high divergence are where the likelihood is flat, has tails that
    /// rise towards ℓ(∞), or has no interior maximum at all. Returns (finite, saturated).
    fn check_against_grid(k: usize, trials: usize, sites: usize, seed: u64) -> (usize, usize) {
        let mut rng = Rng(seed);
        let (mut finite, mut saturated) = (0, 0);
        for trial in 0..trials {
            let model = random_model(&mut rng, k);
            let t = [0.3, 1.0, 2.0, 3.0, 4.5][trial % 5];
            let f = simulate_counts(&mut rng, &model, t, sites);
            let d = ml_distance(&model, &f, RateHet::NONE);
            let (grid_t, grid_ll, limit) = grid_maximum(&model, &f, 0.01, 1000);
            let differing: f64 = (0..k * k).filter(|i| i / k != i % k).map(|i| f[i]).sum();
            if differing == 0.0 {
                assert_eq!(d, 0.0);
                continue;
            }
            let loglik = |t: f64| -> f64 {
                let p = model.transition_matrix(t, RateHet::NONE);
                (0..k * k)
                    .filter(|&i| f[i] > 0.0)
                    .map(|i| f[i] * p[i].ln())
                    .sum()
            };
            // A maximum on the grid's last point means ℓ is still rising at t_max: saturated by
            // definition. A margin keeps knife-edge cases (within rounding of ℓ(∞)) out of the test.
            if grid_t >= 9.995 {
                saturated += 1;
                assert!(
                    d.is_infinite(),
                    "k {k} trial {trial}: still rising at t_max but solver returned {d}"
                );
            } else if grid_ll > limit + 1e-6 {
                finite += 1;
                assert!(
                    d.is_finite(),
                    "k {k} trial {trial}: grid finds a maximum at {grid_t}, solver says saturated"
                );
                // Two maxima of nearly equal height may swap places: accept any maximum whose
                // likelihood is within a hair of the grid's best.
                assert!(
                    (d - grid_t).abs() < 0.011 || loglik(d) >= grid_ll - 1e-3,
                    "k {k} trial {trial}: solver {d}, grid {grid_t}"
                );
            } else if grid_ll < limit - 1e-6 {
                saturated += 1;
                assert!(
                    d.is_infinite(),
                    "k {k} trial {trial}: no maximum beats ℓ(∞) but solver returned {d}"
                );
            }
        }
        (finite, saturated)
    }

    #[test]
    fn test_solver_agrees_with_a_grid_search_on_simulated_pairs() {
        let (finite, saturated) = check_against_grid(4, 800, 100, 0x5EED_0030);
        assert!(
            finite > 300 && saturated > 30,
            "both outcomes must occur: {finite} finite, {saturated} saturated"
        );
        // Saturation is rare at K=20 and 300 sites (about 2%); the extensive test covers it.
        let (finite, _) = check_against_grid(20, 30, 300, 0x5EED_0031);
        assert!(
            finite > 20,
            "K=20: only {finite} pairs had a finite maximum"
        );
    }

    /// Larger sweep: `cargo test --release -p nj --lib ml::tests::extensive -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn extensive_grid_comparison() {
        let (f4, s4) = check_against_grid(4, 60_000, 100, 0xABCD_0001);
        let (f20, s20) = check_against_grid(20, 2_500, 300, 0xABCD_0002);
        println!("K=4: {f4} finite, {s4} saturated; K=20: {f20} finite, {s20} saturated");
    }

    #[test]
    fn test_identical_and_empty_pairs_are_zero() {
        let model = ReversibleModel::equal_rates(20);
        assert_eq!(
            ml_distance(
                &model,
                &dense(20, &[(2, 2, 40.0), (5, 5, 10.0)]),
                RateHet::NONE
            ),
            0.0
        );
        assert_eq!(ml_distance(&model, &vec![0.0; 400], RateHet::NONE), 0.0);
    }

    #[test]
    fn test_saturated_pairs_are_infinite_not_nan() {
        let mut rng = Rng(0x5EED_0011);
        for k in [4usize, 20] {
            let model = ReversibleModel::equal_rates(k);
            let b = (k as f64 - 1.0) / k as f64;
            // At or beyond the closed-form saturation point.
            for frac in [b, b + 0.05, 0.99] {
                let diff = (frac * 1000.0).ceil();
                let f = dense(k, &[(0, 0, 1000.0 - diff), (0, 1, diff)]);
                let d = ml_distance(&model, &f, RateHet::NONE);
                assert!(d.is_infinite() && d > 0.0, "k {k}, p {frac}: got {d}");
            }
            // Random models, nothing but mismatches.
            let model = random_model(&mut rng, k);
            let f = dense(k, &[(0, 1, 500.0), (1, 0, 500.0)]);
            let d = ml_distance(&model, &f, RateHet::NONE);
            assert!(!d.is_nan(), "never NaN");
        }
    }

    #[test]
    fn test_distance_is_monotone_in_the_number_of_differences() {
        let model = ReversibleModel::equal_rates(20);
        let mut last = 0.0;
        for diff in (1..=900).step_by(37) {
            let f = dense(20, &[(1, 1, 1000.0 - diff as f64), (1, 4, diff as f64)]);
            let d = ml_distance(&model, &f, RateHet::NONE);
            if d.is_infinite() {
                break;
            }
            assert!(d > last, "d({diff}) = {d} should exceed {last}");
            last = d;
        }
        assert!(last > 1.0);
    }

    #[test]
    fn test_tiny_distances_are_resolved() {
        // One difference in a million sites.
        let model = ReversibleModel::equal_rates(20);
        let f = dense(20, &[(0, 0, 999_999.0), (0, 1, 1.0)]);
        let want = -(19.0_f64 / 20.0) * (1.0 - 1e-6 / (19.0_f64 / 20.0)).ln();
        let got = ml_distance(&model, &f, RateHet::NONE);
        // Cancellation in the eigen-expansion of P_ab(t) limits the relative accuracy at
        // t = 1e-6 to a few 1e-9, far below the statistical error of a single difference.
        assert!((got - want).abs() < 1e-7 * want, "{got} vs {want}");
    }
}
