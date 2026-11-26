//! SIMD kernels, written once against [`fearless_simd`] (stable Rust, safe API, run-time
//! dispatch to the widest instruction set the CPU offers) next to the scalar reference
//! each one must reproduce.
//!
//! Each kernel is a generic function over `S: Simd` marked `#[simd]`, which makes the body
//! inline into the dispatched, `target_feature`-enabled wrapper. Without that attribute
//! the code is still correct but loses most of its speed.
//!
//! The `*_at` entry points take the [`Level`] explicitly so tests can run every backend
//! (and the scalar fallback) on the same input.

use fearless_simd::{Level, Simd, dispatch, prelude::*};
use fearless_simd_macros::simd;

/// Scalar reference: `(n_diff, n_comparable)` over two equal-length byte sequences, where
/// a byte below `n_states` is a real state and anything else (unknown, gap) is missing.
///
/// Written without data-dependent branches (boolean masks accumulated as `0`/`1`).
#[inline]
pub(crate) fn count_diff_comparable_scalar(a: &[u8], b: &[u8], n_states: u8) -> (usize, usize) {
    let mut n_comparable = 0usize;
    let mut n_diff = 0usize;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let comparable = (x < n_states) & (y < n_states);
        let diff = comparable & (x != y);
        n_comparable += comparable as usize;
        n_diff += diff as usize;
    }
    (n_diff, n_comparable)
}

/// `(n_diff, n_comparable)` over two byte sequences; same contract as
/// [`count_diff_comparable_scalar`].
#[inline]
pub(crate) fn count_diff_comparable(a: &[u8], b: &[u8], n_states: u8) -> (usize, usize) {
    count_diff_comparable_at(Level::new(), a, b, n_states)
}

#[inline]
pub(crate) fn count_diff_comparable_at(level: Level, a: &[u8], b: &[u8], n_states: u8) -> (usize, usize) {
    dispatch!(level, simd => count_diff_comparable_impl(simd, a, b, n_states))
}

/// Lane masks are accumulated into signed byte counters (a true lane is `-1`, so
/// subtracting it adds one) instead of being turned into a scalar bitmask every step;
/// NEON has no movemask, so that conversion is emulated and costly. The counters are
/// flushed to `usize` before they can overflow (`i8::MAX` steps).
#[simd]
fn count_diff_comparable_impl<S: Simd>(simd: S, a: &[u8], b: &[u8], n_states: u8) -> (usize, usize) {
    let n = a.len().min(b.len());
    let w = S::u8s::LEN;
    let limit = S::u8s::splat(simd, n_states);

    let mut n_comparable = 0usize;
    let mut n_diff = 0usize;
    let mut i = 0usize;
    while i + w <= n {
        let zero = S::i8s::splat(simd, 0);
        let (mut acc_c, mut acc_d) = (zero, zero);
        let steps = ((n - i) / w).min(i8::MAX as usize);
        for s in 0..steps {
            let o = i + s * w;
            let va = S::u8s::from_slice(simd, &a[o..o + w]);
            let vb = S::u8s::from_slice(simd, &b[o..o + w]);
            let comparable = va.simd_lt(limit) & vb.simd_lt(limit);
            let diff = comparable & va.simd_ne(vb);
            acc_c = acc_c - comparable.to_vector();
            acc_d = acc_d - diff.to_vector();
        }
        n_comparable += acc_c.as_slice().iter().map(|&x| x as usize).sum::<usize>();
        n_diff += acc_d.as_slice().iter().map(|&x| x as usize).sum::<usize>();
        i += steps * w;
    }

    let (rd, rc) = count_diff_comparable_scalar(&a[i..n], &b[i..n], n_states);
    (n_diff + rd, n_comparable + rc)
}


/// Scalar reference for [`row_minima`]: for every row `r` in `1..m` of the packed lower
/// triangle `data` (row `r` is `data[r(r-1)/2 .. r(r-1)/2 + r]`), the minimum over `t < r`
/// of the Q criterion `scale·d[r][t] − (s[r] + s[t])`; returns the minimum over all rows.
///
/// `out[0]` is set to `+inf` (row 0 has no entries) and `out[r]` to row `r`'s minimum.
/// Every element is computed with the same two separately rounded operations and `min`
/// is exact, so the result does not depend on how a kernel vectorises the loop.
pub(crate) fn row_minima_scalar(data: &[f64], s: &[f64], m: usize, scale: f64, out: &mut [f64]) -> f64 {
    let mut best = f64::INFINITY;
    out[0] = f64::INFINITY;
    for r in 1..m {
        let start = r * (r - 1) / 2;
        let (d, sr) = (&data[start..start + r], s[r]);
        let mut lanes = [f64::INFINITY; 4];
        let (cd, cs) = (d.chunks_exact(4), s[..r].chunks_exact(4));
        let (rd, rs) = (cd.remainder(), cs.remainder());
        for (x, y) in cd.zip(cs) {
            for l in 0..4 {
                lanes[l] = lanes[l].min(scale * x[l] - (sr + y[l]));
            }
        }
        let mut q = lanes[0].min(lanes[1]).min(lanes[2].min(lanes[3]));
        for (x, y) in rd.iter().zip(rs) {
            q = q.min(scale * x - (sr + y));
        }
        out[r] = q;
        best = best.min(q);
    }
    best
}

/// SIMD version of [`row_minima_scalar`] (same contract and bit-identical results).
#[inline]
pub(crate) fn row_minima(data: &[f64], s: &[f64], m: usize, scale: f64, out: &mut [f64]) -> f64 {
    row_minima_at(Level::new(), data, s, m, scale, out)
}

#[inline]
pub(crate) fn row_minima_at(level: Level, data: &[f64], s: &[f64], m: usize, scale: f64, out: &mut [f64]) -> f64 {
    dispatch!(level, simd => row_minima_impl(simd, data, s, m, scale, out))
}

/// One pass over each row slice with two independent minimum accumulators. The
/// multiply and the subtract stay separate operations (no fused multiply-add), so every
/// element rounds exactly as in the scalar reference.
#[simd]
fn row_minima_impl<S: Simd>(simd: S, data: &[f64], s: &[f64], m: usize, scale: f64, out: &mut [f64]) -> f64 {
    let w = S::f64s::LEN;
    let scale_v = S::f64s::splat(simd, scale);
    let inf = S::f64s::splat(simd, f64::INFINITY);

    let mut best = f64::INFINITY;
    out[0] = f64::INFINITY;
    for r in 1..m {
        let start = r * (r - 1) / 2;
        let d = &data[start..start + r];
        let sr = s[r];
        let sr_v = S::f64s::splat(simd, sr);

        let (mut m0, mut m1) = (inf, inf);
        let mut t = 0usize;
        while t + 2 * w <= r {
            let q0 = S::f64s::from_slice(simd, &d[t..t + w]) * scale_v
                - (sr_v + S::f64s::from_slice(simd, &s[t..t + w]));
            let q1 = S::f64s::from_slice(simd, &d[t + w..t + 2 * w]) * scale_v
                - (sr_v + S::f64s::from_slice(simd, &s[t + w..t + 2 * w]));
            m0 = m0.min(q0);
            m1 = m1.min(q1);
            t += 2 * w;
        }
        while t + w <= r {
            let q = S::f64s::from_slice(simd, &d[t..t + w]) * scale_v
                - (sr_v + S::f64s::from_slice(simd, &s[t..t + w]));
            m0 = m0.min(q);
            t += w;
        }
        let mut q = m0.min(m1).reduce_min();
        while t < r {
            q = q.min(scale * d[t] - (sr + s[t]));
            t += 1;
        }
        out[r] = q;
        best = best.min(q);
    }
    best
}

#[cfg(test)]
pub(crate) mod testing {
    //! Every instruction-set level that can run on this machine, plus the scalar fallback.

    use fearless_simd::Level;

    pub(crate) fn levels() -> Vec<(&'static str, Level)> {
        let best = Level::new();
        #[allow(unused_mut)]
        let mut v = vec![("best", best), ("fallback", Level::fallback())];
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if let Some(t) = best.as_sse2() {
                v.push(("sse2", Level::Sse2(t)));
            }
            if let Some(t) = best.as_sse4_2() {
                v.push(("sse4.2", Level::Sse4_2(t)));
            }
            if let Some(t) = best.as_avx2() {
                v.push(("avx2", Level::Avx2(t)));
            }
            if let Some(t) = best.as_avx512() {
                v.push(("avx512", Level::Avx512(t)));
            }
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::testing::levels;
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
    }

    #[test]
    fn test_count_diff_comparable_matches_scalar_at_every_level() {
        let mut rng = Rng(0x5EED_0001);
        // Lengths around vector widths (16/32/64) and beyond the 127-step counter flush.
        let lengths = [0usize, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 255, 1000, 2031, 2032, 2033, 8200, 20_000];
        for n_states in [4u8, 20] {
            for &len in &lengths {
                let a: Vec<u8> = (0..len).map(|_| (rng.next() % (n_states as u64 + 2)) as u8).collect();
                let b: Vec<u8> = a
                    .iter()
                    .map(|&x| if rng.next() % 4 == 0 { (rng.next() % (n_states as u64 + 2)) as u8 } else { x })
                    .collect();
                let want = count_diff_comparable_scalar(&a, &b, n_states);
                for (name, level) in levels() {
                    assert_eq!(
                        count_diff_comparable_at(level, &a, &b, n_states),
                        want,
                        "level {name}, len {len}, n_states {n_states}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_row_minima_match_scalar_at_every_level() {
        let mut rng = Rng(0x5EED_0002);
        for m in [1usize, 2, 3, 4, 5, 9, 17, 18, 40, 41, 130] {
            let data: Vec<f64> = (0..m * m.saturating_sub(1) / 2)
                .map(|_| 0.05 + (rng.next() >> 11) as f64 / (1u64 << 53) as f64)
                .collect();
            let s: Vec<f64> = (0..m).map(|_| 5.0 + (rng.next() >> 11) as f64 / (1u64 << 53) as f64 * 40.0).collect();
            let scale = m.saturating_sub(2) as f64;
            let mut want = vec![0.0; m];
            let want_best = row_minima_scalar(&data, &s, m, scale, &mut want);
            for (name, level) in levels() {
                let mut got = vec![0.0; m];
                let got_best = row_minima_at(level, &data, &s, m, scale, &mut got);
                assert_eq!(got_best.to_bits(), want_best.to_bits(), "level {name}, m {m}: best");
                for r in 0..m {
                    assert_eq!(got[r].to_bits(), want[r].to_bits(), "level {name}, m {m}, row {r}");
                }
            }
        }
    }
}
