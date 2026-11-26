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
}
