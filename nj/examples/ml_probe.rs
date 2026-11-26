//! Cost of the protein ML distance (single thread, synthetic data, random 20-state model).
//!
//! ```text
//! cargo run --release -p nj --example ml_probe -- N L [reps]
//! ```
//!
//! Sequence `i` is the base sequence with a divergence that grows from 2% to `max` across the
//! alignment, so pair distances span the useful range. Prints time per distance matrix and per pair,
//! next to the same alignment under the `Poisson` closed form.

use nj::alphabet::Protein;
use nj::empirical::{ProteinMatrix, ProteinMl};
use nj::models::{Poisson, RateHet};
use nj::msa::MSA;
use nj::reversible::ReversibleModel;
use std::sync::OnceLock;
use std::time::Instant;

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
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct Random20;
impl ProteinMatrix for Random20 {
    fn model() -> &'static ReversibleModel {
        static M: OnceLock<ReversibleModel> = OnceLock::new();
        M.get_or_init(|| {
            let mut r = Rng(0xFEED_FACE_1234_5678);
            let k = 20;
            let mut s = vec![0.0; k * k];
            for a in 0..k {
                for b in 0..a {
                    let v = 0.05 + 2.0 * r.f() * r.f();
                    s[a * k + b] = v;
                    s[b * k + a] = v;
                }
            }
            let pi: Vec<f64> = (0..k).map(|_| 0.3 + r.f()).collect();
            ReversibleModel::new(&s, &pi).unwrap()
        })
    }
}

fn alignment(n: usize, l: usize, max_div: f64) -> MSA<Protein> {
    let mut r = Rng(0x9E37_79B9_7F4A_7C15);
    let letters = b"ARNDCQEGHILKMFPSTWYV";
    let base: Vec<u8> = (0..l).map(|_| letters[(r.next() % 20) as usize]).collect();
    let pairs = (0..n).map(|i| {
        let div = 0.02 + (max_div - 0.02) * i as f64 / (n.max(2) - 1) as f64;
        let s: Vec<u8> = base
            .iter()
            .map(|&c| if r.f() < div { letters[(r.next() % 20) as usize] } else { c })
            .collect();
        (format!("s{i}"), String::from_utf8(s).unwrap())
    });
    MSA::<Protein>::try_from_pairs(pairs).unwrap()
}

/// Best of `reps` timings of `run`, each on a fresh input from `setup` (not timed).
fn best<I, O>(reps: usize, mut setup: impl FnMut() -> I, mut run: impl FnMut(I) -> O) -> (f64, O) {
    let mut out = None;
    let t = (0..reps)
        .map(|_| {
            let input = setup();
            let t = Instant::now();
            out = Some(run(input));
            t.elapsed().as_secs_f64()
        })
        .fold(f64::MAX, f64::min);
    (t, out.unwrap())
}

fn main() {
    let a: Vec<usize> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    if a.len() < 2 {
        eprintln!("usage: ml_probe N L [reps]");
        std::process::exit(2);
    }
    let (n, l, reps) = (a[0], a[1], a.get(2).copied().unwrap_or(3));
    let pairs = (n * (n - 1) / 2) as f64;
    let _ = Random20::model(); // build the eigensystem outside the timing
    for max_div in [0.3, 0.7] {
        let (t_ml, d) = best(reps, || alignment(n, l, max_div), |m| m.into_dist::<ProteinMl<Random20>>(RateHet::NONE));
        let (t_poisson, _) = best(reps, || alignment(n, l, max_div), |m| m.into_dist::<Poisson>(RateHet::NONE));
        let mut inf = 0usize;
        for i in 0..n {
            for j in 0..i {
                inf += !d.get(i, j).is_finite() as usize;
            }
        }
        println!(
            "protein n={n} L={l} divergence<= {max_div:.1}: ML {:8.1} ms ({:5.2} us/pair) | Poisson {:6.1} ms | saturated pairs {inf}",
            t_ml * 1e3,
            t_ml / pairs * 1e6,
            t_poisson * 1e3
        );
    }
}
