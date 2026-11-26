//! Stage timings for the performance work (single thread, synthetic data).
//!
//! ```text
//! cargo run --release -p nj --example perf_probe -- dna 1000 1000 k2p [reps] [bootstrap]
//! ```
//!
//! Prints the distance stage (`distance_matrix`), the whole `nj` run without bootstrap, the NJ stage
//! (their difference) and, if requested, a bootstrap run. Models: pdiff, jc, k2p, tn, tamura (DNA);
//! pdiff, poisson, kimura-protein (protein). Best of `reps` (default 3).

use nj::models::SubstitutionModel as M;
use nj::{DistConfig, NJConfig, SequenceObject, distance_matrix, nj as run_nj};
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

fn synthetic(alphabet: &[u8], n: usize, l: usize, divergence: f64) -> Vec<SequenceObject> {
    let mut r = Rng(0x9E37_79B9_7F4A_7C15);
    let k = alphabet.len() as u64;
    let base: Vec<u8> = (0..l).map(|_| alphabet[(r.next() % k) as usize]).collect();
    (0..n)
        .map(|i| SequenceObject {
            identifier: format!("s{i}"),
            sequence: String::from_utf8(
                base.iter()
                    .map(|&c| if r.f() < divergence { alphabet[(r.next() % k) as usize] } else { c })
                    .collect(),
            )
            .unwrap(),
        })
        .collect()
}

fn best(reps: usize, mut f: impl FnMut()) -> f64 {
    (0..reps)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::MAX, f64::min)
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 4 {
        eprintln!("usage: perf_probe dna|protein N L MODEL [reps] [bootstrap]");
        std::process::exit(2);
    }
    let (n, l): (usize, usize) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let reps: usize = a.get(4).map_or(3, |s| s.parse().unwrap());
    let boot: usize = a.get(5).map_or(0, |s| s.parse().unwrap());
    let model = match a[3].as_str() {
        "pdiff" => M::PDiff,
        "jc" => M::JukesCantor,
        "k2p" => M::Kimura2P,
        "tn" => M::TajimaNei,
        "tamura" => M::Tamura,
        "poisson" => M::Poisson,
        "kimura-protein" => M::KimuraProtein,
        other => panic!("unknown model {other}"),
    };
    let msa = match a[0].as_str() {
        "dna" => synthetic(b"ACGT", n, l, 0.12),
        "protein" => synthetic(b"ARNDCQEGHILKMFPSTWYV", n, l, 0.15),
        other => panic!("unknown alphabet {other}"),
    };
    let t_dist = best(reps, || {
        let _ = distance_matrix(DistConfig {
            msa: msa.clone(),
            substitution_model: model.clone(),
            alphabet: None,
            num_threads: None,
            gamma_shape: None,
            p_invar: None,
        })
        .unwrap();
    });
    let conf = |b: usize| NJConfig {
        msa: msa.clone(),
        n_bootstrap_samples: b,
        substitution_model: model.clone(),
        alphabet: None,
        num_threads: None,
        return_distance_matrix: false,
        return_average_distance: false,
        max_distance: None,
        gamma_shape: None,
        p_invar: None,
    };
    let t_all = best(reps, || {
        let _ = run_nj(conf(0), None).unwrap();
    });
    print!(
        "{} n={n} L={l} {:<14} distances {:8.1} ms | nj() total {:8.1} ms | NJ stage {:8.1} ms",
        a[0],
        a[3],
        t_dist * 1e3,
        t_all * 1e3,
        (t_all - t_dist).max(0.0) * 1e3
    );
    if boot > 0 {
        let t_boot = best(1.max(reps / 3), || {
            let _ = run_nj(conf(boot), None).unwrap();
        });
        print!(" | nj() with {boot} bootstrap {:8.2} s", t_boot);
    }
    println!();
}
