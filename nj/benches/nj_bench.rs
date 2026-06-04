//! Benchmarks for the core Neighbor-Joining pipeline.
//!
//! Covers the three public entry points over synthetic DNA alignments of
//! increasing size: distance-matrix construction, a full NJ run, and a bootstrap
//! NJ run. Run with `cargo bench -p nj` (single-threaded) or
//! `cargo bench -p nj --features parallel` to measure the Rayon-parallel paths.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use nj::models::SubstitutionModel;
use nj::{DistConfig, NJConfig, SequenceObject, distance_matrix, nj as run_nj};
use std::hint::black_box;

/// Tiny deterministic xorshift PRNG so benchmarks are reproducible and pull in
/// no extra dependencies.
struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

/// Builds `n` random DNA sequences of length `len` with a fixed seed.
fn synthetic_dna(n: usize, len: usize) -> Vec<SequenceObject> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut rng = XorShift(0x9E3779B97F4A7C15);
    (0..n)
        .map(|i| {
            let sequence: String = (0..len)
                .map(|_| BASES[(rng.next() % 4) as usize] as char)
                .collect();
            SequenceObject {
                identifier: format!("Seq{i}"),
                sequence,
            }
        })
        .collect()
}

fn dist_config(msa: Vec<SequenceObject>) -> DistConfig {
    DistConfig {
        msa,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
        gamma_shape: None,
        p_invar: None,
    }
}

fn nj_config(msa: Vec<SequenceObject>, n_bootstrap_samples: usize) -> NJConfig {
    NJConfig {
        msa,
        n_bootstrap_samples,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: false,
        return_average_distance: false,
        gamma_shape: None,
        p_invar: None,
    }
}

fn bench_distance_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("distance_matrix");
    for &n in &[25usize, 100, 250] {
        let msa = synthetic_dna(n, 500);
        group.bench_with_input(BenchmarkId::from_parameter(n), &msa, |b, msa| {
            b.iter(|| distance_matrix(black_box(dist_config(msa.clone()))).unwrap());
        });
    }
    group.finish();
}

/// Wide-alignment regime (large L, modest n) where the O(n²·L) distance kernel
/// dominates — the case explicit/auto SIMD on the kernel targets.
fn bench_distance_matrix_wide(c: &mut Criterion) {
    let mut group = c.benchmark_group("distance_matrix_wide");
    let msa = synthetic_dna(100, 5000);
    group.bench_with_input(BenchmarkId::from_parameter(100), &msa, |b, msa| {
        b.iter(|| distance_matrix(black_box(dist_config(msa.clone()))).unwrap());
    });
    group.finish();
}

fn bench_nj(c: &mut Criterion) {
    let mut group = c.benchmark_group("nj");
    for &n in &[25usize, 100, 250] {
        let msa = synthetic_dna(n, 500);
        group.bench_with_input(BenchmarkId::from_parameter(n), &msa, |b, msa| {
            b.iter(|| run_nj(black_box(nj_config(msa.clone(), 0)), None).unwrap());
        });
    }
    group.finish();
}

fn bench_bootstrap(c: &mut Criterion) {
    let mut group = c.benchmark_group("nj_bootstrap_100");
    for &n in &[25usize, 100] {
        let msa = synthetic_dna(n, 500);
        group.bench_with_input(BenchmarkId::from_parameter(n), &msa, |b, msa| {
            b.iter(|| run_nj(black_box(nj_config(msa.clone(), 100)), None).unwrap());
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_distance_matrix,
    bench_distance_matrix_wide,
    bench_nj,
    bench_bootstrap
);
criterion_main!(benches);
