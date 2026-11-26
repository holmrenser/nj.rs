# Performance work log

Branch `perf/nj-simd`, started from `main` at `018e188`. Implements [`performance.md`](performance.md). Each
step is a small commit with the tests green; this file records what was done, what deviated from the plan, and
the measurements taken along the way.

Measurements: `cargo run --release -p nj --example perf_probe -- dna 1000 1000 k2p`, single thread, synthetic data
(12% DNA / 15% protein divergence), best of 3, on aarch64 NEON (Apple-class single-core bandwidth). x86 and WASM
execution were not measured. Timings guide the work; they are not pass/fail gates.

## Plan and status

| # | Step | Status |
|---|---|---|
| 0 | Log, probe example, baseline | done |
| 1 | NJ: contiguous rows, symmetric `Q`, oracle under `cfg(test)` | pending |
| 2a | Missing data: `N`/`X` dropped per pair; DNA models consume a 4×4 pair-count table | pending |
| 2b | `Prepared` alignment, bit-planes for DNA, tiled parallel pair loop | pending |
| 2c | Bootstrap: fuse column gather with packing | pending |
| 3a | `fearless_simd` on stable: port the byte kernel, remove nightly and the `simd` gate | pending |
| 3b | SIMD NJ row-minimum | pending |
| 3c | SIMD popcount loop over the planes | pending |
| 4 | ML solver SIMD | blocked: the ML solver does not exist yet (`substitution-models.md` plan step 4) |

## Decisions and deviations

(filled in as the work proceeds)

## Baseline (`018e188`)

| Case | Distances | `nj()` total | NJ stage |
|---|---|---|---|
| DNA n=500 L=1000 `pdiff` | 38.7 ms | 199.0 ms | 160.3 ms |
| DNA n=1000 L=1000 `pdiff` | 153.9 ms | 1438.3 ms | 1284.4 ms |
| DNA n=1000 L=1000 `jc` | 155.7 ms | 1437.5 ms | 1281.8 ms |
| DNA n=1000 L=1000 `k2p` | 802.3 ms | 2085.4 ms | 1283.1 ms |
| DNA n=1000 L=1000 `tn` | 940.6 ms | 2224.3 ms | 1283.7 ms |
| DNA n=1000 L=1000 `tamura` | 837.9 ms | 2120.9 ms | 1283.0 ms |
| DNA n=2000 L=1000 `k2p` | 3177.6 ms | 13481.8 ms | 10304.2 ms |
| protein n=1000 L=500 `poisson` | 78.7 ms | 1364.2 ms | 1285.5 ms |
| protein n=500 L=300 `kimura-protein` | 12.2 ms | 173.3 ms | 161.2 ms |
| DNA n=200 L=1000 `k2p`, 20 bootstrap | 32.6 ms | 43.1 ms | 10.5 ms (bootstrap run: 0.91 s) |

(Stable toolchain without the `simd` feature, i.e. the autovectorised scalar byte kernel.)

## Commits

- `perf: add progress log and perf_probe example` — this file, `nj/examples/perf_probe.rs`, baseline above.
