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
| 1 | NJ: contiguous rows, symmetric `Q`, oracle under `cfg(test)` | done |
| 2a | Missing data: `N`/`X` dropped per pair; DNA models consume a 4×4 pair-count table | done |
| 2b | `Prepared` alignment, bit-planes for DNA, tiled parallel pair loop | pending |
| 2c | Bootstrap: fuse column gather with packing | pending |
| 3a | `fearless_simd` on stable: port the byte kernel, remove nightly and the `simd` gate | pending |
| 3b | SIMD NJ row-minimum | pending |
| 3c | SIMD popcount loop over the planes | pending |
| 4 | ML solver SIMD | blocked: the ML solver does not exist yet (`substitution-models.md` plan step 4) |

## Decisions and deviations

- **Symmetric `Q` (step 1).** The old `(m−2)·d − S_i − S_j` rounds differently depending on which node is `i`
  (the larger original index). After compaction the slot order no longer follows original order, so that
  association cannot be reproduced by a vectorisable loop. `Q` is now `(m−2)·d − (S_i + S_j)`, which is symmetric
  and layout-independent. The reference implementation uses the same definition, so the oracle comparison stays
  bit-exact. This can change which of two *near-equal* pairs is joined, relative to releases before this branch,
  in the last bit only; the golden and property tests are unaffected.
- **Missing data (step 2a).** `AlphabetEncoding::N_STATES` (4/20) replaces the gap byte in the kernel: a byte
  `>= N_STATES` is missing, which covers the gap and `N`/`X` with one compare. `K2P`, `TajimaNei` and `Tamura`
  became functions of a 4×4 `DnaCounts` table (`from_counts`); `PDiff`/`JC` keep the byte kernel until the planes
  arrive. `TajimaNei` therefore works on 4 states instead of 5 (`N` was a state). This is the accepted decision
  "treat `N`/`X` as missing in the existing models too"; it changes results for alignments containing `N`/`X`
  (CHANGELOG entry added). A sequence of only `N`s still gives distance 0 with no warning (the existing warning
  covers only all-gap sequences).
- **WASM already builds with `+simd128`** (`wasm/.cargo/config.toml`), so the "WASM needs a second SIMD build"
  caveat in `performance.md` does not apply to this repo.

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

## After step 1 (NJ compaction)

| Case | NJ stage before | NJ stage after |
|---|---|---|
| DNA n=500 | 160.3 ms | 5.8 ms (28×) |
| DNA n=1000 | 1284 ms | 44.1 ms (29×) |
| DNA n=2000 | 10304 ms | 329 ms (31×) |
| protein n=1000 | 1285 ms | 43.9 ms (29×) |

A scalar, auto-vectorised row minimum (`.min()` on four lanes) is enough for this; the explicit SIMD kernel
(step 3b) adds on top. Oracle tests: random, integer-tie, duplicate-taxon and degenerate matrices give
bit-identical trees (identifiers, branch-length bits, labels).

## Commits

- `perf: add progress log and perf_probe example` — this file, `nj/examples/perf_probe.rs`, baseline above.
- `perf(nj): compact the matrix, scan rows as slices` — `nj/src/nj.rs` rewritten, previous implementation kept as
  `nj/src/nj_reference.rs` (`cfg(test)`), four oracle tests.
- `perf(models): treat N/X as missing; DNA models consume a 4×4 pair-count table` — `nj/src/counts.rs`,
  `AlphabetEncoding::N_STATES`, `from_counts` on `Kimura2P`/`TajimaNei`/`Tamura`, policy tests, CHANGELOG.
