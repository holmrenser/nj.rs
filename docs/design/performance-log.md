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
| 2b | `Prepared` alignment, bit-planes for DNA, tiled parallel pair loop | done (scalar popcount) |
| 2c | Bootstrap: fuse column gather with packing | not pursued (see below) |
| 3a | `fearless_simd` on stable: port the byte kernel, remove nightly and the `simd` gate | done |
| 3b | SIMD NJ row-minimum | done |
| 3c | SIMD popcount loop over the planes | done |
| 4 | ML solver SIMD | unblocked: the scalar solver exists (`ml.rs`, 4.5–9.5 µs per protein pair, `ml_probe`); SIMD per `performance.md` §2 not started |

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
- **Byte kernel uses lane counters (step 3a).** Compared in the prototype crate first: accumulating the lane
  masks into `i8` counters (flushed every 127 steps) took 0.027 ns/site against 0.038 for extracting a bitmask
  each step (and 0.084 for the nightly `core::simd` kernel it replaces), so the counter variant of decision D5
  was adopted.
- **`simd` Cargo feature kept as a deprecated no-op** so `--features simd` keeps building. The `simd-nightly` CI job
  and the pinned nightly in `release.yml` are gone; wheels build on stable. `rust-version = "1.89"` (the
  `fearless_simd` MSRV; the crate already needed 1.88 for let-chains).
- **Kernel tests run every level** (`kernels::testing::levels`): best, scalar fallback and, on x86, SSE2/SSE4.2/AVX2/AVX-512
  where the CPU offers them. Checked natively on NEON and, cross-compiled, under `qemu-x86_64` with CPU models
  Nehalem, Haswell and `max` (correctness only; AVX-512 was not offered by the emulator).
- **Plane popcount keeps byte counters (step 3c).** First version: one `u64` vector accumulator per cell, `count_ones`
  on `u64` lanes (on NEON that is `cnt` plus three widening adds per cell and step). Second version: popcount the
  words as bytes, accumulate in byte counters for up to 31 steps (at most 8 per step, 248 < 256), fold to 16-bit
  lanes and reduce once per group. The second is faster everywhere measured (n=1000 `k2p` distances 23.7 → 19.8 ms,
  L=10000: 183 → 123 ms) and was kept. The Harley–Seal variant for AVX2 was not written: there is no x86
  hardware here to time it, and the byte-counter loop already avoids the per-step widening that variant targets.
- **2c not pursued.** Measured per replicate (`K2P`, no NJ): gather is 1–10% and packing 2–18% of the distance work
  (smallest at n=1000, L=1000: 1.2% + 2.4%; largest at n=200, L=10000: 8% + 18%). A fused gather-and-pack still does
  one random read and one bit write per site, so it could save only part of that, roughly 1–5% of a replicate. A
  faster `pack` kernel (compare + movemask per plane) is the better lever if bootstrap on long, narrow alignments
  matters; left as a follow-up.
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

## After step 2b (bit-planes, tiled pair loop; scalar popcount)

| Case | Distances before | Distances after | `nj()` total before → after |
|---|---|---|---|
| DNA n=1000 `pdiff` | 153.9 ms | 20.9 ms (7×) | 1438 → 64.7 ms |
| DNA n=1000 `jc` | 155.7 ms | 25.3 ms (6×) | 1438 → 69.4 ms |
| DNA n=1000 `k2p` | 802.3 ms | 27.0 ms (30×) | 2085 → 70.8 ms |
| DNA n=1000 `tn` | 940.6 ms | 29.7 ms (32×) | 2224 → 73.6 ms |
| DNA n=1000 `tamura` | 837.9 ms | 26.3 ms (32×) | 2121 → 70.2 ms |
| DNA n=2000 `k2p` | 3177.6 ms | 105.6 ms (30×) | 13482 → 444.5 ms |
| protein n=1000 `poisson` | 78.7 ms | 65.7 ms (unchanged path) | 1364 → 109.6 ms |
| DNA n=200 `k2p`, 20 bootstrap | — | — | 0.91 s → 0.04 s |

The plane kernel is exact: every pair-count table equals the histogram, and the distance matrix equals a plain
per-pair loop bit for bit (block sizes 1–50 around the 16-row block, lengths around the 64-column word boundary,
both alphabets, with and without rate variation). Protein distances still run through the byte kernel.

## After step 3a (byte kernel on `fearless_simd`, stable)

| Case | Distances before | Distances after |
|---|---|---|
| protein n=1000 L=500 `poisson` | 65.7 ms (scalar byte kernel) | 10.3 ms (6.4×) |

## After step 3b (SIMD row minimum)

| Case | NJ stage, scalar rows | NJ stage, SIMD rows |
|---|---|---|
| DNA n=500 | 5.8 ms | 2.8 ms |
| DNA n=1000 | 44.1 ms | 20.5 ms (2.1×; 63× vs baseline) |
| DNA n=2000 | 329 ms | 157 ms (2.1×; 65× vs baseline) |

The whole row loop is one dispatched kernel (`kernels::row_minima`), so dispatch costs once per NJ iteration, not
per row. It is bit-identical to the scalar reference at every level, and the NJ oracle tests (bit-identical trees)
also pass on x86 SSE4.2 and AVX2 under emulation.

## After step 3c (SIMD popcount over the planes)

| Case | Distances, scalar popcount | Distances, SIMD popcount |
|---|---|---|
| DNA n=1000 L=1000 `pdiff` | 20.9 ms | 15.5 ms |
| DNA n=1000 L=1000 `k2p` | 27.0 ms | 19.8 ms |
| DNA n=1000 L=1000 `tn` | 29.7 ms | 22.9 ms |
| DNA n=2000 L=1000 `k2p` | 105.6 ms | 77.4 ms |
| DNA n=1000 L=10000 `k2p` | not measured | 122.9 ms |

## Final results (`9d5fd54`, single thread)

| Case | Baseline `nj()` | Final `nj()` | Baseline distances | Final distances | Baseline NJ | Final NJ |
|---|---|---|---|---|---|---|
| DNA n=500 `pdiff` | 199 ms | 6.9 ms (29×) | 38.7 ms | 4.1 ms | 160 ms | 2.8 ms |
| DNA n=1000 `pdiff` | 1438 ms | 36.4 ms (40×) | 153.9 ms | 15.6 ms | 1284 ms | 20.9 ms |
| DNA n=1000 `jc` | 1437 ms | 39.0 ms (37×) | 155.7 ms | 18.8 ms | 1282 ms | 20.2 ms |
| DNA n=1000 `k2p` | 2085 ms | 40.0 ms (52×) | 802.3 ms | 19.7 ms | 1283 ms | 20.3 ms |
| DNA n=1000 `tn` | 2224 ms | 43.1 ms (52×) | 940.6 ms | 22.9 ms | 1284 ms | 20.2 ms |
| DNA n=1000 `tamura` | 2121 ms | 40.4 ms (52×) | 837.9 ms | 20.3 ms | 1283 ms | 20.2 ms |
| DNA n=2000 `k2p` | 13482 ms | 233 ms (58×) | 3177.6 ms | 77.3 ms | 10304 ms | 155.8 ms |
| protein n=1000 `poisson` | 1364 ms | 31.0 ms (44×) | 78.7 ms | 10.4 ms | 1286 ms | 20.6 ms |
| protein n=500 `kimura-protein` | 173 ms | 4.9 ms (35×) | 12.2 ms | 2.1 ms | 161 ms | 2.9 ms |
| DNA n=200 `k2p`, 20 bootstrap | 0.91 s | 0.03 s (30×) | | | | |

With the `parallel` feature on 12 cores the distance stage shrinks further (n=2000 `k2p`: 77 → 18 ms, n=1000: 20 → 5 ms),
so NJ, which is serial, is then about 90% of a single `nj()` call; bootstrap replicates still run in parallel.

Verification: 183 tests with and without `parallel` and `cli` (plus 16 CLI tests); Python extension built and its 63
tests pass; the workspace (incl. benches, WASM with `+simd128`) compiles; kernels, NJ oracle and distance-matrix tests also
pass on x86-64 under `qemu-x86_64` (Nehalem → SSE4.2, Haswell → AVX2; the emulator did not offer AVX-512). Not measured:
real x86 and WASM timings.

## Follow-ups (not done)

- ML solver SIMD (step 4): the `X_ab[m]` precompute already exists; padding `K` to the vector width and two accumulators remain (`performance.md` §2.4). The scalar solver costs 4.5–9.5 µs per protein pair, so n = 1000 takes several seconds.
- A SIMD `pack` kernel for bit-plane construction (matters for bootstrap on long, narrow alignments).
- Harley–Seal accumulation on AVX2 (needs x86 hardware to evaluate); AVX-512 correctness.
- Pruned NJ search (RapidNJ/NINJA style) for n ≫ 10⁴; NJ is the serial part in multi-threaded runs.
- Warn when a sequence consists only of `N`/`X` (distance 0 against everything, silently).

## Commits

- `perf: add progress log and perf_probe example` — this file, `nj/examples/perf_probe.rs`, baseline above.
- `perf(nj): compact the matrix, scan rows as slices` — `nj/src/nj.rs` rewritten, previous implementation kept as
  `nj/src/nj_reference.rs` (`cfg(test)`), four oracle tests.
- `perf(models): treat N/X as missing; DNA models consume a 4×4 pair-count table` — `nj/src/counts.rs`,
  `AlphabetEncoding::N_STATES`, `from_counts` on `Kimura2P`/`TajimaNei`/`Tamura`, policy tests, CHANGELOG.
- `perf(distances): bit-plane pair counts and a tiled pair loop` — `nj/src/prepared.rs`, `ModelCalculation::distance_prepared`,
  `DistMat::from_prepared` (row blocks written in place; no `(i, j)` / result vectors), bit-exactness tests.
- `perf(simd): port the byte kernel to fearless_simd and drop nightly` — `nj/src/kernels.rs`, `fearless_simd` dependency,
  `simd` feature is a no-op, CI/release on stable, `rust-version = "1.89"`.
- `perf(nj): SIMD row-minimum kernel` — `kernels::row_minima` (+ scalar reference and per-level test), used by `select_min_q_pair`.
- `perf(distances): SIMD popcount over the bit-planes` — `kernels::plane_counts` (byte counters), `Prepared` stores the SIMD level.
