# Design: performance priorities

Status: proposal. Tells us what to do first and why. Companion to [`substitution-models.md`](substitution-models.md).

## 1. Priorities

| # | Change | Why this slot | Expected effect |
|---|---|---|---|
| 1 | NJ: contiguous rows | NJ is the largest cost (≈96% of a `PDiff` run at n=1000, ≈60% for `K2P` at L=1000); its loop is bound by control flow, not arithmetic | ≈20× on the search |
| 2 | DNA pair counts: bit-slicing, `Prepared`, tiled pair loop | Counting is the next cost; the histogram cannot be vectorised, bit-slicing can | ≈7× on counting, all DNA models share it |
| 3 | `fearless_simd`: byte kernel, NJ row-min, popcount loop | SIMD only pays on the regular loops that 1 and 2 create | ≈2–3× on top, proportional to vector width |
| 4 | ML solver SIMD | Short vectors and per-entry reductions limit it | ≤2–3×, only for protein ML |

Steps 1 and 2 need no SIMD and no new dependency. Do them first, in scalar Rust.

## 2. What each step does and why

### 1. NJ (`nj.rs`, `distance_matrix.rs`)
Each iteration minimises `Q_ij = (m−2)·d_ij − S_i − S_j` over all active pairs, `n³/6` elements in total. Per element
that is one 8-byte load, a multiply, a subtract and a min. Today it costs 7.7 ns because each element also tests an
active-set bit, rebuilds the flat triangular index and branches. No vector width fixes that; contiguity does.

- Keep active rows contiguous (swap-remove the joined row); scan each row as a slice with a row-minimum reduction.
- Define `Q_ij = (m−2)·d_ij − (S_i + S_j)`. The old `(m−2)·d_ij − S_i − S_j` rounds differently depending on which
  node is `i`, and slot order no longer follows original order; addition is commutative, so the symmetric form
  depends only on the pair. The oracle uses the same definition.
- Keep multiply and subtract separate (no `mul_add`) so results are bit-identical at every instruction-set level.
- Keep today's tie-breaking (first minimum by original `(i, j)`): record per-row minima, rescan only rows equal to
  the global minimum and compare original ids, since slot order no longer follows them. This matters for
  alignments with duplicate sequences.
- Keep the row-sum updates bit-identical to today's. `update_distances` accumulates `S_i` and `S_k` one `k` at a
  time in ascending original index, so the rounding depends on that order. After a swap-remove the slot order
  differs, so visit `k` in ascending original id (through a sorted list of active ids, O(m) per iteration). Without
  this the `Q` values can differ in the last bit from the oracle, which changes which of several near-equal pairs
  is chosen, most visibly with duplicate sequences.
- Keep the current implementation under `cfg(test)` as the oracle. Keep f64.
- **Scope:** exact unpruned NJ for n ≲ 10⁴ (matrix `4·n²` bytes: 0.4 GB at 10⁴, 40 GB at 10⁵). Once the matrix
  outgrows the cache the loop is bandwidth-bound and SIMD adds nothing; the way beyond is pruning
  ([RapidNJ], [NINJA], [Clausen 2022]), which is out of scope.

### 2. DNA pair counts
The per-pair table `F_ab` is built today by a histogram: a read-modify-write per site at a data-dependent address,
which has no vector form. Store four bitmasks per sequence (A, C, G, T) instead; then
`F_ab = Σ_words popcount(mask1_a & mask2_b)`, 56 operations per 64 sites against ≈384 for the histogram.

- Gaps, `N` and ambiguity codes set no bit, so they leave `F` automatically (the missing-data decision).
- `PDiff`/`JC` come from the same planes: `comparable = popcount(valid1 & valid2)`, `diff = comparable − Σ_a F_aa`.
- DNA only. The cost grows as `K²`, so bit-slicing loses to the histogram above `K ≈ 10`; protein stays on bytes.
- Planes are built once per alignment and once per bootstrap replicate, with the column gather fused into the
  packing (`prepare_resampled`).
- Replace the parallel `from_msa`, which materialises every `(i, j)` (16 B) and every result (24 B), about 2 GB of
  transient memory at n=10⁴, by row-block tiles written straight into the packed triangle.

### 3. SIMD with `fearless_simd` (stable Rust)
A vector instruction does the same work on W lanes (`f64`: 2 on NEON/SSE, 4 on AVX2, 8 on AVX-512). The gain is
proportional to W, so one binary needs run-time selection of the widest level; `fearless_simd` provides that with
a safe API on stable and a scalar fallback, replacing nightly `core::simd`. Do this in order:

1. **Toolchain.** Port the byte kernel (`n_diff`, `n_comparable`; still used for protein) and remove the `simd`
   feature gate, the pinned nightly in `release.yml` and the `simd-nightly` CI job. MSRV becomes 1.89. In the port,
   accumulate lane masks into `u8` counters, flushed every 255 steps, instead of extracting a bitmask each step.
2. **NJ row-minimum** on the contiguous rows from step 1.
3. **Plane popcount loop:** one vector accumulator per `(a, b)` cell for the whole loop, reduced once. AVX2 has no
   vector popcount; use Harley–Seal carry-save accumulation there via `kernel!` ([Muła 2018]).

Every kernel needs the `#[simd]` attribute; without it the function is not inlined into the dispatched code and
silently loses its speed.

### 4. ML solver (after the solver exists)
Per Newton iteration each nonzero `F_ab` needs three length-`K` dot products (with `e`, `λe`, `λ²e`) sharing
`x_k = U_ak·W_kb`. The iterations are sequential, and per-entry horizontal sums cap the vector gain: for `K=20`
the model gives 1.7× at W=2, 2.8× at W=4, 3.3× at W=8.

- Precompute `X_ab[k] = U_ak·W_kb` once per model (64 KB for `K=20`): removes a quarter of the operations at any W.
- Pad `K` to a multiple of W; use at least two independent accumulators per entry.
- Start Newton from the closed-form distance: saves about one of four iterations.
- The DNA solver (16 entries) is already small next to counting; do not vectorise it.

## 3. Not doing

- Bit-slicing for protein (above the `K ≈ 10` crossover).
- f32 NJ matrix: halves the bytes but changes rounding and ties against the oracle.
- Site-pattern compression: popcount needs unweighted columns, and weights would multiply the cost by the number
  of weight bits. Revisit if redundant alignments dominate (the repo's 16S example has 978 unique of 2,444 columns).
- Vectorising across pairs; pruned NJ search (scope above).

## 4. Interface change

Packed planes cannot be reached from `distance(a: &[Symbol], b: &[Symbol])`. Models work on a prepared
alignment and a pair of indices:

```rust
pub struct Prepared<A: AlphabetEncoding> { /* byte symbols; for DNA also 4 bit-planes per sequence */ }
impl<A: AlphabetEncoding> MSA<A> {
    fn prepare(&self) -> Prepared<A>;
    fn prepare_resampled(&self, cols: &[usize]) -> Prepared<A>;       // bootstrap: gather + pack in one pass
}
impl<A: AlphabetEncoding> Prepared<A> {
    fn pair_counts(&self, i: usize, j: usize) -> PairCounts;          // exact integers; K×K
}
pub trait DistanceModel<A: AlphabetEncoding>: Send + Sync {
    fn distance(&self, data: &Prepared<A>, i: usize, j: usize) -> f64;   // INFINITY if saturated
}
```

`DistMat::from_msa` iterates tiles of index pairs and calls `distance(&prepared, i, j)`; `run_one_replicate` builds
the replicate with `prepare_resampled`. `substitution-models.md` §4 reflects this.

## 5. Validation

Correctness at every instruction-set level (scalar fallback, SSE2/SSE4.2/AVX2, NEON):

- pair counts equal the scalar histogram exactly, including columns with gaps and `N`;
- NJ output equals the retained current implementation on golden, property and random-matrix tests, including
  duplicate sequences;
- solver dot products agree with the scalar version to 1e-12 relative, and `K2P`/`TN`/`Tamura` distances are unchanged to 1e-12.

Prototypes on aarch64 NEON agreed with the expectations above (NJ 7.7 → 0.3–0.4 ns per element with contiguity, 0.11–0.14
with SIMD; bit-slicing 6–9× scalar, 13–14× NEON; solver 1.5–1.7×). x86 was not measured.

## References

- RapidNJ: Simonsen M, Mailund T, Pedersen CNS (2008). Rapid neighbour-joining. *WABI*, LNCS 5251:113–122.
- NINJA: Wheeler TJ (2009). Large-scale neighbor-joining with NINJA. *WABI*, LNCS 5724:375–389.
- Clausen PTLC (2022). Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining.
  *Bioinformatics* 39(1):btac774.
- Muła W, Kurz N, Lemire D (2018). Faster population counts using AVX2 instructions. *Computer Journal* 61(1):111–120.

[RapidNJ]: https://link.springer.com/doi/10.1007/978-3-540-87361-7_10
[NINJA]: https://unpaywall.org/10.1007%2F978-3-642-04241-6_31
[Clausen 2022]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9805563/
[Muła 2018]: https://arxiv.org/abs/1611.07612
