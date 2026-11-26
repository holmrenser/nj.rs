# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Build

```bash
# Build all targets (library, CLI, Python bindings, WASM)
make all

# Build only the Rust library and CLI
make nj
# Equivalent: cd nj && CARGO_TARGET_DIR=target cargo build --release && cargo build --release --features cli --bin nj

# Build Python bindings (requires maturin and uv)
make python

# Build WASM bindings
make wasm
```

### Tests

```bash
# Run all tests (Rust + WASM + Python)
make test

# Rust unit tests only
cargo test

# Run a single Rust test by name
cargo test test_nj_simple

# Python binding tests (requires maturin develop first)
make -C python test   # uv run pytest

# WASM binding tests (requires wasm build first)
make -C wasm test     # both node and web targets
```

### Versioning

```bash
# Bump patch version across all workspace members and sync to pyproject.toml and package.json
make bump
```

## Architecture

This is a Cargo workspace with three members:

- **`nj/`** — Core library + CLI binary. The `cli` feature flag gates the `clap` dependency; the default build produces only the library.
- **`python/`** — PyO3 bindings; thin wrapper that calls `nj::nj()` and exposes it to Python as `nj_py`.
- **`wasm/`** — wasm-bindgen bindings; thin wrapper targeting JavaScript/TypeScript as `@holmrenser/nj`.

### Data flow (nj crate)

```
FASTA input / SequenceObject list
        │
        ▼
  MSA<A>  (msa.rs)          — typed by alphabet A (DNA or Protein)
        │  .into_dist::<M>()
        ▼
  DistMat  (dist.rs)        — compact lower-triangular flat Vec, O(n²/2) memory
        │  .neighbor_joining()
        ▼
  NJState::run()  (nj.rs)   — iterative NJ with BitVec active-set and incremental row sums
        │
        ▼
  TreeNode  (tree.rs)       — recursive binary tree; .to_newick() produces output
```

### Key types

| Type | File | Purpose |
|---|---|---|
| `MSA<A>` | `msa.rs` | Multiple sequence alignment parameterised by `AlphabetEncoding` |
| `DistMat` | `dist.rs` | Lower-triangular distance matrix; `from_msa::<M, A>()` computes pairwise distances |
| `NJState` | `nj.rs` | Stateful NJ runner; holds `DistMat`, active `BitVec`, precomputed row sums |
| `TreeNode` | `tree.rs` | Recursive binary tree node; leaf vs. internal distinguished by `children: Option<[Box<TreeNode>; 2]>` |
| `NJConfig` | `config.rs` | Public API config struct consumed by `nj()` |

### Public API entry point

`nj::nj(NJConfig, Option<Box<dyn Fn(NJEvent)>>) -> Result<NJResult, NJError>` in `lib.rs` is the main entry point for all three distribution targets; `NJResult` holds the Newick string plus the optional distance matrix and average distance. `distance_matrix(DistConfig)` and `average_distance(DistConfig)` cover distance-only computation. `nj()` auto-detects the alphabet (`DNA` vs `Protein`) and dispatches to the correct generic `run_nj::<A, M>()` monomorphisation.

`lib.rs` holds only the public functions, re-exports and tests; the implementation lives in private modules:

| Module | Purpose |
|---|---|
| `validate.rs` | `validate_msa`, `validate_rate_params`, `detect_alphabet`, all-gap warning |
| `runner.rs` | `build_thread_pool`, `run_nj` / `run_distance_matrix` / `run_average_distance`, and the `dispatch_run!` macro |
| `bootstrap.rs` | Replicate runs (sequential and Rayon-parallel), split counting, support annotation |

### Substitution models (`models.rs`)

Models implement `ModelCalculation<A: AlphabetEncoding>`. DNA-only: `JukesCantor`, `Kimura2P`, `TajimaNei`, `Tamura`. Protein-only: `Poisson`, `KimuraProtein`. `PDiff` works for both alphabets. Model–alphabet compatibility is enforced at runtime inside `nj()` (the `dispatch_run!` macro in `runner.rs`).

`distance()` takes a `RateHet { gamma_shape: Option<f64>, p_invar: f64 }` value (threaded through `into_dist`/`from_msa`/`pairwise_distance_with`) that applies among-site rate variation to each model's `−c·ln(arg)` term: gamma rate heterogeneity becomes `c·α·(arg^(−1/α) − 1)` via the shared `corrected_term` helper, and invariant sites divide observed proportions by `(1 − p_invar)` then scale the result by `(1 − p_invar)`. `RateHet::NONE` reproduces the classic formulas exactly; `PDiff` ignores it. The `gamma_shape`/`p_invar` config fields are validated in `validate_rate_params` (`NJError::InvalidGammaShape` / `InvalidPInvar`).

### Bootstrap support

Bootstrap runs N replicate NJ trees on column-resampled MSAs, tallies splits using `BitVec`-derived keys in a `HashMap` counter, then maps support values back onto the main tree's internal nodes (`add_bootstrap_to_tree`).

NJ output is unrooted but stored with an arbitrary root, so splits are keyed by `split_key` — the side that excludes taxon 0 — making keys independent of rooting. `collect_splits` records each split once per tree (the root's two children describe the same split), and trivial splits (a side with fewer than two taxa) are never labelled.

### Performance structure

- **Missing data.** A symbol byte `>= AlphabetEncoding::N_STATES` (4 DNA, 20 protein) is missing: gaps and `N`/`X` are excluded per pair in every model (pairwise deletion).
- **`prepared.rs` / `counts.rs`.** `DistMat::from_msa` builds a `Prepared` alignment; for DNA it holds four bit-planes (A, C, G, T) per sequence, and `dna_counts(i, j)` returns the exact 4×4 pair-count table (`DnaCounts`) as `Σ popcount(plane_a(i) & plane_b(j))`. The DNA models are functions of that table (`from_counts`); `ModelCalculation::distance` (symbol slices, histogram) is the reference path, `distance_prepared` the fast one. Protein uses the byte kernel.
- **`DistMat::from_prepared`** fills blocks of 16 rows in place (Rayon over blocks with `parallel`).
- **`kernels.rs`.** All SIMD kernels (`fearless_simd`, stable Rust, run-time dispatch; scalar fallback) next to their scalar references: the `(n_diff, n_comparable)` byte kernel, `row_minima` (NJ search) and `plane_counts`. Each is `#[simd]`-annotated (without it the code loses its inlining) and tested at every level via `kernels::testing::levels()`.
- **`nj.rs`.** Active nodes occupy slots `0..m` of the packed matrix (swap-remove on join). `Q = (m−2)·d − (S_i + S_j)`, ties go to the smallest original `(i, j)`, and row sums update in ascending original index; `nj_reference.rs` (`cfg(test)`) is the straightforward oracle and the tests require bit-identical trees.
- Design and results: `docs/design/performance.md`, `docs/design/performance-log.md`.

### CLI (`main.rs`)

Enabled only with `--features cli`. Parses FASTA via `clap`, builds a `NJConfig`, calls `nj()`, and prints the Newick string (or a JSON object when `--distance-matrix` / `--average-distance` is set).
