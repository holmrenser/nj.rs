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

`nj::nj(NJConfig) -> Result<String, String>` in `lib.rs` is the single public entry point for all three distribution targets. It auto-detects the alphabet (`DNA` vs `Protein`) and dispatches to the correct generic `run_nj::<A, M>()` monomorphisation.

### Substitution models (`models.rs`)

Models implement `ModelCalculation<A: AlphabetEncoding>`. DNA-only: `PDiff`, `JukesCantor`, `Kimura2P`. Protein-only: `Poisson`. `PDiff` works for both alphabets. Model–alphabet compatibility is enforced at runtime inside `nj()`.

### Bootstrap support

Bootstrap runs N replicate NJ trees on column-resampled MSAs, tallies clade membership using `BitVec`-keyed `HashMap` counters (`count_clades`), then maps support values back onto the main tree's internal nodes (`add_bootstrap_to_tree`).

### CLI (`main.rs`)

Enabled only with `--features cli`. Parses FASTA via `clap`, builds a `NJConfig`, calls `nj()`, and prints the Newick string.
