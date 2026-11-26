# Changelog

All notable changes to this project are documented here. The format is loosely
based on [Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/). The three published artifacts —
the `nj` crate (crates.io), `nj_py` (PyPI), and `@holmrenser/nj` (npm) — share a
single version, bumped together via `make bump-{patch,minor,major}`.

## [Unreleased]

### Added
- `max_distance` option (default 10.0 substitutions per site) on `NJConfig`, the
  CLI (`--max-distance`) and the Python `nj()` function: the distance given to a
  sequence pair whose distance is not finite because the sequences are saturated
  for the model.

### Fixed
- Saturated sequence pairs no longer produce `NaN`/`inf` branch lengths. With, for
  example, `JukesCantor` and two sequences differing at every site, the distance is
  infinite; the NJ update then turned the infinite row sums into `NaN`, and the
  Newick output contained `:NaN` and `:inf` (and the JSON distance matrix `null`).
  Non-finite distances are now replaced by `max_distance` before the tree is built,
  a warning event reports how many pairs were affected, and the run continues.
  Bootstrap replicates are handled the same way (a single warning reports how many
  replicates were affected).

### Changed
- `nj()` now returns the capped distance matrix and average when
  `return_distance_matrix` / `return_average_distance` are set, i.e. the values the
  tree was built from. The standalone `distance_matrix()` and `average_distance()`
  are unchanged and still return raw estimates, including `inf`.
- `DistMat::neighbor_joining` returns an error for a matrix containing a non-finite
  distance instead of silently producing a `NaN` tree; use the new
  `DistMat::cap_non_finite` first.
- **Unknown characters are missing data.** `N` (DNA) and `X` (protein) are now
  excluded per pair, exactly like gaps, in every model. Before, they counted as a
  distinct state: one `N` against a base was a mismatch and `N` against `N` a match,
  so a pair differing only by an `N` had a non-zero distance. Alignments containing
  `N` or `X` now give smaller distances. Ambiguity codes were already read as `N`.
- DNA distances are computed from bit-planes (four bitmasks per sequence, A/C/G/T) and a
  pair-count table: `Kimura2P`, `TajimaNei` and `Tamura` are about 30× faster, `PDiff`
  and `JukesCantor` about 6×, with identical results. The parallel pair loop writes rows
  in blocks straight into the matrix instead of materialising every `(i, j)` pair and result.
- Neighbor-joining is about 30× faster: the active rows of the matrix are kept
  contiguous instead of testing an active-set bit and rebuilding the triangular
  index for every element (n=1000: 1.3 s to 44 ms; n=2000: 10.3 s to 0.33 s). The
  Q criterion is now computed as `(m−2)·d − (S_i + S_j)`, which is symmetric in the
  two nodes; the value can differ from `(m−2)·d − S_i − S_j` in the last bit, so
  which of two pairs that tie to within rounding is joined first may change.

## [0.0.26] - 2026-11-20

### Fixed
- Bootstrap support values are now independent of tree rooting. Splits were
  previously keyed by the clade below each node, so the same split could be
  tallied under different keys when replicate trees were rooted differently, and
  the root's two children double-counted its split. Splits are now keyed by the
  side excluding the first taxon and counted once per replicate. Trivial splits
  (one side with a single taxon) no longer receive a support label.

### Changed
- Internal: `lib.rs` split into `bootstrap`, `validate` and `runner` modules, with
  no public API change. README and Python type hints/docstrings updated to cover
  all seven substitution models and rate-variation options.
- Releases now rotate the `[Unreleased]` changelog section automatically and use
  it as the GitHub release notes.

## [0.0.25] - 2026-06-04

### Added
- Among-site rate-variation corrections for the distance models: a gamma
  rate-heterogeneity shape parameter (`gamma_shape`, Jin & Nei 1990) and a
  proportion of invariant sites (`p_invar`). For distance methods both are exact
  closed-form corrections to each model's `−ln(...)` step, so no discrete rate
  categories are needed; they apply to every correctable model (`JukesCantor`,
  `Kimura2P`, `TajimaNei`, `Tamura`, `Poisson`, `KimuraProtein`) and compose as
  `+I+Γ`. The gamma correction converges back to the uncorrected model as the
  shape parameter grows; `PDiff` (raw p-distance) is unaffected. Exposed across
  the CLI (`-g/--gamma-shape`, `-i/--p-invar`), Python, and WASM bindings, with
  `NJError::InvalidGammaShape` / `NJError::InvalidPInvar` for out-of-range values.

## [0.0.24] - 2026-06-04

### Added
- Three substitution models: `TajimaNei` (Tajima-Nei 1984, DNA — corrects
  Jukes-Cantor for unequal base frequencies), `Tamura` (Tamura 1992, DNA —
  Kimura two-parameter with a GC-content correction), and `KimuraProtein`
  (Kimura 1983, protein — empirical correction to the Poisson distance). All are
  closed-form and reduce to the model they extend under the appropriate
  conditions (Tajima-Nei → Jukes-Cantor at equal frequencies, Tamura →
  Kimura-2P at GC = 0.5). Exposed across the CLI, Python, and WASM bindings.

## [0.0.23] - 2026-06-04

### Fixed
- Protein sequences are now encoded case-insensitively; lowercase residues were
  previously mapped to `X`, silently corrupting protein alignments.
- Newick output now single-quotes identifiers containing whitespace or reserved
  characters (`()[]{}',:;`), producing valid Newick for names with spaces, etc.
- Bootstrap replicate failures now propagate as errors instead of panicking a
  worker thread.

### Added
- New `simd` Cargo feature providing an explicit `core::simd` distance kernel
  (~5-7x faster distance computation). It requires a nightly toolchain
  (unstable `portable_simd`); published wheels enable it, while crate/source
  builds fall back to the autovectorized scalar path on stable.
- Duplicate sequence identifiers are rejected with `NJError::DuplicateIdentifier`
  (they previously corrupted bootstrap clade counting).
- A `Warning` log event is emitted for all-gap sequences (whose distances are 0).
- Language bindings surface error kind: Python raises `ValueError` for usage
  errors and `RuntimeError` for internal failures; WASM throws an `Error` whose
  `name` is the stable error code. Failing user callbacks are logged rather than
  silently dropped.
- Criterion benchmarks (`make bench`) and proptest property + golden-tree tests.
- `scripts/check-versions.sh` and a CI job that verify the three package
  versions stay in sync.
- `SequenceObject::is_empty()` companion to `len()`.

### Changed
- The `cli` feature now implies `parallel`, so the installed `nj` binary runs
  bootstrap and distance computation multi-threaded by default; the bare library
  stays Rayon-free unless `parallel` is requested. Cap workers with `-t/--num-threads`.
- Python wheels are built once per platform as a stable-ABI (`abi3-py310`) wheel
  covering CPython 3.10+, replacing the per-minor-version build matrix.

## [0.0.22] and earlier

See the Git tag history (`git tag`) and the GitHub releases page for prior
versions.
