# Changelog

All notable changes to this project are documented here. The format is loosely
based on [Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/). The three published artifacts —
the `nj` crate (crates.io), `nj_py` (PyPI), and `@holmrenser/nj` (npm) — share a
single version, bumped together via `make bump-{patch,minor,major}`.

## [Unreleased]

### Fixed
- Protein sequences are now encoded case-insensitively; lowercase residues were
  previously mapped to `X`, silently corrupting protein alignments.
- Newick output now single-quotes identifiers containing whitespace or reserved
  characters (`()[]{}',:;`), producing valid Newick for names with spaces, etc.
- Bootstrap replicate failures now propagate as errors instead of panicking a
  worker thread.

### Added
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

### Changed
- Python wheels are built once per platform as a stable-ABI (`abi3-py310`) wheel
  covering CPython 3.10+, replacing the per-minor-version build matrix.

## [0.0.22] and earlier

See the Git tag history (`git tag`) and the GitHub releases page for prior
versions.
