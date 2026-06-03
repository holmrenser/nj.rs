# Contributing

Thanks for your interest in contributing to `nj.rs`! This is a Cargo workspace
that ships three artifacts from one core crate:

- `nj/` — core Rust library + CLI (the `cli` and `parallel` features are optional)
- `python/` — PyO3 bindings published to PyPI as `nj_py`
- `wasm/` — wasm-bindgen bindings published to npm as `@holmrenser/nj`

## Building & testing

```bash
make all      # build library, CLI, Python, and WASM
make test     # run Rust (incl. cli + parallel), WASM (node + web), and Python tests
make bench    # run the criterion benchmarks
```

You can also test pieces individually:

```bash
cargo test -p nj                     # core unit + integration tests
cargo test -p nj --features parallel # exercise the parallel bootstrap path
make -C python test                  # maturin develop + pytest (needs maturin + uv)
make -C wasm test                    # wasm-pack build + node/web tests (needs wasm-pack + node)
```

Please add tests with any behavioural change. Cross-cutting changes to the core
API should be reflected in the Python (`python/tests/`) and WASM
(`wasm/tests/`) suites where relevant.

## Documentation single source of truth

The root `README.md` is the **only** committed README. It is copied into
`python/` and `wasm/` at release time and is git-ignored there, so never commit
`python/README.md` or `wasm/README.md`.

## Releasing

All three packages share one version. Bump them together — never edit the
version strings by hand:

```bash
make bump-patch   # or bump-minor / bump-major
```

`release.sh` (invoked by `make bump-*`) requires `cargo-workspaces`, `toml`,
`jq`, and `sponge` on `PATH`. It bumps the Cargo workspace version, syncs
`python/pyproject.toml` and `wasm/package.json`, runs `scripts/check-versions.sh`
to confirm they agree, then creates a `Release vX.Y.Z` commit and tag.

Pushing a `vX.Y.Z` tag triggers `.github/workflows/release.yml`, which runs the
full test matrix, re-verifies the version, and publishes to crates.io, PyPI
(OIDC trusted publishing), and npm (OIDC + provenance).

You can sanity-check version consistency at any time:

```bash
bash scripts/check-versions.sh
```
