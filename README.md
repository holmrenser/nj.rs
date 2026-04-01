# nj.rs
Neighbor-Joining phylogenetic tree inference in Rust, with Python and WASM bindings.

[![Release](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml)
[![Tests](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml)


[![Crates.io Version](https://img.shields.io/crates/v/nj?logo=rust&label=crates.io%3A%20nj)](https://crates.io/crates/nj)
[![PyPI - Version](https://img.shields.io/pypi/v/nj_py?logo=pypi&logoColor=blue&label=pypi%3A%20nj_py)](https://pypi.org/project/nj_py/)
[![NPM Version](https://img.shields.io/npm/v/%40holmrenser%2Fnj?style=flat&logo=npm&label=npm%3A%20%40holmrenser%2Fnj)](https://www.npmjs.com/package/@holmrenser/nj)

Takes an aligned FASTA file as input and outputs a Newick string. Alphabet (DNA/protein) is auto-detected. Supports optional bootstrap support values on internal nodes.

**Substitution models:** `PDiff` (DNA + protein), `JukesCantor` (DNA), `Kimura2P` (DNA), `Poisson` (protein)

## CLI

```bash
cargo install nj --features cli

nj sequences.fasta
nj --substitution-model kimura2-p --n-bootstrap-samples 100 sequences.fasta > tree.nwk
```

A progress bar is shown on stderr when bootstrapping.

## Rust

```toml
[dependencies]
nj = "0.0.12"
```

```rust
use nj::{NJConfig, SequenceObject, nj, models::SubstitutionModel};

let newick = nj(NJConfig {
    msa: vec![
        SequenceObject { identifier: "A".into(), sequence: "ACGT".into() },
        SequenceObject { identifier: "B".into(), sequence: "ACCT".into() },
    ],
    n_bootstrap_samples: 100,
    substitution_model: SubstitutionModel::JukesCantor,
}, Some(Box::new(|current, total| println!("{current}/{total}"))))?;
```

## Python

```bash
pip install nj_py
```

```python
from tqdm import tqdm
from nj_py import nj

config = {
    "msa": [
        {"identifier": "A", "sequence": "ACGT"},
        {"identifier": "B", "sequence": "ACCT"},
    ],
    "n_bootstrap_samples": 100,
    "substitution_model": "JukesCantor",
}

with tqdm(total=100) as progress_bar:
    newick = nj(config, on_progress=lambda current, total: progress_bar.update(1))
```

## JavaScript / WASM

```bash
npm install @holmrenser/nj
```

```js
import { nj } from '@holmrenser/nj';

const config = {
    msa: [
        { identifier: 'A', sequence: 'ACGT' },
        { identifier: 'B', sequence: 'ACCT' },
    ],
    n_bootstrap_samples: 100,
    substitution_model: 'JukesCantor',
};

const newick = nj(config, (current, total) => {
    progressBar.value = current / total * 100;
});
```
