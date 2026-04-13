# nj.rs
Neighbor-Joining phylogenetic tree inference in Rust, with Python and WASM bindings.

[![Release](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml)
[![Tests](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml)


[![Crates.io Version](https://img.shields.io/crates/v/nj?logo=rust&label=crates.io%3A%20nj)](https://crates.io/crates/nj)
[![PyPI - Version](https://img.shields.io/pypi/v/nj_py?logo=pypi&logoColor=blue&label=pypi%3A%20nj_py)](https://pypi.org/project/nj_py/)
[![NPM Version](https://img.shields.io/npm/v/%40holmrenser%2Fnj?style=flat&logo=npm&label=npm%3A%20%40holmrenser%2Fnj)](https://www.npmjs.com/package/@holmrenser/nj)

Takes a mutiple sequence alignment and returns a Newick string. Alphabet (DNA/protein) is auto-detected. Supports optional bootstrap support values on internal nodes. Optionally, only the distance matrix or average distance can be computed. Wrappers use a plugin system for implementing progress/error logging.

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
nj = "0.0.18"
```

```rust
use nj::{NJConfig, NJEvent, SequenceObject, nj, parse_fasta};
use nj::models::SubstitutionModel;

// Parse from a FASTA string
let sequences = parse_fasta(">A\nACGTACGT\n>B\nACCTACGT\n>C\nTCGTACGT\n")?;

// Run Neighbor-Joining
let newick = nj(
    NJConfig {
        msa: sequences,
        n_bootstrap_samples: 100,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
    },
    Some(Box::new(|event| {
        if let NJEvent::BootstrapProgress { completed, total } = event {
            eprintln!("{completed}/{total}");
        }
    })),
)?;
```

Distance-only computation (no tree, no bootstrap):

```rust
use nj::{DistConfig, distance_matrix};
use nj::models::SubstitutionModel;

let result = distance_matrix(DistConfig {
    msa: sequences,
    substitution_model: SubstitutionModel::JukesCantor,
    alphabet: None,
    num_threads: None,
})?;
// result.names — Vec<String> of taxon names
// result.matrix — n×n Vec<Vec<f64>>, symmetric, diagonal zero
```

`average_distance` takes the same `DistConfig` and returns the mean of all pairwise distances as an `f64`.

## Python

```bash
pip install nj_py
```

```python
from nj_py import nj, distance_matrix

msa = [
    {"identifier": "A", "sequence": "ACGTACGT"},
    {"identifier": "B", "sequence": "ACCTACGT"},
    {"identifier": "C", "sequence": "TCGTACGT"},
]

def on_event(event):
    if event["type"] == "BootstrapProgress":
        print(f"{event['completed']}/{event['total']}")

newick = nj(msa, substitution_model="JukesCantor", n_bootstrap_samples=100, on_event=on_event)
```

Distance-only computation:

```python
result = distance_matrix(msa, substitution_model="JukesCantor")
# result["names"] — list of taxon names
# result["matrix"] — n×n list of lists, symmetric, diagonal zero
```

## JavaScript / WASM

```bash
npm install @holmrenser/nj
```

```js
import { nj, distance_matrix } from '@holmrenser/nj';

const msa = [
    { identifier: 'A', sequence: 'ACGTACGT' },
    { identifier: 'B', sequence: 'ACCTACGT' },
    { identifier: 'C', sequence: 'TCGTACGT' },
];

const newick = nj(
    { msa, n_bootstrap_samples: 100, substitution_model: 'JukesCantor' },
    (event) => {
        if (event.type === 'BootstrapProgress') {
            progressBar.value = event.completed / event.total * 100;
        }
    }
);

const { names, matrix } = distance_matrix({ msa, substitution_model: 'JukesCantor' });
```
