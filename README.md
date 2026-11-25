# nj.rs
Neighbor-Joining phylogenetic tree inference in Rust, with Python and WASM bindings.

[![Release](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/release.yml)
[![Tests](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml/badge.svg)](https://github.com/holmrenser/nj.rs/actions/workflows/test.yml)


[![Crates.io Version](https://img.shields.io/crates/v/nj?logo=rust&label=crates.io%3A%20nj)](https://crates.io/crates/nj)
[![PyPI - Version](https://img.shields.io/pypi/v/nj_py?logo=pypi&logoColor=blue&label=pypi%3A%20nj_py)](https://pypi.org/project/nj_py/)
[![NPM Version](https://img.shields.io/npm/v/%40holmrenser%2Fnj?style=flat&logo=npm&label=npm%3A%20%40holmrenser%2Fnj)](https://www.npmjs.com/package/@holmrenser/nj)

Takes a multiple sequence alignment and returns a Newick tree. The alphabet (DNA/protein) is auto-detected. Features:

- Seven distance models, with optional among-site rate variation (gamma shape and/or proportion of invariant sites)
- Optional bootstrap support values on internal nodes (parallel when built with the `parallel` feature)
- Distance-only modes: the full pairwise distance matrix or just the average distance
- A progress/log event callback in every binding

## Substitution models

| Model | DNA | Protein | Rate variation (`gamma_shape`, `p_invar`) |
|-------|:---:|:-------:|:-----------------------------------------:|
| `PDiff` | ✓ | ✓ | ignored |
| `JukesCantor` | ✓ | | ✓ |
| `Kimura2P` | ✓ | | ✓ |
| `TajimaNei` | ✓ | | ✓ |
| `Tamura` | ✓ | | ✓ |
| `Poisson` | | ✓ | ✓ |
| `KimuraProtein` | | ✓ | ✓ |

Using a model that does not match the alphabet returns an error. `gamma_shape` (α) must be `> 0` and `p_invar` must be in `[0, 1)`; they combine as `+I+Γ`.

**Saturated sequence pairs.** When two sequences are too diverged for a model, its distance is infinite, which neighbor-joining cannot use. `nj` sets such pairs to `max_distance` (default `10` substitutions per site; the CLI flag, config field and Python argument are all called `max_distance`), emits a warning event naming how many pairs were affected, and carries on, including in bootstrap replicates. A capped value means "at least this far". The standalone `distance_matrix` and `average_distance` return the raw estimates, including infinity.

## CLI

```bash
cargo install nj --features cli

nj sequences.fasta
nj --substitution-model kimura2-p --n-bootstrap-samples 100 sequences.fasta > tree.nwk
cat sequences.fasta | nj -m jukes-cantor -g 0.5 -i 0.2    # reads stdin when no file is given

# Bootstrap and distance computation run in parallel by default (the `cli`
# feature enables threading). Cap the worker count with -t/--num-threads:
nj -t 4 --n-bootstrap-samples 1000 sequences.fasta > tree.nwk
```

| Option | Description |
|--------|-------------|
| `[FASTA]` | Input alignment (stdin if omitted) |
| `-o, --output <FILE>` | Write output to a file instead of stdout |
| `-b, --n-bootstrap-samples <N>` | Bootstrap replicates, `0` to disable (default `100`) |
| `-m, --substitution-model <MODEL>` | `p-diff` (default), `jukes-cantor`, `kimura2-p`, `tajima-nei`, `tamura`, `poisson`, `kimura-protein` |
| `-a, --alphabet <ALPHABET>` | Override alphabet auto-detection |
| `-g, --gamma-shape <ALPHA>` | Gamma rate-heterogeneity shape |
| `-i, --p-invar <FRACTION>` | Proportion of invariant sites |
| `--max-distance <D>` | Distance given to saturated pairs, which have no finite estimate (default `10`) |
| `--distance-matrix`, `--average-distance` | Also compute these; output becomes a JSON object with a `newick` field |
| `-t, --num-threads <N>` | Worker thread cap |
| `-v, --verbose` | Print stage log messages to stderr |

A progress bar is shown on stderr when bootstrapping.

## Rust

```toml
[dependencies]
nj = "0.0.25"
```

```rust
use nj::{NJConfig, NJEvent, nj, parse_fasta};
use nj::models::SubstitutionModel;

// Parse from a FASTA string
let sequences = parse_fasta(">A\nACGTACGT\n>B\nACCTACGT\n>C\nTCGTACGT\n")?;

// Run Neighbor-Joining
let result = nj(
    NJConfig {
        msa: sequences,
        n_bootstrap_samples: 100,
        substitution_model: SubstitutionModel::JukesCantor,
        alphabet: None,
        num_threads: None,
        return_distance_matrix: false,
        return_average_distance: false,
        max_distance: None,
        gamma_shape: None,
        p_invar: None,
    },
    Some(Box::new(|event| {
        if let NJEvent::BootstrapProgress { completed, total } = event {
            eprintln!("{completed}/{total}");
        }
    })),
)?;
println!("{}", result.newick);
// result.distance_matrix / result.average_distance are `Some` when requested
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
    gamma_shape: None,
    p_invar: None,
})?;
// result.names — Vec<String> of taxon names
// result.matrix — n×n Vec<Vec<f64>>, symmetric, diagonal zero
```

`average_distance` takes the same `DistConfig` and returns the mean of all pairwise distances as an `f64`.

Cargo features: `parallel` (Rayon-based bootstrap and distance computation), `cli` (the binary; implies `parallel`), and `simd` (explicit SIMD kernel, nightly only).

## Python

```bash
pip install nj_py
```

```python
from nj_py import nj, distance_matrix, average_distance

msa = [
    {"identifier": "A", "sequence": "ACGTACGT"},
    {"identifier": "B", "sequence": "ACCTACGT"},
    {"identifier": "C", "sequence": "TCGTACGT"},
]

def on_event(event):
    if event["type"] == "BootstrapProgress":
        print(f"{event['completed']}/{event['total']}")

result = nj(
    msa,
    substitution_model="JukesCantor",
    n_bootstrap_samples=100,
    gamma_shape=0.5,   # optional
    p_invar=0.2,       # optional
    on_event=on_event,
)
print(result["newick"])
```

`nj` returns a dict with a `newick` key. Pass `return_distance_matrix=True` and/or `return_average_distance=True` to also get `distance_matrix` / `average_distance` keys.

Distance-only computation:

```python
result = distance_matrix(msa, substitution_model="JukesCantor")
# result["names"] — list of taxon names
# result["matrix"] — n×n list of lists, symmetric, diagonal zero

mean = average_distance(msa, substitution_model="JukesCantor")
```

## JavaScript / WASM

```bash
npm install @holmrenser/nj
```

```js
import { nj, distance_matrix, average_distance } from '@holmrenser/nj';

const msa = [
    { identifier: 'A', sequence: 'ACGTACGT' },
    { identifier: 'B', sequence: 'ACCTACGT' },
    { identifier: 'C', sequence: 'TCGTACGT' },
];

const { newick } = nj(
    { msa, n_bootstrap_samples: 100, substitution_model: 'JukesCantor' },
    (event) => {
        if (event.type === 'BootstrapProgress') {
            progressBar.value = event.completed / event.total * 100;
        }
    }
);

const { names, matrix } = distance_matrix({ msa, substitution_model: 'JukesCantor' });
const mean = average_distance({ msa, substitution_model: 'JukesCantor' });
```

`gamma_shape`, `p_invar`, `max_distance`, `return_distance_matrix` and `return_average_distance` are accepted in the config object as in Python.
