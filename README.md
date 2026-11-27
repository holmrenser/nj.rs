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
| `Gtr` (ML distance, explicit rates) | ✓ | | ✓ |
| `Hky` (ML distance, explicit `kappa`) | ✓ | | ✓ |
| `F84` (ML distance, explicit `kappa`) | ✓ | | ✓ |
| `EmpiricalGtr` (closed form, estimated per pair) | ✓ | | ✓ |
| `Poisson` | | ✓ | ✓ |
| `KimuraProtein` | | ✓ | ✓ |
| `Lg`, `Wag`, `Jtt`, `JttDcmut`, `Dayhoff`, `DayhoffDcmut` (ML distance, empirical matrix) | | ✓ | ✓ |
| `Custom` (ML distance, your rate matrix) | | ✓ | ✓ |

**GTR (maximum likelihood, DNA).** `Gtr` gives, for each pair of sequences, the distance `t` that maximises the likelihood of its aligned base pairs under a general time-reversible model, with gamma rates and invariant sites if requested. The six exchangeabilities are **explicit**, in the order `AC, AG, AT, CG, CT, GT` (`gtr_rates`; CLI `--gtr-rates 1,4,1,1,4,1`; Python `gtr_rates=[1, 4, 1, 1, 4, 1]`), because any built-in default would be arbitrary: estimate them with a maximum-likelihood program on a tree (IQ-TREE, RAxML, PAML) and pass them in. Only their ratios matter. The base frequencies `A, C, G, T` (`base_frequencies`, `--base-frequencies`) default to those pooled over all unambiguous bases of the alignment, and, like the exchangeabilities, are fixed model parameters rather than fitted per pair. Special cases: HKY85 is `[1, κ, 1, 1, κ, 1]` (κ the transition/transversion *rate* ratio, as in PAML; see `Hky` below), Tamura–Nei `[1, κ₁, 1, 1, κ₂, 1]`, F81 all rates equal, and Jukes–Cantor all rates equal with frequencies `0.25`. This differs from `Kimura2P`, `TajimaNei` and `Tamura`, which estimate their extra parameters from each pair. Gaps and `N` are dropped per pair, saturated pairs go through `max_distance`, and `Gtr` is DNA-only. HKY85 distances reproduce PAML `baseml` (fixed κ, per pair) to its printed precision. Cost: about 0.15 s for 1000 sequences of 1000 sites.

**HKY85.** `Hky` is the ML distance under the model of Hasegawa, Kishino and Yano (1985) with an explicit `kappa >= 0` (`--kappa`, `kappa=`), the ratio of the transition rate to the transversion rate; `kappa = 1` is F81. It is `Gtr` with rates `[1, κ, 1, 1, κ, 1]`, named so that the common case needs no list of six numbers. Base frequencies are as for `Gtr`. Distances reproduce PAML `baseml` (`model = 4`, fixed `κ`) to its printed precision, and, unlike `Kimura2P`, `κ` is fixed rather than estimated from each pair.

**F84.** `F84` is the ML distance under Felsenstein's (1984) model with an explicit `kappa >= 0` (`--kappa`, `kappa=`): transversions have rate 1 and the two kinds of transition get the extra rates `κ/π_R` (A↔G) and `κ/π_Y` (C↔T), with `π_R = π_A + π_G` and `π_Y = π_C + π_T`, so `κ = 0` is F81. This is the `κ` of PAML's F84, not the transition/transversion rate ratio of HKY85 (F84 has two different transition rates unless `π_R = π_Y`). Base frequencies are as for `Gtr`; PAML `baseml` reproduces the distances to its printed precision.

**Empirical GTR.** `EmpiricalGtr` needs no parameters: it is the closed-form distance of a general time-reversible model whose rates *and* frequencies are all estimated from each pair (Lanave et al. 1984; Tavaré 1986; Rodríguez et al. 1990). It is the maximum-likelihood distance of a GTR model with all rates free whenever the estimated rate matrix is valid (non-negative rates, the usual case at low divergence), and agrees with PAML `baseml`'s REV model there. The price is variance: it fails for a growing share of pairs beyond about one substitution per site, and sooner for short sequences (about 2% of pairs at 300 sites and one substitution per site, 26% at 100 sites, 0% at 1000 sites). A pair without an estimate is infinite and `nj` caps it at `max_distance` with a warning, so use it for alignments of a few hundred sites or more at low divergence, and `Gtr` with rates fitted on a tree otherwise. DNA only.

**Empirical protein matrices (maximum likelihood).** `Lg`, `Wag`, `Jtt`, `JttDcmut`, `Dayhoff` and `DayhoffDcmut` give, for every pair of sequences, the distance `t` that maximises the likelihood of its aligned residue pairs under the named replacement matrix, with gamma rates and invariant sites if requested. They differ from the closed-form `Poisson` and `KimuraProtein` in using the matrix's unequal exchange rates and residue frequencies. `JttDcmut` and `DayhoffDcmut` are the versions recomputed by Kosiol and Goldman, who recommend them over the originals. Please cite the matrix you use:

- **LG**: Le SQ, Gascuel O (2008) An improved general amino acid replacement matrix. *Mol Biol Evol* 25:1307–1320. [doi:10.1093/molbev/msn067](https://doi.org/10.1093/molbev/msn067)
- **WAG**: Whelan S, Goldman N (2001) A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Mol Biol Evol* 18:691–699. [doi:10.1093/oxfordjournals.molbev.a003851](https://doi.org/10.1093/oxfordjournals.molbev.a003851)
- **JTT**: Jones DT, Taylor WR, Thornton JM (1992) The rapid generation of mutation data matrices from protein sequences. *Comput Appl Biosci* 8:275–282. [doi:10.1093/bioinformatics/8.3.275](https://doi.org/10.1093/bioinformatics/8.3.275)
- **Dayhoff**: Dayhoff MO, Schwartz RM, Orcutt BC (1978) A model of evolutionary change in proteins. In: *Atlas of Protein Sequence and Structure*, Vol. 5, pp. 345–352. National Biomedical Research Foundation, Washington DC
- **JTT-DCMut, DCMut-Dayhoff**: Kosiol C, Goldman N (2005) Different versions of the Dayhoff rate matrix. *Mol Biol Evol* 22:193–199. [doi:10.1093/molbev/msi005](https://doi.org/10.1093/molbev/msi005)

`JttDcmut` and `DayhoffDcmut` reproduce the distances of PHYLIP 3.697 `protdist` (JTT and PAM models) to its printed precision, with and without gamma rates. The data files, their sources and checksums are in [`nj/data/matrices/`](nj/data/matrices/README.md); they were taken from the list at <http://giphy.pasteur.fr/empirical-models-of-amino-acid-substitution/>.

**Custom protein matrices.** `Custom` is the same maximum-likelihood distance under a reversible rate matrix you supply (another matrix, or one of the above from its authors): for every pair of sequences it fits the single distance `t` that maximises `Σ F_ab ln P_ab(t)` (`F_ab`: columns with residues `a` and `b`; `P(t) = exp(Qt)`), so the matrix is fixed and only `t` is estimated. Gamma rates and invariant sites enter the likelihood, and gaps and `X` are dropped per pair. Pass the matrix as text in PAML layout in the `matrix` option (CLI: `--matrix FILE`): the lower triangle of the symmetric exchangeability matrix by rows (1, 2, …, 19 numbers), then the 20 amino-acid frequencies, in the order `ARNDCQEGHILKMFPSTWYV`; anything after those numbers is ignored. `matrix` is required for `Custom` and an error for any other model. PAML-format matrix files can be passed as they are. Cost: about 5–10 µs per pair (n = 1000: several seconds), against milliseconds for the closed-form models. Theory and solver: [`docs/design/substitution-models.md`](docs/design/substitution-models.md), [`docs/theory/distance-estimation.md`](docs/theory/distance-estimation.md).

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
| `-m, --substitution-model <MODEL>` | `p-diff` (default), `jukes-cantor`, `kimura2-p`, `tajima-nei`, `tamura`, `gtr`, `hky`, `f84`, `empirical-gtr`, `poisson`, `kimura-protein`, `lg`, `wag`, `jtt`, `jtt-dcmut`, `dayhoff`, `dayhoff-dcmut`, `custom` |
| `-a, --alphabet <ALPHABET>` | Override alphabet auto-detection |
| `--gtr-rates <AC,AG,AT,CG,CT,GT>` | Exchangeabilities for `-m gtr` (required there) |
| `--kappa <K>` | `kappa` for `-m hky` and `-m f84` (required there) |
| `--base-frequencies <A,C,G,T>` | Base frequencies for `-m gtr`, `-m hky` and `-m f84` (default: pooled from the alignment) |
| `--matrix <FILE>` | Rate matrix (PAML layout) for `-m custom` |
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

Cargo features: `parallel` (Rayon-based bootstrap and distance computation) and `cli` (the binary; implies `parallel`). SIMD kernels are always on (stable Rust, instruction set chosen at run time); `simd` is a deprecated no-op.

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
