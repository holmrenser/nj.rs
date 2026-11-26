# Design: generalised substitution models and ML distances

Status: agreed 2026-11-25; the protein ML engine and the `Custom` model are implemented on branch `feat/protein-ml`
(§6). Built-in matrices wait on §5.6; the DNA presets (K2P/HKY/GTR as ML) are not started. Why ML distances belong in a distance method, how other tools do it, and the statistical
guarantees: [`theory/distance-estimation.md`](../theory/distance-estimation.md). Performance work (NJ restructure,
packed pair counts, SIMD) and the `Prepared` alignment it requires: [`performance.md`](performance.md).

## 1. Goal

Add a nucleotide **GTR** distance and **ML distances for protein matrices** (Dayhoff/DCMut, JTT, WAG, LG),
and make further models cheap to add. Today a model touches ~8 files (enum, impl, `dispatch_run!`, lib docs,
README, `CLAUDE.md`, Python `Literal` and docstrings, changelog, generated TS), and
`ModelCalculation::distance` is a stateless function, so parametrised models and cached eigensystems do not fit.

**Non-goals:** tree search or ML tree inference, per-site mixtures, codon models, BLOSUM-based distances
(a scoring matrix, not a rate matrix; theory page §4).

## 2. Decisions

| Question | Decision |
|---|---|
| Model parameters | `π` pooled over the whole alignment (or user-supplied); `κ` and GTR rates are **explicit, no defaults**; protein matrices are fixed tables; only `t` is fitted per pair |
| Gaps, `N`, `X`, ambiguity | Dropped per pair in **every** model, existing ones included |
| Saturated pairs | Replaced by `max_distance` (**default 10.0** substitutions/site) with a warning; never abort |
| Solver | Bracketed, safeguarded Newton; `t_max = 10`; `ℓ(t̂) > ℓ(∞)` required |
| Matrix sourcing and licence | **Open** (§5.6); blocks only the protein steps |

Two behaviour changes need changelog entries: `N`/`X` no longer count as differences, and saturated distances
become finite with a warning instead of `NaN`/`inf` in the output.

## 3. The estimator

Reversible model `Q_ij = s_ij π_j`, mean rate 1, `A = Π^½ Q Π^-½ = R Λ Rᵀ` (symmetric eigensolver, once per
model; Jacobi, [Press 2007 §11.1]). For a pair, `F_ab` counts columns with state `a` in one sequence and `b` in the other.

```
P(t) = Π^-½ R e^{Λt} Rᵀ Π^½
ℓ(t) = Σ_ab F_ab ln P_ab(t),   ℓ′ = Σ F_ab P′_ab / P_ab,   ℓ″ = Σ F_ab (P″_ab P_ab − P′_ab²) / P_ab²
d̂ = argmax ℓ(t)
```

Only nonzero `F_ab` are visited. The equal-rates model reproduces the closed form for any `K`: Jukes–Cantor for `K=4`,
and `−(19/20)·ln(1 − 20p/19)` for `K=20`. The library's `Poisson`, `−ln(1 − p)`, is the `K → ∞` limit and differs from
the 20-state form by about 0.5% at `p = 0.18`, so it is not an ML special case. K2P, Tajima–Nei and Tamura keep their
closed forms (they estimate extra parameters per pair).

**Rate variation, closed form** [Waddell & Steel 1997; Gu & Li 1996]: gamma (shape `α`) replaces `e^{λt}` by
`(1 − λt/α)^{-α}`; invariant sites (`p_inv`) mix `P` with the identity. `Q` is normalised on the variable sites, so
the reported distance is `(1 − p_inv)·t̂` (substitutions per site over all sites, as `RateHet` does today).
`α` and `p_inv` stay user-supplied.

**Empirical GTR** (optional, §5.1): with `F` normalised and symmetrised and `π` its marginal,
`d = −tr(Π ln(Π⁻¹F))` [Lanave 1984; Rodríguez 1990; Tavaré 1986], computed from `S = Π^-½ F Π^-½ = V M Vᵀ` as
`d = −Σ_i π_i (V ln M Vᵀ)_ii`; any eigenvalue `≤ 0` means saturated.

## 4. Architecture

| Layer | Responsibility |
|---|---|
| `PairCounts` | `F` for a pair of sequences, taken from the prepared alignment (packed bit-planes for DNA, byte histogram for protein), dropping columns outside the `K` states. `PDiff` needs only `n_diff`/`n_comparable`. |
| `ReversibleModel` | `π`, `s_ij`, cached eigensystem; constructors `jc`, `k2p(κ)`, `hky(κ)`, `gtr(s)`, protein tables; `P(t)`, `P′(t)`, `P″(t)`, gamma/`p_inv` hooks. |
| Estimator | Closed form (today's), per-pair ML, optional empirical GTR; maps `F` and a model to `d`. |

Models become values: the distance-matrix, NJ and bootstrap code is generic over `D: DistanceModel<A>` and takes the
model as `&D`. Implemented as

```rust
pub trait DistanceModel<A: AlphabetEncoding>: Send + Sync {
    fn distance(&self, data: &Prepared<A>, i: usize, j: usize, rates: RateHet) -> f64;   // INFINITY if saturated
}
pub struct Static<M>(…);   // adapts a stateless `ModelCalculation` type, so closed-form models are unchanged
```

`dispatch_run!` builds the value per arm: `Static::<Kimura2P>::new()` for the closed-form models, or a
`CustomProtein` parsed from the `matrix` option. No `Box<dyn …>`: each model monomorphises the pair loop, and the
`PDiff` and DNA timings did not change. The `build(&msa, cfg)` constructor from the first draft is only needed once a
model depends on the alignment (pooled `π`); today's models carry their parameters.

`ModelConfig` carries `κ`, the six GTR rates and optionally `π`, `α` and `p_inv`. A `dyn` call per pair is
negligible next to the per-pair work, and `dispatch_run!` collapses to one function (about 14 instantiations
to 2); benchmark the `PDiff` SIMD path for regressions.

**Cheap additions.** One registry (name, CLI name, alphabets, description) generates the enum, the compatibility
check, `list_models()` in all bindings, the README and `lib.rs` model tables and the Python `Literal`, with CI
drift checks; docstrings and `CLAUDE.md` point to `list_models()` instead of listing models. Target: a new
model is one impl, one registry entry and a changelog line.

## 5. Decision details

### 5.1 Parameters

- `π`: pooled over all unambiguous residues. Under stationarity every sequence has marginal `π`, so the pooled
  frequency is unbiased and far less noisy than a per-pair estimate; PHYLIP and phangorn do the same.
- `κ` and GTR rates are explicit. PHYLIP presets its ratio and phangorn takes `Q` as an argument; any built-in
  default would be arbitrary. The expected workflow is to estimate them with an ML program on a tree.
- Per-pair composition stays available through Tajima–Nei and Tamura; compositional heterogeneity calls for
  LogDet [Lockhart 1994].
- **Empirical GTR** (12 free parameters per pair) is opt-in and DNA-only. It fails whenever an eigenvalue of
  `S` is non-positive: in simulation 38% of pairs (K=4, L=100, `t`=1.2) and 80% (K=20, L=300, `t`=1).
- **Deferred:** estimating exchangeabilities from the alignment. Candidates are Dayhoff/JTT-style counting from
  close pairs (`F_ab ≈ π_a π_b s_ab t` at small `t`; biased at high divergence [Müller & Vingron 2000]) and a
  pairwise composite likelihood [Varin 2011]; no precedent for the latter was found.

### 5.2 Gaps, `N`, `X` and ambiguity

Today gaps are excluded but `N` counts as a state, so one `N` in ten sites gives `PDiff` 0.1 (measured). PHYLIP
drops `N`, `X`, `?` and `-` as unknown. A fully unknown character contributes `Σ_ab π_a P_ab(t) = 1`, i.e. no
information, so dropping it is exact. Partial IUPAC codes (`R`, `Y`, …) do carry information; dropping them is
conservative, and unbiased if ambiguity is unrelated to the state. A set-likelihood `Σ_{a∈A,b∈B} π_a P_ab(t)`
is a later enhancement; it needs the alphabets to keep IUPAC codes, which are currently collapsed to `N`.

### 5.3 Saturated distances

A saturated `JukesCantor` pair currently makes NJ print `NaN` and `inf` branch lengths (measured): infinite row
sums become `NaN` in the NJ update. Instead:

1. No non-finite value reaches NJ (a bug fix, independent of ML).
2. Non-finite distances become `max_distance` (default 10.0, a field of `NJConfig`; the standalone distance
   functions keep returning raw estimates), with a warning event giving the number of affected pairs (names for
   the main matrix). A cap is a censoring surrogate: "at least this far".
3. Never abort. Bootstrap replicates can saturate even when the original data do not (PHYLIP documents this).
4. Later: treat them as missing and use NJ*/BIONJ*, which handle incomplete matrices in O(n³) [Criscuolo & Gascuel 2008].

### 5.4 Solver

- No differing column: `d̂ = 0`.
- `ℓ′(0⁺) = +∞`, but `ℓ′` can change sign again in the tail (theory page §2). Start from the closed-form estimate,
  expand upward (×2) until `ℓ′ < 0`, then run safeguarded Newton inside the bracket. Saturated if no sign change
  occurs by `t_max = 10` or `ℓ(t̂) ≤ ℓ(∞)`. Judging saturation from the sign of `ℓ′(t_max)` alone would be wrong:
  in simulations at `t`=3 with random 20-state models, all 18 pairs whose tail rose towards `ℓ(∞)` still had an
  interior global maximum.
- Do not assume unimodality: in 19,200 random trials (K=4 and 20, L=100) 16 had two interior maxima, always
  with the second at `t ≥ 11`, beyond `t_max`, with a gain over `ℓ(∞)` of at most 0.013 log-likelihood units.
- Tolerance `|Δt| ≤ 10⁻¹² + 10⁻⁹·t`, far below the standard error (`≈ √(t/L)`). Accuracy at very small `t` is limited by
  cancellation in the eigen-expansion of `P_ab(t)` (about 4·10⁻⁹ relative at `t = 10⁻⁶`).
  `t_max = 10`: Jukes–Cantor at `p = 0.7499` is already 6.7.
- Cost control: one `exp(λt)` vector per iteration, shortcut for identical pairs.

### 5.5 Cost

The solver sits in the O(n²) distance stage, never in NJ's O(n³) loop, but runs once per bootstrap replicate.
Scalar prototype on simulated pairs: DNA 0.23–0.37 µs per pair; protein 2.6–12.6 µs at L=300–1000 (68–292
nonzero entries, ≈4 iterations). The implementation measures 4.5–9.5 µs per pair for a random 20-state model
(`ml_probe`, L=300–1000, divergence up to 30–70%), after working on unordered pairs and starting Newton directly. At n=1000 that is roughly 1.3–6 s per distance matrix, against ≈1.4 s for a
whole `PDiff` run (almost all NJ; measured), times (B+1) with bootstrap. Warm starts did not reduce iterations.
Confirm with a criterion benchmark in the repo before the protein step.

### 5.6 Matrix data and licence (open)

nj.rs is MIT. PAML, the usual source of `.dat` files, is GPL-3 (earlier releases academic-only) according to
Debian metadata; the terms of the individual files were not established. Until decided: do not copy `.dat`
files from PAML or IQ-TREE. Prefer the original publications or the authors' own distribution, record
provenance and terms in each data file, and check invariants (frequencies sum to 1, rows of `Q` sum to 0, mean
rate 1). This is for the maintainer to settle; it is not legal advice.

## 6. Plan

0. **Done:** sanitise non-finite distances before NJ (§5.3): `DistMat::cap_non_finite`, `NJConfig::max_distance`,
   warning events for the main matrix and for bootstrap replicates, and an error from `DistMat::neighbor_joining`
   if a non-finite value ever reaches it.
1. Registry, `list_models()`, generated docs and typings, drift checks. **Not started.**
2. `DistanceModel` values and the `Static` adapter; `dispatch_run!` builds them. **Done** (§4); the alignment-dependent
   `build(&msa, cfg)` and the registry are still open.
3. `Prepared` alignment, pair counts and the missing-data policy (§5.2). **Done**, together with the packed DNA counts
   ([`performance.md`](performance.md), [`performance-log.md`](performance-log.md)).
4. `ReversibleModel` and solver (§5.4), gamma and `p_inv`. **Done** (`reversible.rs`, `ml.rs`), exposed for protein through the
   `Custom` model. **Open:** DNA presets (JC, K2P, HKY/F84, GTR with explicit parameters) and a DNA `Custom`.
5. Protein matrices: `Custom` accepts any matrix as PAML text. **Open:** built-in LG, WAG, JTT, Dayhoff (§5.6): the engine
   embeds a matrix with one `protein_matrix!` line once its data file is cleared.
6. Later: empirical GTR, estimating exchangeabilities, set-likelihood for partial codes, NJ*.

## 7. Validation

- Recover a known `t` from `F = Π·exp(Qt)` for random `Q` (DNA and protein, with gamma and `p_inv`); ML equals
  the closed form for JC and the 20-state equal-rates model, also with gamma and `p_inv`; `α → ∞` and `p_inv → 0` reproduce the uncorrected result.
- `d ≥ 0`, symmetry, identical sequences give `0`.
- Saturation: a saturated `JukesCantor` pair yields no `NaN`/`inf` in the Newick; a constructed count matrix with
  two interior maxima exercises the bracket scan and the `ℓ(t̂) > ℓ(∞)` check.
- Missing data: appending a column with a gap, `N` or `X` in either sequence changes no distance.
- Solver against an independent reference: a dense grid maximum of the full (unmerged) likelihood on 62,500 simulated pairs
  (K=4 and 20, high divergence, short sequences), with no disagreement (`ml::tests`, extensive run `--ignored`).
- Cross-check against a reference tool (PHYLIP `dnadist`/`protdist` or PAML [Yang 2007]) on a small alignment.

## References

Full details and tool documentation links are on the theory page; only works cited here are listed.

- Lanave C et al. (1984). *J Mol Evol* 20:86–93. · Rodríguez F et al. (1990). *J Theor Biol* 142:485–501. · Tavaré S (1986). *Lectures on Mathematics in the Life Sciences* 17:57– (end page not verified).
- Waddell PJ, Steel MA (1997). *Mol Phylogenet Evol* 8:398–414. · Gu X, Li W-H (1996). *PNAS* 93:4671–4676.
- Lockhart PJ et al. (1994). *Mol Biol Evol* 11:605–612.
- Müller T, Vingron M (2000). *J Comput Biol* 7:761–776.
- Criscuolo A, Gascuel O (2008). *BMC Bioinformatics* 9:166.
- Varin C, Reid N, Firth D (2011). *Statistica Sinica* 21:5–42.
- Yang Z (2007). PAML 4. *Mol Biol Evol* 24:1586–1591.
- Press WH et al. (2007). *Numerical Recipes*, 3rd ed., §11.1.
