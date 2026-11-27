# Design: generalised substitution models and ML distances

Status: agreed 2026-11-25; the protein ML engine, the `Custom` model and six built-in matrices (branch `feat/protein-ml`) and
the DNA `Gtr` model (branch `feat/dna-gtr`), `F84` and `EmpiricalGtr` (branch `feat/dna-gtr-extras`) and `Hky` (branch `feat/hky-and-ml-no-alloc`) are implemented (§6). Estimating exchangeabilities is not. Why ML distances belong in a distance method, how other tools do it, and the statistical
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
| Matrix sourcing and licence | Source chosen by the maintainer, caveats recorded (§5.6) |

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
- **Empirical GTR** (nine free parameters per pair: three frequencies, five relative rates and the distance) is
  opt-in and DNA-only (`EmpiricalGtr`). It fails whenever an eigenvalue of `S` is non-positive (the pair is then
  reported saturated and capped by `max_distance`). Share of simulated pairs that fail, an HKY-like model
  (κ = 4, π = 0.3/0.2/0.2/0.3), 400 pairs per cell; the ML distance with *explicit* rates saturated in ≤ 4% of
  pairs, in the 100-site, `t` = 2.5 cell only:

  | sites | t = 0.3 | 0.6 | 1.0 | 1.5 | 2.5 |
  |---|---|---|---|---|---|
  | 100 | 0% | 0.5% | 26% | 65% | 89% |
  | 300 | 0% | 0% | 2% | 34% | 83% |
  | 1000 | 0% | 0% | 0% | 5.5% | 69% |
  | 5000 | 0% | 0% | 0% | 0% | 46% |

  Where it exists it is unbiased (median distance ≈ `t`). Use it for alignments of a few hundred sites or more
  with distances below about one substitution per site; otherwise `Gtr` with rates fitted on a tree. Reproduce
  with `cargo test --release -p nj --lib empirical_gtr_failure -- --ignored --nocapture`. The same formula
  would also run for 20 states, but fails more often there (80% at L=300, `t`=1 in the early simulations), so it is
  not exposed for protein.
- **Embeddability.** The closed form equals the ML estimate of a GTR model with all rates free only when
  `ln(Π⁻¹J)` has non-negative off-diagonal entries. Otherwise PAML's constrained ML (REV) and the closed form differ:
  of 28 simulated pairs (600 sites, distances up to 2) 23 were valid and agreed with `baseml` to 4·10⁻⁶; of the other
  five, three had a negative implied rate and differed by 2–9% (and two had no closed form while `baseml` gave a
  finite distance). The closed form is the one in the literature [Lanave 1984; Rodríguez 1990; Tavaré 1986] and is
  kept; the constrained estimate would need a numerical fit of all rates.
- **HKY85** (`Hky`, explicit `kappa`) is `Gtr` with `[1, κ, 1, 1, κ, 1]`; κ is the transition/transversion *rate* ratio
  (Hasegawa, Kishino & Yano 1985, *J Mol Evol* 22:160–174); `κ = 1` is F81. Checked against `baseml` `model = 4`.
- **F84** (`F84`, explicit `kappa`) is `Gtr` with exchangeabilities `[1, 1 + κ/π_R, 1, 1, 1 + κ/π_Y, 1]`; `κ = 0` is F81.
  It is *not* HKY85's rate ratio (the two transition rates differ unless `π_R = π_Y`). Its κ is the one PAML's `baseml`
  uses for `model = 3`, which was verified numerically rather than taken from a reference.
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
- Stopping: Newton converges quadratically, with error `≈ c·s²` after a step of relative size `s` and `c` of order 1
  (measured over consecutive Newton steps of simulated pairs: median 1.2–1.9, maximum 2.9). A Newton step of at most
  `3·10⁻⁵·t` is applied and ends the iteration, leaving a relative error of at most about 2.6·10⁻⁹, far below the
  standard error (`≈ √(t/L)`). Bisection and doubling steps stop only at `10⁻⁸`. A Newton step of a few ulps counts as
  converged: the slope's rounding noise is as large as such a step, and treating it as unusable sent the solver
  doubling and bisecting for ~30 evaluations. The start is the equal-rates closed form including gamma and `p_inv`.
  Accuracy at very small `t` is limited by cancellation in the eigen-expansion of `P_ab(t)` (about 4·10⁻⁹ relative at
  `t = 10⁻⁶`). `t_max = 10`: Jukes–Cantor at `p = 0.7499` is already 6.7.
- Cost control: one `exp(λt)` vector per iteration (none for the zero eigenvalue); the log-likelihood is summed once at
  the end (the iteration needs only `ℓ′` and `ℓ″`); shortcut for identical pairs.

### 5.5 Cost

The solver sits in the O(n²) distance stage, never in NJ's O(n³) loop, but runs once per bootstrap replicate.
Scalar prototype on simulated pairs: DNA 0.23–0.37 µs per pair; protein 2.6–12.6 µs at L=300–1000 (68–292
nonzero entries, ≈4 iterations). The implementation measures 3.6–6.1 µs per pair for a random 20-state model
(`ml_probe`, L=300–1000, divergence up to 30–70%), after working on unordered pairs, starting Newton directly and the stopping rule of §5.4; DNA
takes 0.12 µs per pair beyond counting (`performance-log.md`). At n=1000 that is roughly 1.3–6 s per distance matrix, against ≈1.4 s for a
whole `PDiff` run (almost all NJ; measured), times (B+1) with bootstrap. Warm starts did not reduce iterations.
Confirm with a criterion benchmark in the repo before the protein step.

### 5.6 Matrix data and licence

nj.rs is MIT. PAML, the usual source of `.dat` files, is GPL-3 (earlier releases academic-only) according to Debian
metadata; the terms of the individual files were not established, so no file was copied from PAML or IQ-TREE. The
authors' own servers (`ebi.ac.uk/goldman-srv`, `atgc.lirmm.fr`) could not be reached from the development sandbox. The
maintainer chose the list at `giphy.pasteur.fr/empirical-models-of-amino-acid-substitution/` as the source (2026-11-26):
it gives PAML-layout files with, for each, the reference and the origin of the data, and states no terms of use. The
six files (LG, WAG, JTT, JTT-DCMut, Dayhoff, DCMut-Dayhoff) are embedded unmodified in `nj/data/matrices/`, with
references, stated data sources, SHA-256 sums and the licence caveat in that directory's `README.md`. Invariants are
tested for every matrix (frequencies sum to 1, `P(t)` is stochastic and reversible, mean rate 1, residue order checked
against chemically conservative and dissimilar pairs). The JTT file's stated source is the GPL-licensed MrBayes
(`src/model.c`); whether that matters for an MIT crate is for the maintainer to settle; this is not legal advice.

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
   `Custom` model and for DNA through `Gtr` (`gtr.rs`: six explicit exchangeabilities, base frequencies pooled from the
   alignment unless given). HKY85, Tamura–Nei and F81 are `Gtr` with particular rates; JC and K2P keep their closed forms,
   which estimate κ per pair. `Hky` and `F84` are named presets of `Gtr` (explicit κ) and `EmpiricalGtr` the closed-form GTR
   distance (§5.1). A DNA `Custom` would be `Gtr` again (a reversible 4×4 matrix with frequencies); only non-reversible
   matrices would add something, and the eigendecomposition does not support them.
5. Protein matrices. **Done:** `Custom` takes any matrix as PAML text, and LG, WAG, JTT, JTT-DCMut, Dayhoff and DCMut-Dayhoff
   are built in (§5.6): one data file and one `protein_matrix!` line each.
6. Later: estimating exchangeabilities, set-likelihood for partial codes, NJ*.

## 7. Validation

- Recover a known `t` from `F = Π·exp(Qt)` for random `Q` (DNA and protein, with gamma and `p_inv`); ML equals
  the closed form for JC and the 20-state equal-rates model, also with gamma and `p_inv`; `α → ∞` and `p_inv → 0` reproduce the uncorrected result.
- `d ≥ 0`, symmetry, identical sequences give `0`.
- Saturation: a saturated `JukesCantor` pair yields no `NaN`/`inf` in the Newick; a constructed count matrix with
  two interior maxima exercises the bracket scan and the `ℓ(t̂) > ℓ(∞)` check.
- Missing data: appending a column with a gap, `N` or `X` in either sequence changes no distance.
- Solver against an independent reference: a dense grid maximum of the full (unmerged) likelihood on 62,500 simulated pairs
  (K=4 and 20, high divergence, short sequences), with no disagreement (`ml::tests`, extensive run `--ignored`).
- **Cross-check against PHYLIP 3.697 `protdist` (done).** `JttDcmut` and `DayhoffDcmut` reproduce `protdist`'s JTT and PAM
  distances to its printed precision (six decimals, ≤ 1·10⁻⁵ relative) on an 8-taxon, 400-residue alignment with up to 53%
  of columns resampled, for uniform rates and gamma shapes 1 and 0.5, over distances up to 9.6 substitutions per site.
  The one pair whose `protdist` distance is 12.6 lies beyond `t_max` and is reported saturated (§5.3). The original `Jtt`
  and `Dayhoff` files differ from `protdist` by 0.2–0.5%, which is the size of the data-version differences Kosiol and
  Goldman (2005) report, so `protdist` evidently uses the recomputed matrices. A five-taxon excerpt is kept as a
  regression test (`nj/tests/protdist_reference.rs`).
- **Cross-check of `Gtr` against PAML 4.9j `baseml` (done).** HKY85 with a fixed κ, one two-sequence run per pair (so
  `baseml`'s frequencies are those of the pair, as `Gtr`'s default pools over the sequences present): all 28 pairs of an
  8-taxon, 600-site alignment (composition skewed, 70:30 transition bias, distances up to 2.2) at κ = 0.7, 3 and 10 agree to
  `baseml`'s printed precision (≤ 5·10⁻⁶ absolute). With gamma shape 0.5 `baseml`'s discrete gamma approaches the exact
  continuous value used here as its category count grows (8.3% apart at 5 categories, 2.7% at 10, 0.4% at 40). `baseml`'s
  `runmode = -2` is not implemented for DNA and its TN93 could not be fixed per pair, so AG and CT are separated by a test
  with fully distinct exchangeabilities against an unrelated matrix exponential instead. A ten-pair excerpt is kept as
  `nj/tests/paml_reference.rs`.
- **F84 and empirical GTR against `baseml` (done).** F84 with a fixed κ (`model = 3`): all 28 pairs at κ = 0.5, 2 and 5
  agree to ≤ 5·10⁻⁶, so κ means the same thing in both. Empirical GTR against REV (`model = 7`, all rates free, observed
  frequencies, which on a two-sequence tree is the saturated symmetric model): the 23 pairs whose implied matrix is a
  valid rate matrix agree to ≤ 4·10⁻⁶ (distances up to 1.98); the other five do not, as explained under *Embeddability*
  (§5.1). Both are in `nj/tests/paml_reference.rs` (REV only for the valid pairs).
- Exact tables: the empirical GTR distance returns `t` exactly from the expected pair table of GTR models built with an
  independent matrix exponential, and agrees with the ML solver under gamma and invariant sites.

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
