# How nj.rs estimates evolutionary distances

## Summary (README-ready)

nj.rs is a **distance method**: it first estimates an evolutionary distance for every pair of
sequences, then fits a tree to that matrix with neighbor-joining (NJ). The distance between two
sequences is an estimate of the expected number of substitutions per site on the path connecting
them, under a Markov model of sequence change. Observed differences (`PDiff`) understate it
because sites mutate more than once, so each model corrects for these multiple hits.

"Maximum likelihood" can appear in a distance method without making it an ML tree method. The
likelihood is that of *two* sequences, and its only unknown is one branch length. The familiar
closed-form distances (Jukes–Cantor, Kimura two-parameter) are exactly such maximum-likelihood
estimates; models with no closed form (general reversible matrices, empirical protein matrices)
are solved numerically for the same quantity. The resulting tree is an NJ tree, not an ML tree.
With a correct model, ML distances converge to the true path lengths as sequences get longer, and
NJ recovers the true topology once all distances are accurate enough (Atteson 1999) [12].

Planned extensions (GTR, empirical protein matrices, a model registry) are described in
[`docs/design/substitution-models.md`](../design/substitution-models.md).

## 1. What is being estimated

Sequence change at one site is a continuous-time Markov chain with rate matrix `Q`, equilibrium
frequencies `π`, and `P(t) = exp(Qt)`. `Q` is scaled so the mean rate is 1, so `t` is *expected
substitutions per site*. For a time-reversible model, `π_a P_ab(t) = π_b P_ba(t)`: two sequences
descended from a common ancestor are equivalent to one sequence evolving for time `t = d`, so the
position of the ancestor does not matter. A pairwise distance is therefore the branch length of a
two-leaf tree, which is exactly how PHYLIP describes it: "The distance for each pair of species
estimates the total branch length between the two species" [P1].

For aligned sites that are independent and identically distributed, all the information about `d`
in a pair is the K×K **count matrix** `F`, with `F_ab` = number of columns with state `a` in one
sequence and `b` in the other (columns with gaps or unknown states are dropped: pairwise
deletion). The log-likelihood factorises through `F`,

```
ℓ(t) = Σ_ab F_ab · ln( π_a · P_ab(t) )          (F is sufficient for t)
```

## 2. A family of estimators of one quantity

The same `d` can be estimated from different summaries of `F`. nj.rs currently offers:

| Model | Statistic used | Estimator |
|---|---|---|
| `PDiff` | proportion of differing sites `p` | none (observed distance, biased low) |
| `JukesCantor` (DNA), `Poisson` (protein) | `p` | closed-form inverse of `E[p](t)`; this **is** the ML estimate, because `p` is sufficient when all substitutions are equally likely |
| `Kimura2P` | transition and transversion proportions `P`, `Q` | closed form; equals the ML estimate when `t` and the transition/transversion ratio are both estimated per pair |
| `TajimaNei`, `Tamura` | `p`/`P`/`Q` plus per-pair base composition | analytic estimators from expected divergence under the model; not claimed to be ML |
| `KimuraProtein` | `p` | `d = −ln(1 − p − 0.2p²)`, an empirical approximation of PAM distance [P2][3] |

**Method of moments vs ML.** For `JukesCantor`, `d = −¾ ln(1 − 4p/3)` is obtained by solving
`p = E[p](t) = ¾(1 − e^{−4t/3})` for `t`. In general `E[p](t) = 1 − Σ_a π_a P_aa(t)` is strictly
increasing in `t` for any irreducible reversible model, so it can always be inverted (until `p` reaches
`1 − Σπ_a²`, where the distance is undefined and nj.rs returns infinity). That inversion is a
moment estimator. It is consistent, but it throws away information whenever `p` is not
sufficient. The ML estimator uses all of `F`; for `JukesCantor`, `Poisson` and per-pair
`Kimura2P` the two coincide exactly.

**ML under a fixed model.** For a general reversible `Q` (HKY/F84, GTR, WAG/LG/JTT/Dayhoff) there is
no closed form. The ML distance is the root of the one-dimensional score equation `ℓ′(t) = 0`.
With `Q = Π^{-½} A Π^{½}`, `A = R Λ Rᵀ`, one eigendecomposition makes every evaluation cheap:
`P(t) = Π^{-½} R e^{Λt} Rᵀ Π^{½}` [9]. This is the same estimator as the closed forms above,
computed numerically. It is also what existing software does (section 5).

**A closed form when everything is estimated per pair (empirical GTR).** If the rates and frequencies are
not fixed but all estimated from the pair, there is again a closed form [17][18][19]. Let `J` be the
table `F` as joint frequencies, symmetrised as `(J + Jᵀ)/2` (a reversible model has `J = Π P(t)`, which is
symmetric), and `π` its marginal. A GTR model predicts `Π⁻¹J = exp(Qt)`, so `Qt = ln(Π⁻¹J)` and the
expected number of substitutions per site is `d = −tr(Π ln(Π⁻¹J))`. With `S = Π^{-½} J Π^{-½} = V M Vᵀ`
(symmetric, so one small eigenproblem) this is `d = −Σ_k w_k ln M_k`, `w_k = Σ_i π_i V_ik²`. It is the ML
estimate for a model with a free symmetric `J` (ML estimates are invariant under reparametrisation), so it
is the ML distance of a GTR model with *all* rates free, provided `ln(Π⁻¹J)` is a rate matrix, that is,
has non-negative off-diagonal entries. If it has negative ones the constrained ML estimate differs from
this formula (we measured 2–9% on simulated pairs at distances of 1 to 2); if an eigenvalue of `S` is not
positive, `ln` does not exist and the pair has no estimate. Gamma rates replace `−ln M` by
`α(M^{−1/α} − 1)`, and invariant sites replace `M` by `(M − p_inv)/(1 − p_inv)`. The price of estimating
nine parameters from one pair is variance: the formula fails for a rapidly growing share of pairs beyond
about one substitution per site, and sooner for short sequences (see `docs/design/substitution-models.md`
§5.1 for measured rates).

**When the estimate is infinite.** If at least one column differs, `ℓ(0⁺) = −∞` and `ℓ′(0⁺) = +∞`, so
the maximiser is positive (with no differing column, `d̂ = 0`). As `t → ∞`, `P_ab(t) → π_b` and

```
ℓ′(t) ≈ λ₂ e^{λ₂ t} · C(F),     C(F) = Σ_ab F_ab (Π₂)_ab / π_b
```

where `λ₂ < 0` is the slowest non-zero eigenvalue and `Π₂` the projector onto its eigenspace. If
`C(F) > 0`, the likelihood approaches its limit `ℓ(∞) = Σ F_ab ln π_b` from above, so a finite
maximum exists. For Jukes–Cantor on `K` states this reduces to `p < (K−1)/K`, the familiar
saturation point. If `C(F) ≤ 0` the maximum may still be finite, but the likelihood rises again
towards `ℓ(∞)` in the tail, so the first maximum has to be found by scanning from small `t`. Even
when `C(F) > 0` the gain over `ℓ(∞)` can be negligible, so a finite estimate is not necessarily a
useful one. PHYLIP describes the practical consequence: two sequences may be so dissimilar that the
distance "would have to be infinite, as the likelihood rises indefinitely as the estimated divergence
time increases", and bootstrapped data can create such pairs even when the original data have none
[P1][P2].

**Rate variation across sites.** Gamma-distributed rates (shape `α`) and a proportion `p_inv` of
invariant sites modify the eigenvalues: `e^{λt}` becomes `(1 − λt/α)^{−α}`, and `P` is mixed
with the identity for invariant sites. The corrections therefore stay closed-form for any
reversible model [16][22]. The reported distance must say which scale it uses: with `p_inv`, the
variable-site `t` times `(1 − p_inv)` gives substitutions per site over *all* sites, which is the
convention of the existing `RateHet`.

## 3. Why "ML in a distance method" is not a contradiction

1. **Different objects.** In ML phylogenetics, the likelihood is over the *whole alignment and
   the whole tree*, and the tree is the parameter [6]. Here the likelihood is over *one pair* and
   the parameter is a single branch length. Nothing about the tree topology is estimated by it.
2. **The classical distances already are ML.** PHYLIP states that its Jukes–Cantor distance "is a
   maximum likelihood estimate for that model", that its Kimura two-parameter calculation "is a
   maximum likelihood estimate", and that the F84 distance "is in effect a maximum likelihood
   estimate of the divergence time" [P1]. Nobody considers JC + NJ contradictory; replacing a
   closed-form root by a numerical one changes the algebra, not the logic.
3. **Two stages with separate guarantees.** Stage 1 estimates `d_ij` for each pair; stage 2 fits a
   tree to the matrix. If the model is correct, ML estimators are consistent and asymptotically
   efficient for each `d_ij` (under standard regularity conditions [32]). Atteson showed that NJ returns the
   correct topology whenever every distance error is below half the shortest edge length of the true
   tree [12]; since a finite set of pairwise errors tends to zero as sequences lengthen,
   NJ with consistent distances is a consistent tree estimator. NJ itself is a greedy
   balanced-minimum-evolution algorithm [14].
4. **What is lost.** The tree is not an ML tree. Each distance ignores all other sequences, and NJ
   ignores that estimates share sequences and so have correlated errors. Yang and Rannala conclude
   that distance methods suit large sets of similar sequences, while likelihood and Bayesian methods
   usually have more power and robustness, especially for deep phylogenies [10].
5. **Both approaches depend on the model.** ML distances are consistent only if the model is right.
   Compositional heterogeneity between lineages breaks reversible-model distances; LogDet/paralinear
   distances [20] avoid assuming a model at the cost of efficiency. nj.rs does not claim ML tree
   inference; bootstrap support comes from resampling columns and repeating the whole procedure.

**Variance and weighting.** Asymptotically `Var(d̂) ≈ 1 / (L · I(t))` with Fisher information
`I(t) = Σ_ab π_a P_ab (∂ ln P_ab / ∂t)²`. For Jukes–Cantor this is `p(1−p) / (L (1 − 4p/3)²)`
(the delta-method form; closed-form variances for related corrections go back to [31]). It diverges as `p`
approaches saturation, so long distances are noisy. This motivates weighted least squares [7],
BIONJ's variance model for NJ [13], and WEIGHBOR's likelihood-based weighting of distances [15].

## 4. Protein matrices

- Empirical matrices (Dayhoff/PAM, JTT, WAG, LG) are rate matrices `Q` with fixed exchangeabilities
  and frequencies [23][25][26][27]. Dayhoff's original model was published as probability matrices,
  not rates; the derivation of a rate matrix matters, and DCMut is the version PHYLIP uses [24][P2].
- **BLOSUM is not such a matrix.** It is a log-odds *scoring* matrix at a fixed clustering level
  [28]. Every scoring matrix has implicit target frequencies [29], but that is a joint distribution
  at one divergence, not a model of change through time. Distances computed directly from BLOSUM
  scores (as Biopython's `DistanceCalculator` does, rescaled to 0–1 [B1]) are dissimilarity scores,
  not corrections for multiple hits.
- `Poisson` is the equal-rate special case of the same model family: it is exact when all exchange
  rates are equal, and understates distances otherwise, because observed differences accumulate
  more slowly when substitutions concentrate on a few frequent, easily reversed exchanges.
  `KimuraProtein` is a one-parameter empirical patch for Dayhoff-like data.

## 5. How other software does this

| Software | What it does (as documented) |
|---|---|
| **PHYLIP `dnadist`** | JC, K2P and F84 distances are ML estimates of branch length; base frequencies come from the whole data set, the transition/transversion ratio is preset, gamma and rate categories are supported; output feeds `Neighbor`/`Fitch` [P1] |
| **PHYLIP `protdist`** | Dayhoff (DCMut), JTT, PMB: the distance is the value "which makes the likelihood highest", found by bisection; Kimura's formula is offered as a "rough-and-ready" PAM approximation [P2] |
| **MEGA** | Dayhoff and JTT distances "iteratively computed by a maximum-likelihood method", optionally with gamma rates [M1] |
| **phangorn (R) `dist.ml`** | ML pairwise distances for DNA and many amino-acid models (WAG, LG, JTT, ...), with user-supplied `Q` and base frequencies, gamma, then NJ [R1] |
| **FastME 2.0** | Distance-based (balanced minimum evolution) with DNA and protein distance models including LG, WAG, JTT; the paper's method table marks the matrix models as ML-estimated (read via a secondary summary, paper text not retrievable) [F1] |
| **ape `dist.dna` (R)** | Closed-form distances: JC69, K80, F81, K81, F84, T92, TN93, GG95, LogDet/paralinear [A1] |
| **Biopython `DistanceCalculator`** | Identity or scoring-matrix (e.g. BLOSUM62) dissimilarity, not a model-based correction [B1] |

Common design points: in PHYLIP and phangorn, model parameters (frequencies, ratios) come from the
*whole* alignment or the user while only the branch length is fitted per pair, and in all of them the
distances feed a separate tree-building step.

## 6. Limitations to state in the documentation

- Distances are valid only under the chosen model; model misfit biases them in ways NJ cannot detect.
- Pairwise deletion means different pairs use different sites.
- Variance grows quickly with divergence; saturated pairs have no finite estimate.
- Gamma shape and `p_inv` are user-supplied, not estimated.
- The tree is not an ML tree; ML only describes how each distance is computed.

## References

Numbers refer to this page; `[P*]`, `[M1]`, `[R1]`, `[F1]`, `[A1]`, `[B1]` are tool documentation, read 2026-11-25.

[1] Jukes TH, Cantor CR (1969). Evolution of protein molecules. In Munro HN (ed.), *Mammalian Protein Metabolism*, Academic Press, pp. 21–132.
[2] Kimura M (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. *J Mol Evol* 16:111–120. doi:10.1007/BF01731581
[3] Kimura M (1983). *The Neutral Theory of Molecular Evolution*. Cambridge University Press.
[4] Tajima F, Nei M (1984). Estimation of evolutionary distance between nucleotide sequences. *Mol Biol Evol* 1:269–285.
[5] Tamura K (1992). Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases. *Mol Biol Evol* 9:678–687.
[6] Felsenstein J (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *J Mol Evol* 17:368–376. doi:10.1007/BF01734359
[7] Felsenstein J (1984). Distance methods for inferring phylogenies: a justification. *Evolution* 38:16–24. doi:10.2307/2408542
[8] Felsenstein J (2004). *Inferring Phylogenies*. Sinauer Associates.
[9] Yang Z (2006). *Computational Molecular Evolution*. Oxford University Press (ch. 1 nucleotide models, ch. 2 amino acid and codon models, ch. 4 maximum likelihood).
[10] Yang Z, Rannala B (2012). Molecular phylogenetics: principles and practice. *Nat Rev Genet* 13:303–314. doi:10.1038/nrg3186
[11] Saitou N, Nei M (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol Biol Evol* 4:406–425.
[12] Atteson K (1999). The performance of neighbor-joining methods of phylogenetic reconstruction. *Algorithmica* 25:251–278. doi:10.1007/PL00008277
[13] Gascuel O (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Mol Biol Evol* 14:685–695. doi:10.1093/oxfordjournals.molbev.a025808
[14] Gascuel O, Steel M (2006). Neighbor-joining revealed. *Mol Biol Evol* 23:1997–2000.
[15] Bruno WJ, Socci ND, Halpern AL (2000). Weighted neighbor joining: a likelihood-based approach to distance-based phylogeny reconstruction. *Mol Biol Evol* 17:189–197. doi:10.1093/oxfordjournals.molbev.a026231
[16] Waddell PJ, Steel MA (1997). General time-reversible distances with unequal rates across sites: mixing Γ and inverse Gaussian distributions with invariant sites. *Mol Phylogenet Evol* 8:398–414.
[17] Lanave C, Preparata G, Saccone C, Serio G (1984). A new method for calculating evolutionary substitution rates. *J Mol Evol* 20:86–93. doi:10.1007/BF02101990
[18] Rodríguez F, Oliver JL, Marín A, Medina JR (1990). The general stochastic model of nucleotide substitution. *J Theor Biol* 142:485–501.
[19] Tavaré S (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences* 17:57–.
[20] Lockhart PJ, Steel MA, Hendy MD, Penny D (1994). Recovering evolutionary trees under a more realistic model of sequence evolution. *Mol Biol Evol* 11:605–612.
[21] Gu X, Li W-H (1996). A general additive distance with time-reversibility and rate variation among nucleotide sites. *PNAS* 93:4671–4676. doi:10.1073/pnas.93.10.4671
[22] Jin L, Nei M (1990). Limitations of the evolutionary parsimony method of phylogenetic analysis. *Mol Biol Evol* 7:82–102 (gamma distances).
[23] Dayhoff MO, Schwartz RM, Orcutt BC (1978). A model of evolutionary change in proteins. *Atlas of Protein Sequence and Structure* 5(Suppl 3):345–352.
[24] Kosiol C, Goldman N (2005). Different versions of the Dayhoff rate matrix. *Mol Biol Evol* 22:193–199. doi:10.1093/molbev/msi005
[25] Jones DT, Taylor WR, Thornton JM (1992). The rapid generation of mutation data matrices from protein sequences. *CABIOS* 8:275–282.
[26] Whelan S, Goldman N (2001). A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Mol Biol Evol* 18:691–699.
[27] Le SQ, Gascuel O (2008). An improved general amino acid replacement matrix. *Mol Biol Evol* 25:1307–1320.
[28] Henikoff S, Henikoff JG (1992). Amino acid substitution matrices from protein blocks. *PNAS* 89:10915–10919.
[29] Altschul SF (1991). Amino acid substitution matrices from an information theoretic perspective. *J Mol Biol* 219:555–565.
[30] Müller T, Vingron M (2000). Modeling amino acid replacement. *J Comput Biol* 7:761–776. doi:10.1089/10665270050514918
[31] Kimura M, Ohta T (1972). On the stochastic model for estimation of mutational distance between homologous proteins. *J Mol Evol* 2:87–90. doi:10.1007/BF01653945
[32] Lehmann EL, Casella G (1998). *Theory of Point Estimation*, 2nd ed. Springer. doi:10.1007/b98854

[P1] PHYLIP 3.67 manual, `dnadist`: <https://gensoft.pasteur.fr/docs/phylip/3.67/doc/dnadist.html>
[P2] PHYLIP 3.67 manual, `protdist`: <https://gensoft.pasteur.fr/docs/phylip/3.67/doc/protdist.html>
[M1] MEGA help, Dayhoff and JTT models: <https://www.megasoftware.net/web_help_7/hc_dayhoff_model.htm>
[R1] phangorn reference (`dist.ml`): <https://search.r-project.org/CRAN/refmans/phangorn/html/dist.hamming.html>
[F1] Lefort V, Desper R, Gascuel O (2015). FastME 2.0. *Mol Biol Evol* 32:2798–2800. <https://academic.oup.com/mbe/article/32/10/2798/1212138>
[A1] ape reference (`dist.dna`): <https://search.r-project.org/CRAN/refmans/ape/html/dist.dna.html>
[B1] Biopython `Bio.Phylo.TreeConstruction`: <https://biopython.org/docs/dev/api/Bio.Phylo.TreeConstruction.html>
