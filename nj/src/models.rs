//! Substitution models for computing pairwise evolutionary distances.
//!
//! Each model implements [`ModelCalculation::distance`], mapping a pair of
//! aligned, pre-encoded sequences of equal length to an evolutionary distance.
//! Columns where either sequence has a gap or an unknown character (`N` for DNA, `X`
//! for protein) are excluded (pairwise deletion); the number of comparable columns is
//! the denominator in every formula. An unknown character carries no information about
//! the distance, so it is neither a match nor a mismatch.
//!
//! The count-based models (`PDiff`, `JukesCantor`, `Poisson`, `KimuraProtein`)
//! share one byte kernel ([`diff_and_comparable`], in [`crate::kernels`]) that
//! dispatches to the widest SIMD level at run time. `Kimura2P`, `TajimaNei`, and `Tamura` need
//! per-mismatch classification (transition vs. transversion) or base-composition
//! tallies; they are functions of the 4×4 pair-count table [`DnaCounts`], exposed as
//! `from_counts` on each model.

use serde::{Deserialize, Serialize};

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol, Protein, ProteinSymbol};
use crate::counts::DnaCounts;
use crate::prepared::Prepared;

/// Among-site rate-variation parameters applied on top of a substitution model.
///
/// Both corrections act on the model's final `−c·ln(arg)` step (a no-op for the
/// uncorrected [`PDiff`]):
///
/// * **Gamma rate heterogeneity** (`gamma_shape` = shape `α`): each `−c·ln(arg)`
///   term becomes `c·α·(arg^(−1/α) − 1)` (Jin & Nei 1990). Converges back to
///   `−c·ln(arg)` as `α → ∞`. `None` disables the correction.
/// * **Invariant sites** (`p_invar` = proportion of invariant sites): observed
///   divergence proportions are divided by `(1 − p_invar)` before forming `arg`,
///   and the resulting distance is multiplied by `(1 − p_invar)`.
///
/// The two compose: `d = (1 − p_invar) · Σ c·α·(arg(p/(1−p_invar))^(−1/α) − 1)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RateHet {
    /// Gamma distribution shape parameter `α` (> 0), or `None` for uniform rates.
    pub gamma_shape: Option<f64>,
    /// Proportion of invariant sites in `[0, 1)`; `0.0` disables the correction.
    pub p_invar: f64,
}

impl RateHet {
    /// No rate-variation correction (uniform rates, no invariant sites). Yields
    /// exactly the classic model formulas.
    pub const NONE: RateHet = RateHet {
        gamma_shape: None,
        p_invar: 0.0,
    };
}

impl Default for RateHet {
    fn default() -> Self {
        RateHet::NONE
    }
}

/// Applies the rate-variation correction to a single `−coef·ln(arg)` distance term.
///
/// With no gamma shape this is the classic `−coef·ln(arg)`; with shape `α` it is
/// the gamma-corrected `coef·α·(arg^(−1/α) − 1)`. Callers must guarantee
/// `arg > 0` (the per-model saturation guards do this), so the result is finite
/// and non-negative.
#[inline]
fn corrected_term(coef: f64, arg: f64, gamma_shape: Option<f64>) -> f64 {
    match gamma_shape {
        Some(alpha) => coef * alpha * (arg.powf(-1.0 / alpha) - 1.0),
        None => -coef * arg.ln(),
    }
}

/// Trait for substitution-model distance calculations over a given alphabet.
///
/// A model maps two aligned, pre-encoded sequences of equal length to an
/// evolutionary distance, optionally applying the among-site rate variation in
/// `rates` ([`RateHet`]). Columns where either sequence has a gap are excluded
/// from the comparable-site count (pairwise deletion). Implementations return
/// [`f64::INFINITY`] when the model's formula is undefined for the observed
/// divergence (saturation) — the NJ algorithm handles infinite distances
/// gracefully — and `0.0` when there are no comparable (non-gap) columns.
pub trait ModelCalculation<A: AlphabetEncoding> {
    /// Computes the pairwise distance between two aligned sequences.
    fn distance(s1: &[A::Symbol], s2: &[A::Symbol], rates: RateHet) -> f64;

    /// Distance between sequences `i` and `j` of a prepared alignment.
    ///
    /// The default works on the symbol slices. DNA models override it to read the
    /// pair counts from the bit-planes, which gives the same result (the counts are
    /// exact integers) without a pass over the symbols.
    fn distance_prepared(data: &Prepared<'_, A>, i: usize, j: usize, rates: RateHet) -> f64 {
        Self::distance(data.sequence(i), data.sequence(j), rates)
    }
}

/// Computes the pairwise distance between two aligned sequences using model `M`,
/// applying the rate-variation correction in `rates`.
///
/// `s1` and `s2` must have equal length (i.e. already be aligned). Thin wrapper
/// over [`ModelCalculation::distance`]; the stable call site used by
/// [`crate::distance_matrix`].
#[inline(always)]
pub fn pairwise_distance_with<M, A>(s1: &[A::Symbol], s2: &[A::Symbol], rates: RateHet) -> f64
where
    M: ModelCalculation<A>,
    A: AlphabetEncoding,
{
    M::distance(s1, s2, rates)
}

/// Computes the pairwise distance with no rate-variation correction
/// ([`RateHet::NONE`]). Convenience wrapper over [`pairwise_distance_with`].
#[inline(always)]
pub fn pairwise_distance<M, A>(s1: &[A::Symbol], s2: &[A::Symbol]) -> f64
where
    M: ModelCalculation<A>,
    A: AlphabetEncoding,
{
    pairwise_distance_with::<M, A>(s1, s2, RateHet::NONE)
}

/// Counts `(n_diff, n_comparable)` over two equal-length encoded sequences.
///
/// A column is comparable when both symbols are real states: neither a gap nor the
/// alphabet's unknown symbol (pairwise deletion). `n_diff` is the number of comparable
/// columns where the two symbols differ. Shared by every count-based model — they
/// differ only in how [`ModelCalculation::distance`] maps these two counts to a
/// distance.
#[inline]
fn diff_and_comparable<A: AlphabetEncoding>(s1: &[A::Symbol], s2: &[A::Symbol]) -> (usize, usize) {
    crate::kernels::count_diff_comparable(A::as_bytes(s1), A::as_bytes(s2), A::N_STATES as u8)
}

/// p-distance (proportion of differing sites).
///
/// `d = n_diff / aln_len`
///
/// The simplest possible distance: the raw fraction of alignment columns where
/// the two sequences differ, ignoring gaps. Valid for both DNA and protein.
/// Does not correct for multiple substitutions at the same site.
pub struct PDiff;

impl PDiff {
    /// p-distance from the counts of differing and comparable columns.
    fn from_diff(n_diff: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        n_diff as f64 / n_comparable as f64
    }
}

impl ModelCalculation<DNA> for PDiff {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], _rates: RateHet) -> f64 {
        // Raw p-distance has no log-correction step, so rate-variation
        // parameters do not apply and are ignored.
        let (n_diff, n_comparable) = diff_and_comparable::<DNA>(s1, s2);
        Self::from_diff(n_diff, n_comparable)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, _rates: RateHet) -> f64 {
        let c = data.dna_counts(i, j);
        Self::from_diff(c.n_diff(), c.n_comparable())
    }
}

impl ModelCalculation<Protein> for PDiff {
    fn distance(s1: &[ProteinSymbol], s2: &[ProteinSymbol], _rates: RateHet) -> f64 {
        // Raw p-distance has no log-correction step, so rate-variation
        // parameters do not apply and are ignored.
        let (n_diff, n_comparable) = diff_and_comparable::<Protein>(s1, s2);
        Self::from_diff(n_diff, n_comparable)
    }
}

/// Jukes-Cantor (1969) distance for DNA.
///
/// `d = -0.75 · ln(1 - (4/3) · p)`
///
/// Assumes equal base frequencies and equal substitution rates among all four
/// nucleotides. Corrects for multiple hits at the same site. Returns
/// [`f64::INFINITY`] when `p ≥ 0.75` (the formula is undefined at saturation).
pub struct JukesCantor;

impl JukesCantor {
    /// Distance from the counts of differing and comparable columns.
    fn from_diff(n_diff: usize, n_comparable: usize, rates: RateHet) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        let scale = 1.0 - rates.p_invar;
        let p = (n_diff as f64 / n_comparable as f64) / scale;
        let arg = 1.0 - (4.0 / 3.0) * p;
        if arg <= 0.0 {
            f64::INFINITY // distance undefined (saturation)
        } else {
            scale * corrected_term(0.75, arg, rates.gamma_shape)
        }
    }
}

impl ModelCalculation<DNA> for JukesCantor {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], rates: RateHet) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<DNA>(s1, s2);
        Self::from_diff(n_diff, n_comparable, rates)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        let c = data.dna_counts(i, j);
        Self::from_diff(c.n_diff(), c.n_comparable(), rates)
    }
}

/// Kimura two-parameter (1980) distance for DNA.
///
/// `d = -0.5 · ln(1 - 2p - q) - 0.25 · ln(1 - 2q)`
///
/// Distinguishes transitions (A↔G, C↔T) from transversions (all other
/// substitutions), allowing the two classes to have different rates. Returns
/// [`f64::INFINITY`] when either denominator is non-positive, which occurs at
/// high divergence or when transversion frequency alone reaches 0.5.
pub struct Kimura2P;

impl Kimura2P {
    /// Distance from a pair-count table (`P` and `Q` are the transition and transversion
    /// proportions of the comparable columns).
    pub fn from_counts(c: &DnaCounts, rates: RateHet) -> f64 {
        let n_comparable = c.n_comparable();
        if n_comparable == 0 {
            return 0.0;
        }
        let scale = 1.0 - rates.p_invar;
        let p = (c.transitions() as f64 / n_comparable as f64) / scale;
        let q = (c.transversions() as f64 / n_comparable as f64) / scale;
        let denom1 = 1.0 - 2.0 * p - q;
        let denom2 = 1.0 - 2.0 * q;
        if denom1 <= 0.0 || denom2 <= 0.0 {
            f64::INFINITY // distance undefined (saturation)
        } else {
            scale
                * (corrected_term(0.5, denom1, rates.gamma_shape)
                    + corrected_term(0.25, denom2, rates.gamma_shape))
        }
    }
}

impl ModelCalculation<DNA> for Kimura2P {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], rates: RateHet) -> f64 {
        Self::from_counts(&DnaCounts::from_symbols(s1, s2), rates)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        Self::from_counts(&data.dna_counts(i, j), rates)
    }
}

/// Poisson distance for protein sequences.
///
/// `d = -ln(1 - p)`
///
/// Assumes all amino acid substitutions occur at equal rates (Poisson process).
/// Corrects for multiple hits. Returns [`f64::INFINITY`] when `p ≥ 1.0`.
pub struct Poisson;

impl ModelCalculation<Protein> for Poisson {
    fn distance(s1: &[ProteinSymbol], s2: &[ProteinSymbol], rates: RateHet) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<Protein>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        let scale = 1.0 - rates.p_invar;
        let p = (n_diff as f64 / n_comparable as f64) / scale;
        let arg = 1.0 - p;
        if arg <= 0.0 {
            f64::INFINITY // distance undefined (saturation)
        } else {
            scale * corrected_term(1.0, arg, rates.gamma_shape)
        }
    }
}

/// Kimura (1983) protein distance.
///
/// `d = -ln(1 - p - 0.2 · p²)`
///
/// An empirical correction to the [`Poisson`] distance: the extra `0.2 · p²`
/// term is Kimura's fit to PAM-based distances, giving a better approximation of
/// amino acid divergence at moderate `p` while staying a fast closed form
/// (e.g. ClustalW uses it). Returns [`f64::INFINITY`] once the log argument is
/// non-positive (`p ≳ 0.8541`, the model's saturation point).
pub struct KimuraProtein;

impl ModelCalculation<Protein> for KimuraProtein {
    fn distance(s1: &[ProteinSymbol], s2: &[ProteinSymbol], rates: RateHet) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<Protein>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        let scale = 1.0 - rates.p_invar;
        let p = (n_diff as f64 / n_comparable as f64) / scale;
        let arg = 1.0 - p - 0.2 * p * p;
        if arg <= 0.0 {
            f64::INFINITY // distance undefined (saturation)
        } else {
            scale * corrected_term(1.0, arg, rates.gamma_shape)
        }
    }
}

/// Tajima-Nei (1984) distance for DNA.
///
/// `d = -b · ln(1 - p/b)`,  `b = ½ · (1 - Σ gᵢ² + p²/h)`
///
/// Generalises [`JukesCantor`] to unequal base frequencies: `gᵢ` are the base
/// frequencies and `h = Σ_{i<j} xᵢⱼ² / (2 gᵢ gⱼ)` aggregates the per-pair
/// difference frequencies `xᵢⱼ`, all estimated from the two sequences. Reduces
/// exactly to Jukes-Cantor when base frequencies are equal and substitutions are
/// uniform. Returns [`f64::INFINITY`] at saturation (`p ≥ b`).
pub struct TajimaNei;

impl TajimaNei {
    /// Distance from a pair-count table: base frequencies `gᵢ` come from both sequences,
    /// the difference frequencies `xᵢⱼ` from the off-diagonal cells.
    pub fn from_counts(c: &DnaCounts, rates: RateHet) -> f64 {
        const K: usize = 4;
        let n_comparable = c.n_comparable();
        let n_diff = c.n_diff();
        if n_comparable == 0 {
            return 0.0;
        }
        let n = n_comparable as f64;
        let p = n_diff as f64 / n;
        if p == 0.0 {
            return 0.0;
        }
        // Base frequencies g_i (total base observations = 2 · n_comparable).
        let base_count = c.base_counts();
        let mut g = [0.0f64; K];
        let mut sum_g2 = 0.0;
        for k in 0..K {
            g[k] = base_count[k] as f64 / (2.0 * n);
            sum_g2 += g[k] * g[k];
        }
        // h = Σ_{i<j} x_ij² / (2 g_i g_j); pairs touching a zero-frequency base
        // have x_ij = 0 and are skipped to avoid 0/0.
        let mut h = 0.0;
        for i in 0..K {
            for j in (i + 1)..K {
                let x_ij = c.unordered(i, j) as f64 / n;
                if x_ij > 0.0 && g[i] > 0.0 && g[j] > 0.0 {
                    h += (x_ij * x_ij) / (2.0 * g[i] * g[j]);
                }
            }
        }
        if h <= 0.0 {
            return f64::INFINITY;
        }
        // `b = ½(1 - Σgᵢ² + p²/h)` is invariant under the invariant-sites scaling:
        // dividing both `p` and the per-pair `xᵢⱼ` (hence `h`) by `(1 - p_invar)`
        // leaves `p²/h` unchanged, so `b` is computed from the observed `p`.
        let b = 0.5 * (1.0 - sum_g2 + (p * p) / h);
        let scale = 1.0 - rates.p_invar;
        let arg = 1.0 - (p / scale) / b;
        if b <= 0.0 || arg <= 0.0 {
            f64::INFINITY // distance undefined (saturation)
        } else {
            scale * corrected_term(b, arg, rates.gamma_shape)
        }
    }
}

impl ModelCalculation<DNA> for TajimaNei {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], rates: RateHet) -> f64 {
        Self::from_counts(&DnaCounts::from_symbols(s1, s2), rates)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        Self::from_counts(&data.dna_counts(i, j), rates)
    }
}

/// Tamura (1992) three-parameter distance for DNA.
///
/// `d = -h · ln(1 - P/h - Q) - ½ · (1 - h) · ln(1 - 2Q)`,  `h = 2θ(1-θ)`
///
/// Extends [`Kimura2P`] with a correction for GC-content bias: `P` and `Q` are
/// the transition and transversion proportions and `θ` is the GC content,
/// estimated from both sequences. Reduces exactly to Kimura two-parameter when
/// `θ = 0.5`. Transition/transversion classification matches [`Kimura2P`]; `N`
/// contributes to neither GC content nor transitions. Returns [`f64::INFINITY`]
/// when either logarithm argument is non-positive (saturation or degenerate GC).
pub struct Tamura;

impl Tamura {
    /// Distance from a pair-count table: `P`/`Q` are the transition/transversion
    /// proportions and `θ` the G+C content over both sequences.
    pub fn from_counts(c: &DnaCounts, rates: RateHet) -> f64 {
        let n_comparable = c.n_comparable();
        let (ti, tv) = (c.transitions(), c.transversions());
        if n_comparable == 0 || ti + tv == 0 {
            return 0.0; // no comparable sites, or all comparable sites identical
        }
        let n = n_comparable as f64;
        let scale = 1.0 - rates.p_invar;
        let big_p = (ti as f64 / n) / scale;
        let big_q = (tv as f64 / n) / scale;
        // θ (GC content) is a composition statistic, not a divergence proportion,
        // so it is not rescaled by the invariant-sites correction.
        let base_count = c.base_counts();
        let gc = base_count[1] + base_count[2]; // C and G over both sequences
        let theta = gc as f64 / (2.0 * n);
        let h = 2.0 * theta * (1.0 - theta);
        let arg1 = 1.0 - big_p / h - big_q;
        let arg2 = 1.0 - 2.0 * big_q;
        if h <= 0.0 || arg1 <= 0.0 || arg2 <= 0.0 {
            f64::INFINITY // distance undefined (saturation or degenerate GC)
        } else {
            scale
                * (corrected_term(h, arg1, rates.gamma_shape)
                    + corrected_term(0.5 * (1.0 - h), arg2, rates.gamma_shape))
        }
    }
}

impl ModelCalculation<DNA> for Tamura {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol], rates: RateHet) -> f64 {
        Self::from_counts(&DnaCounts::from_symbols(s1, s2), rates)
    }

    fn distance_prepared(data: &Prepared<'_, DNA>, i: usize, j: usize, rates: RateHet) -> f64 {
        Self::from_counts(&data.dna_counts(i, j), rates)
    }
}

/// Available substitution models.
///
/// | Variant | Alphabet | Formula |
/// |---------|----------|---------|
/// | `PDiff` | DNA, Protein | `p` |
/// | `JukesCantor` | DNA only | `-0.75 · ln(1 - 4p/3)` |
/// | `Kimura2P` | DNA only | `-0.5 · ln(1-2p-q) - 0.25 · ln(1-2q)` |
/// | `TajimaNei` | DNA only | `-b · ln(1 - p/b)` |
/// | `Tamura` | DNA only | `-h · ln(1 - P/h - Q) - 0.5(1-h)·ln(1-2Q)` |
/// | `Poisson` | Protein only | `-ln(1 - p)` |
/// | `KimuraProtein` | Protein only | `-ln(1 - p - 0.2 p²)` |
///
/// Model–alphabet compatibility is enforced at runtime in [`crate::nj`].
#[derive(Clone, Debug, ts_rs::TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum SubstitutionModel {
    /// p-distance: raw proportion of differing sites. No multiple-hit correction.
    PDiff,
    /// Jukes-Cantor (1969): single-rate DNA model with multiple-hit correction.
    JukesCantor,
    /// Kimura two-parameter (1980): separates transition and transversion rates.
    Kimura2P,
    /// Tajima-Nei (1984): DNA model correcting for unequal base frequencies.
    TajimaNei,
    /// Tamura (1992): Kimura two-parameter with a GC-content correction.
    Tamura,
    /// Poisson: equal-rate protein model with multiple-hit correction.
    Poisson,
    /// Kimura (1983): empirical protein distance (`Poisson` plus a `0.2 p²` term).
    KimuraProtein,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{DNA, Protein, dna, protein};

    // --- PDiff DNA ---

    #[test]
    fn test_pdiff_dna_two_differences() {
        assert_eq!(
            pairwise_distance::<PDiff, DNA>(&dna!("ACGTA"), &dna!("AGGTC")),
            0.4
        );
    }

    #[test]
    fn test_pdiff_dna_identical() {
        let s = dna!("ACGT");
        assert_eq!(pairwise_distance::<PDiff, DNA>(&s, &s), 0.0);
    }

    #[test]
    fn test_pdiff_dna_gapped_positions_excluded_from_denominator() {
        // pos 0: A vs T — difference; pos 1: Gap vs C — not comparable (excluded);
        // pos 2: G vs G — same. n_comparable=2, diffs=1, so 1/2.
        assert!(
            (pairwise_distance::<PDiff, DNA>(&dna!("A-G"), &dna!("TCG")) - 1.0 / 2.0).abs() < 1e-12
        );
    }

    // --- PDiff Protein ---

    #[test]
    fn test_pdiff_protein_one_difference() {
        assert_eq!(
            pairwise_distance::<PDiff, Protein>(&protein!("ACDE"), &protein!("ACDF")),
            0.25
        );
    }

    #[test]
    fn test_pdiff_protein_identical() {
        let s = protein!("ARN");
        assert_eq!(pairwise_distance::<PDiff, Protein>(&s, &s), 0.0);
    }

    #[test]
    fn test_pdiff_protein_gaps_not_counted_as_differences() {
        assert_eq!(
            pairwise_distance::<PDiff, Protein>(&protein!("A-D"), &protein!("ARD")),
            0.0
        );
    }

    // --- JukesCantor ---

    #[test]
    fn test_jukes_cantor_dna() {
        let p = 0.4_f64;
        let expected = -0.75 * (1.0 - (4.0 / 3.0) * p).ln();
        assert!(
            (pairwise_distance::<JukesCantor, DNA>(&dna!("ACGTA"), &dna!("AGGTC")) - expected)
                .abs()
                < 1e-6
        );
    }

    #[test]
    fn test_jukes_cantor_identical() {
        let s = dna!("ACGT");
        assert_eq!(pairwise_distance::<JukesCantor, DNA>(&s, &s), 0.0);
    }

    #[test]
    fn test_jukes_cantor_saturated_returns_infinity() {
        // p >= 0.75 means the formula is undefined
        assert_eq!(
            pairwise_distance::<JukesCantor, DNA>(&dna!("AAAA"), &dna!("CGTC")),
            f64::INFINITY
        );
    }

    // --- Kimura2P ---

    #[test]
    fn test_kimura2p_identical() {
        let s = dna!("ACGT");
        assert_eq!(pairwise_distance::<Kimura2P, DNA>(&s, &s), 0.0);
    }

    #[test]
    fn test_kimura2p_pure_transitions() {
        // A↔G and C↔T are transitions; p = 1.0, q = 0.0 → denom1 = -1 ≤ 0 → infinity
        assert_eq!(
            pairwise_distance::<Kimura2P, DNA>(&dna!("ACAC"), &dna!("GTGT")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_kimura2p_pure_transversions() {
        // A↔C, A↔T, G↔C, G↔T are transversions; p = 0, q = 1.0 → denom2 = -1 ≤ 0 → infinity
        assert_eq!(
            pairwise_distance::<Kimura2P, DNA>(&dna!("AAGG"), &dna!("CTCT")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_kimura2p_mixed() {
        // 1 transition (A→G) and 1 transversion (A→C) out of 4 positions
        let p = 0.25_f64;
        let q = 0.25_f64;
        let expected = -0.5 * (1.0 - 2.0 * p - q).ln() - 0.25 * (1.0 - 2.0 * q).ln();
        assert!(
            (pairwise_distance::<Kimura2P, DNA>(&dna!("AATT"), &dna!("GCTT")) - expected).abs()
                < 1e-12
        );
    }

    #[test]
    fn test_kimura2p_saturated_transversions_returns_infinity() {
        // q = 1.0 → denom2 = 1 - 2*1 = -1 ≤ 0 → infinity
        assert_eq!(
            pairwise_distance::<Kimura2P, DNA>(&dna!("AG"), &dna!("CT")),
            f64::INFINITY
        );
    }

    // --- Poisson (Protein) ---

    #[test]
    fn test_poisson_identical() {
        let s = protein!("ARND");
        assert_eq!(pairwise_distance::<Poisson, Protein>(&s, &s), 0.0);
    }

    #[test]
    fn test_poisson_one_difference() {
        // p = 0.25 → d = -ln(0.75)
        let expected = -(1.0_f64 - 0.25).ln();
        assert!(
            (pairwise_distance::<Poisson, Protein>(&protein!("ARND"), &protein!("ARNE"))
                - expected)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn test_poisson_fully_different_returns_infinity() {
        // p = 1.0 → infinity
        assert_eq!(
            pairwise_distance::<Poisson, Protein>(&protein!("AR"), &protein!("DE")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_poisson_gaps_not_counted_as_differences() {
        // only 0 real differences out of 3 positions → p = 0 → d = 0
        assert_eq!(
            pairwise_distance::<Poisson, Protein>(&protein!("A-N"), &protein!("ARN")),
            0.0
        );
    }

    // --- KimuraProtein ---

    #[test]
    fn test_kimura_protein_identical() {
        let s = protein!("ARND");
        assert_eq!(pairwise_distance::<KimuraProtein, Protein>(&s, &s), 0.0);
    }

    #[test]
    fn test_kimura_protein_one_difference() {
        // p = 0.25 → d = -ln(1 - 0.25 - 0.2·0.25²)
        let p = 0.25_f64;
        let expected = -(1.0 - p - 0.2 * p * p).ln();
        assert!(
            (pairwise_distance::<KimuraProtein, Protein>(&protein!("ARND"), &protein!("ARNE"))
                - expected)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn test_kimura_protein_exceeds_poisson() {
        // The 0.2·p² term makes Kimura's correction larger than Poisson's for p>0.
        let (a, b) = (protein!("ARNDC"), protein!("ARKEG")); // 3 diffs of 5 → p = 0.6
        let kimura = pairwise_distance::<KimuraProtein, Protein>(&a, &b);
        let poisson = pairwise_distance::<Poisson, Protein>(&a, &b);
        assert!(kimura > poisson, "kimura={kimura}, poisson={poisson}");
    }

    #[test]
    fn test_kimura_protein_saturated_returns_infinity() {
        // p = 1.0 → arg = 1 - 1 - 0.2 = -0.2 ≤ 0 → infinity
        assert_eq!(
            pairwise_distance::<KimuraProtein, Protein>(&protein!("AR"), &protein!("DE")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_kimura_protein_gaps_not_counted_as_differences() {
        assert_eq!(
            pairwise_distance::<KimuraProtein, Protein>(&protein!("A-N"), &protein!("ARN")),
            0.0
        );
    }

    // --- TajimaNei ---

    #[test]
    fn test_tajima_nei_identical() {
        let s = dna!("ACGT");
        assert_eq!(pairwise_distance::<TajimaNei, DNA>(&s, &s), 0.0);
    }

    #[test]
    fn test_tajima_nei_hand_computed() {
        // ACGT vs AGGT: 1 diff (C↔G) of 4. Base counts A=2,C=1,G=3,T=2 over the
        // 8 observed bases → Σgᵢ² = 0.28125; the single C/G pair gives
        // h = 0.25²/(2·0.125·0.375) = 2/3, so p²/h = 0.09375 and
        // b = ½(1 - 0.28125 + 0.09375) = 0.40625.
        let p = 0.25_f64;
        let b = 0.40625_f64;
        let expected = -b * (1.0 - p / b).ln();
        assert!(
            (pairwise_distance::<TajimaNei, DNA>(&dna!("ACGT"), &dna!("AGGT")) - expected).abs()
                < 1e-12
        );
    }

    #[test]
    fn test_tajima_nei_equals_jukes_cantor_when_uniform() {
        // Tajima-Nei reduces exactly to Jukes-Cantor when base frequencies are
        // equal AND every base-pair difference is equally frequent. The 6 leading
        // columns realise each of the 6 distinct unordered base pairs exactly once
        //   {A,C} {A,G} {A,T} {C,G} {C,T} {G,T}
        // which on its own gives every base count 3; the 12 trailing identical
        // columns (3 each of A/C/G/T) keep gᵢ = 0.25 while lowering p to 1/3 so
        // Jukes-Cantor stays below saturation.
        let s1 = dna!("AGTCTGAAACCCGGGTTT");
        let s2 = dna!("CAAGCTAAACCCGGGTTT");
        let tn = pairwise_distance::<TajimaNei, DNA>(&s1, &s2);
        let jc = pairwise_distance::<JukesCantor, DNA>(&s1, &s2);
        assert!(
            (tn - jc).abs() < 1e-9,
            "Tajima-Nei ({tn}) should match Jukes-Cantor ({jc}) under uniform composition"
        );
    }

    #[test]
    fn test_tajima_nei_saturated_returns_infinity() {
        // AA vs CC: p = 1, g_A = g_C = 0.5, h = 2 → b = 0.5 → 1 - p/b = -1 → infinity
        assert_eq!(
            pairwise_distance::<TajimaNei, DNA>(&dna!("AA"), &dna!("CC")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_tajima_nei_gaps_excluded() {
        // pos1 gapped (excluded); remaining sites identical → p = 0 → d = 0
        assert_eq!(
            pairwise_distance::<TajimaNei, DNA>(&dna!("A-GT"), &dna!("ACGT")),
            0.0
        );
    }

    // --- Tamura ---

    #[test]
    fn test_tamura_identical() {
        let s = dna!("ACGT");
        assert_eq!(pairwise_distance::<Tamura, DNA>(&s, &s), 0.0);
    }

    #[test]
    fn test_tamura_equals_kimura2p_at_gc_half() {
        // Pooled GC content is exactly 0.5 (each sequence has 4 of 8 bases G/C),
        // at which θ = 0.5, h = 0.5, and Tamura's formula collapses to Kimura-2P.
        // One transition (A↔G) and one transversion (C↔A); GC stays balanced.
        let s1 = dna!("ACGTACGT");
        let s2 = dna!("GAGTACGT");
        let tamura = pairwise_distance::<Tamura, DNA>(&s1, &s2);
        let k2p = pairwise_distance::<Kimura2P, DNA>(&s1, &s2);
        assert!(
            (tamura - k2p).abs() < 1e-12,
            "tamura={tamura} should equal k2p={k2p} at GC=0.5"
        );
    }

    #[test]
    fn test_tamura_gc_correction_differs_from_kimura2p() {
        // θ ≠ 0.5: AATT vs GCTT has 1 transition, 1 transversion, GC = 0.25 → h = 0.375.
        let big_p = 0.25_f64;
        let big_q = 0.25_f64;
        let h = 0.375_f64;
        let expected =
            -h * (1.0 - big_p / h - big_q).ln() - 0.5 * (1.0 - h) * (1.0 - 2.0 * big_q).ln();
        let tamura = pairwise_distance::<Tamura, DNA>(&dna!("AATT"), &dna!("GCTT"));
        assert!((tamura - expected).abs() < 1e-12, "tamura={tamura}");
        // And it must genuinely differ from Kimura-2P when GC ≠ 0.5.
        let k2p = pairwise_distance::<Kimura2P, DNA>(&dna!("AATT"), &dna!("GCTT"));
        assert!((tamura - k2p).abs() > 1e-6, "tamura={tamura}, k2p={k2p}");
    }

    #[test]
    fn test_tamura_saturated_returns_infinity() {
        // Pure transitions → P = 1, arg1 = 1 - 1/h - 0 < 0 → infinity
        assert_eq!(
            pairwise_distance::<Tamura, DNA>(&dna!("ACAC"), &dna!("GTGT")),
            f64::INFINITY
        );
    }

    #[test]
    fn test_tamura_gaps_excluded() {
        assert_eq!(
            pairwise_distance::<Tamura, DNA>(&dna!("A-GT"), &dna!("ACGT")),
            0.0
        );
    }

    // --- Rate heterogeneity (+Γ) and invariant sites (+I) ---

    fn gamma(alpha: f64) -> RateHet {
        RateHet {
            gamma_shape: Some(alpha),
            p_invar: 0.0,
        }
    }

    fn pinv(p: f64) -> RateHet {
        RateHet {
            gamma_shape: None,
            p_invar: p,
        }
    }

    #[test]
    fn test_rate_het_none_matches_uncorrected() {
        // RateHet::NONE must reproduce the classic model exactly, including a
        // p_invar of 0.0 written out explicitly.
        let (s1, s2) = (dna!("ACGTA"), dna!("AGGTC"));
        let classic = pairwise_distance::<JukesCantor, DNA>(&s1, &s2);
        assert_eq!(JukesCantor::distance(&s1, &s2, RateHet::NONE), classic);
        assert_eq!(JukesCantor::distance(&s1, &s2, pinv(0.0)), classic);
    }

    #[test]
    fn test_jukes_cantor_gamma_reference_value() {
        // d_Γ = (3α/4)·[(1 − 4p/3)^(−1/α) − 1] (Jin & Nei 1990). p = 2/5.
        let p = 0.4_f64;
        let alpha = 0.5_f64;
        let expected = 0.75 * alpha * ((1.0 - (4.0 / 3.0) * p).powf(-1.0 / alpha) - 1.0);
        let got = JukesCantor::distance(&dna!("ACGTA"), &dna!("AGGTC"), gamma(alpha));
        assert!((got - expected).abs() < 1e-12, "got={got}, expected={expected}");
        // Gamma rate variation inflates the distance relative to the uniform-rate model.
        let plain = pairwise_distance::<JukesCantor, DNA>(&dna!("ACGTA"), &dna!("AGGTC"));
        assert!(got > plain, "gamma={got} should exceed plain JC={plain}");
    }

    #[test]
    fn test_gamma_large_alpha_approaches_uncorrected() {
        // As α → ∞ the gamma correction converges back to −c·ln(arg).
        let (s1, s2) = (dna!("ACGTA"), dna!("AGGTC"));
        let plain = pairwise_distance::<JukesCantor, DNA>(&s1, &s2);
        let big_alpha = JukesCantor::distance(&s1, &s2, gamma(1e6));
        assert!(
            (plain - big_alpha).abs() < 1e-3,
            "plain={plain}, big_alpha={big_alpha}"
        );
    }

    #[test]
    fn test_jukes_cantor_invariant_reference_value() {
        // d_I = (1 − p_inv)·[−¾·ln(1 − (4/3)·(p/(1 − p_inv)))]. p = 2/5, p_inv = 0.2.
        let p = 0.4_f64;
        let p_inv = 0.2_f64;
        let scale = 1.0 - p_inv;
        let expected = scale * (-0.75 * (1.0 - (4.0 / 3.0) * (p / scale)).ln());
        let got = JukesCantor::distance(&dna!("ACGTA"), &dna!("AGGTC"), pinv(p_inv));
        assert!((got - expected).abs() < 1e-12, "got={got}, expected={expected}");
    }

    #[test]
    fn test_kimura2p_gamma_matches_jin_nei() {
        // AATT vs GCTT: 1 transition, 1 transversion → P = Q = 0.25.
        // Jin & Nei (1990): d_Γ = (α/2)[(1−2P−Q)^(−1/α) + ½(1−2Q)^(−1/α) − 3/2].
        let (big_p, big_q, alpha) = (0.25_f64, 0.25_f64, 0.7_f64);
        let arg1 = 1.0 - 2.0 * big_p - big_q;
        let arg2 = 1.0 - 2.0 * big_q;
        let expected =
            0.5 * alpha * (arg1.powf(-1.0 / alpha) + 0.5 * arg2.powf(-1.0 / alpha) - 1.5);
        let got = Kimura2P::distance(&dna!("AATT"), &dna!("GCTT"), gamma(alpha));
        assert!((got - expected).abs() < 1e-12, "got={got}, expected={expected}");
    }

    #[test]
    fn test_poisson_invariant_and_gamma_compose() {
        // Poisson +I+Γ: d = (1−p_inv)·α·[(1 − p/(1−p_inv))^(−1/α) − 1]. p = 1/4.
        let p = 0.25_f64;
        let (alpha, p_inv) = (0.6_f64, 0.1_f64);
        let scale = 1.0 - p_inv;
        let expected = scale * alpha * ((1.0 - p / scale).powf(-1.0 / alpha) - 1.0);
        let rates = RateHet {
            gamma_shape: Some(alpha),
            p_invar: p_inv,
        };
        let got = Poisson::distance(&protein!("ARND"), &protein!("ARNE"), rates);
        assert!((got - expected).abs() < 1e-12, "got={got}, expected={expected}");
    }

    #[test]
    fn test_corrections_preserve_saturation() {
        // Saturated p stays infinite regardless of the gamma shape.
        assert_eq!(
            JukesCantor::distance(&dna!("AAAA"), &dna!("CGTC"), gamma(0.5)),
            f64::INFINITY
        );
        // Invariant-sites rescaling can push an otherwise-finite p past saturation
        // (p = 0.4, p_inv = 0.5 → effective p = 0.8 ≥ 0.75).
        assert_eq!(
            JukesCantor::distance(&dna!("ACGTA"), &dna!("AGGTC"), pinv(0.5)),
            f64::INFINITY
        );
    }

    #[test]
    fn test_pdiff_ignores_rate_variation() {
        let (s1, s2) = (dna!("ACGTA"), dna!("AGGTC"));
        let classic = <PDiff as ModelCalculation<DNA>>::distance(&s1, &s2, RateHet::NONE);
        let with_rates = <PDiff as ModelCalculation<DNA>>::distance(
            &s1,
            &s2,
            RateHet {
                gamma_shape: Some(0.5),
                p_invar: 0.3,
            },
        );
        assert_eq!(classic, with_rates);
    }

    // --- unknown characters (N / X) are missing data, like gaps ---

    #[test]
    fn test_unknown_dna_characters_change_no_distance() {
        // Extra columns where either side is N or a gap (including N/N, which is not a
        // match) must not alter any model's distance.
        let (a, b) = (dna!("ACGTACGTACGT"), dna!("ACGAACGTCCGA"));
        let (an, bn) = (dna!("ACGTACGTACGTNN-A"), dna!("ACGAACGTCCGAAN-N"));
        macro_rules! check {
            ($($m:ty),*) => {$(
                let plain = pairwise_distance::<$m, DNA>(&a, &b);
                let padded = pairwise_distance::<$m, DNA>(&an, &bn);
                assert_eq!(plain.to_bits(), padded.to_bits(), stringify!($m));
            )*};
        }
        check!(PDiff, JukesCantor, Kimura2P, TajimaNei, Tamura);
    }

    #[test]
    fn test_unknown_protein_characters_change_no_distance() {
        let (a, b) = (protein!("ARNDCQEGHILK"), protein!("ARNDCQEGHIAK"));
        let (ax, bx) = (protein!("ARNDCQEGHILKXX-A"), protein!("ARNDCQEGHIAKAX-X"));
        macro_rules! check {
            ($($m:ty),*) => {$(
                let plain = pairwise_distance::<$m, Protein>(&a, &b);
                let padded = pairwise_distance::<$m, Protein>(&ax, &bx);
                assert_eq!(plain.to_bits(), padded.to_bits(), stringify!($m));
            )*};
        }
        check!(PDiff, Poisson, KimuraProtein);
    }

    #[test]
    fn test_unknown_is_neither_match_nor_mismatch() {
        // One N against a base used to count as a difference.
        assert_eq!(pairwise_distance::<PDiff, DNA>(&dna!("ACGTACGTAC"), &dna!("ACGTNCGTAC")), 0.0);
        assert_eq!(
            pairwise_distance::<PDiff, Protein>(&protein!("ARNDCQEGHI"), &protein!("ARNDXQEGHI")),
            0.0
        );
    }

    #[test]
    fn test_pair_counts_agree_with_the_byte_kernel() {
        // The two ways of counting differences (byte kernel, 4×4 table) must agree.
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for len in [0usize, 1, 5, 31, 32, 33, 100, 257] {
            let s1: Vec<DnaSymbol> = (0..len).map(|_| DNA::encode(b"ACGTN-"[(next() % 6) as usize])).collect();
            let s2: Vec<DnaSymbol> = (0..len).map(|_| DNA::encode(b"ACGTN-"[(next() % 6) as usize])).collect();
            let c = DnaCounts::from_symbols(&s1, &s2);
            assert_eq!((c.n_diff(), c.n_comparable()), diff_and_comparable::<DNA>(&s1, &s2), "len {len}");
        }
    }
}
