//! Substitution models for computing pairwise evolutionary distances.
//!
//! Each model implements [`ModelCalculation::distance`], mapping a pair of
//! aligned, pre-encoded sequences of equal length to an evolutionary distance.
//! Columns where either sequence has a gap are excluded (pairwise deletion); the
//! number of comparable (non-gap) columns is the denominator in every formula.
//!
//! The count-based models (`PDiff`, `JukesCantor`, `Poisson`) share a single
//! branchless byte kernel ([`diff_and_comparable`]) that the compiler
//! autovectorizes; with the `simd` Cargo feature it dispatches to an explicit
//! [`core::simd`] kernel. `Kimura2P` must classify each mismatch as a transition
//! or transversion and keeps its own scalar pass.

use serde::{Deserialize, Serialize};

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol, Protein, ProteinSymbol};

/// Trait for substitution-model distance calculations over a given alphabet.
///
/// A model maps two aligned, pre-encoded sequences of equal length to an
/// evolutionary distance. Columns where either sequence has a gap are excluded
/// from the comparable-site count (pairwise deletion). Implementations return
/// [`f64::INFINITY`] when the model's formula is undefined for the observed
/// divergence (saturation) — the NJ algorithm handles infinite distances
/// gracefully — and `0.0` when there are no comparable (non-gap) columns.
pub trait ModelCalculation<A: AlphabetEncoding> {
    /// Computes the pairwise distance between two aligned sequences.
    fn distance(s1: &[A::Symbol], s2: &[A::Symbol]) -> f64;
}

/// Computes the pairwise distance between two aligned sequences using model `M`.
///
/// `s1` and `s2` must have equal length (i.e. already be aligned). Thin wrapper
/// over [`ModelCalculation::distance`]; retained as the stable call site used by
/// [`crate::distance_matrix`].
#[inline(always)]
pub fn pairwise_distance<M, A>(s1: &[A::Symbol], s2: &[A::Symbol]) -> f64
where
    M: ModelCalculation<A>,
    A: AlphabetEncoding,
{
    M::distance(s1, s2)
}

/// Counts `(n_diff, n_comparable)` over two equal-length encoded sequences.
///
/// `n_comparable` is the number of columns where neither symbol is a gap
/// (pairwise deletion); `n_diff` is the number of those comparable columns where
/// the two symbols differ. Shared by every count-based model — they differ only
/// in how [`ModelCalculation::distance`] maps these two counts to a distance.
#[inline]
fn diff_and_comparable<A: AlphabetEncoding>(s1: &[A::Symbol], s2: &[A::Symbol]) -> (usize, usize) {
    count_diff_comparable(A::as_bytes(s1), A::as_bytes(s2), A::GAP_BYTE)
}

/// Branchless scalar kernel: counts `(n_diff, n_comparable)` over raw bytes.
///
/// Written without data-dependent branches (boolean masks accumulated as `0`/`1`)
/// so the compiler can autovectorize the two reductions. Also used for the
/// SIMD path's sub-`N` remainder.
#[inline]
fn count_diff_comparable_scalar(a: &[u8], b: &[u8], gap: u8) -> (usize, usize) {
    let mut n_comparable = 0usize;
    let mut n_diff = 0usize;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let comparable = (x != gap) & (y != gap);
        let diff = comparable & (x != y);
        n_comparable += comparable as usize;
        n_diff += diff as usize;
    }
    (n_diff, n_comparable)
}

/// SIMD kernel (`simd` feature; nightly [`core::simd`]): same contract as
/// [`count_diff_comparable_scalar`], processing `N` bytes per step. Counts set
/// lanes per chunk via the mask bitmask, so the per-chunk popcounts never
/// overflow regardless of sequence length.
#[cfg(feature = "simd")]
#[inline]
fn count_diff_comparable(a: &[u8], b: &[u8], gap: u8) -> (usize, usize) {
    use core::simd::{Simd, cmp::SimdPartialEq};

    const N: usize = 32;
    type V = Simd<u8, N>;

    let gapv = V::splat(gap);
    let mut n_comparable = 0usize;
    let mut n_diff = 0usize;

    let mut ca = a.chunks_exact(N);
    let mut cb = b.chunks_exact(N);
    for (ka, kb) in ca.by_ref().zip(cb.by_ref()) {
        let va = V::from_slice(ka);
        let vb = V::from_slice(kb);
        let comparable = va.simd_ne(gapv) & vb.simd_ne(gapv);
        let diff = comparable & va.simd_ne(vb);
        n_comparable += comparable.to_bitmask().count_ones() as usize;
        n_diff += diff.to_bitmask().count_ones() as usize;
    }

    let (rd, rc) = count_diff_comparable_scalar(ca.remainder(), cb.remainder(), gap);
    (n_diff + rd, n_comparable + rc)
}

/// Scalar fallback used when the `simd` feature is disabled (stable toolchains).
#[cfg(not(feature = "simd"))]
#[inline]
fn count_diff_comparable(a: &[u8], b: &[u8], gap: u8) -> (usize, usize) {
    count_diff_comparable_scalar(a, b, gap)
}

/// p-distance (proportion of differing sites).
///
/// `d = n_diff / aln_len`
///
/// The simplest possible distance: the raw fraction of alignment columns where
/// the two sequences differ, ignoring gaps. Valid for both DNA and protein.
/// Does not correct for multiple substitutions at the same site.
pub struct PDiff;

impl ModelCalculation<DNA> for PDiff {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol]) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<DNA>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        n_diff as f64 / n_comparable as f64
    }
}

impl ModelCalculation<Protein> for PDiff {
    fn distance(s1: &[ProteinSymbol], s2: &[ProteinSymbol]) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<Protein>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        n_diff as f64 / n_comparable as f64
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

impl ModelCalculation<DNA> for JukesCantor {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol]) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<DNA>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        let p = n_diff as f64 / n_comparable as f64;
        if p >= 0.75 {
            f64::INFINITY // distance undefined
        } else {
            -0.75 * (1.0 - (4.0 / 3.0) * p).ln()
        }
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

impl ModelCalculation<DNA> for Kimura2P {
    fn distance(s1: &[DnaSymbol], s2: &[DnaSymbol]) -> f64 {
        // Kimura2P classifies each mismatch as a transition (A↔G, C↔T) or a
        // transversion, so it keeps its own scalar pass rather than the shared
        // count kernel. Gap columns are excluded (pairwise deletion).
        let mut n_comparable = 0usize;
        let mut ti = 0usize;
        let mut tv = 0usize;
        for (&a, &b) in s1.iter().zip(s2.iter()) {
            if a != DnaSymbol::Gap && b != DnaSymbol::Gap {
                n_comparable += 1;
                if a != b {
                    match (a, b) {
                        (DnaSymbol::A, DnaSymbol::G)
                        | (DnaSymbol::G, DnaSymbol::A)
                        | (DnaSymbol::C, DnaSymbol::T)
                        | (DnaSymbol::T, DnaSymbol::C) => ti += 1, // transition
                        _ => tv += 1,                              // transversion
                    }
                }
            }
        }
        if n_comparable == 0 {
            return 0.0;
        }
        let p = ti as f64 / n_comparable as f64;
        let q = tv as f64 / n_comparable as f64;
        let denom1 = 1.0 - 2.0 * p - q;
        let denom2 = 1.0 - 2.0 * q;
        if denom1 <= 0.0 || denom2 <= 0.0 {
            f64::INFINITY // distance undefined
        } else {
            -0.5 * denom1.ln() - 0.25 * denom2.ln()
        }
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
    fn distance(s1: &[ProteinSymbol], s2: &[ProteinSymbol]) -> f64 {
        let (n_diff, n_comparable) = diff_and_comparable::<Protein>(s1, s2);
        if n_comparable == 0 {
            return 0.0;
        }
        let p = n_diff as f64 / n_comparable as f64;
        if p >= 1.0 {
            f64::INFINITY // distance undefined
        } else {
            -(1.0 - p).ln()
        }
    }
}

/// Available substitution models.
///
/// | Variant | Alphabet | Formula |
/// |---------|----------|---------|
/// | `PDiff` | DNA, Protein | `p` |
/// | `JukesCantor` | DNA only | `-0.75 · ln(1 - 4p/3)` |
/// | `Kimura2P` | DNA only | `-0.5 · ln(1-2p-q) - 0.25 · ln(1-2q)` |
/// | `Poisson` | Protein only | `-ln(1 - p)` |
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
    /// Poisson: equal-rate protein model with multiple-hit correction.
    Poisson,
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
}
