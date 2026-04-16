//! Substitution models for computing pairwise evolutionary distances.
//!
//! All models implement [`ModelCalculation`] via a three-phase accumulator
//! pattern over aligned column pairs: [`init`](ModelCalculation::init) → per-column
//! [`accumulate`](ModelCalculation::accumulate) → [`finalize`](ModelCalculation::finalize).
//! Gap positions are skipped during accumulation, but the full alignment length
//! is always used as the denominator in the final distance formula.

use serde::{Deserialize, Serialize};

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol, Protein, ProteinSymbol};

/// Trait for substitution model calculations over a given alphabet.
///
/// Implementations follow a three-phase accumulator pattern so the distance
/// computation is a single linear pass over aligned columns:
///
/// 1. [`init`](ModelCalculation::init) — create a zeroed accumulator
/// 2. [`accumulate`](ModelCalculation::accumulate) — update the accumulator for one aligned column
/// 3. [`finalize`](ModelCalculation::finalize) — convert the accumulator to a distance
///
/// Gap positions (`DnaSymbol::Gap` / `ProteinSymbol::Gap`) are ignored in all
/// implementations: a column where either sequence has a gap is not counted as
/// a difference, but the full alignment length is still used as the denominator.
pub trait ModelCalculation<A: AlphabetEncoding> {
    /// Accumulator type — e.g. `usize` for simple difference counts, or
    /// `(usize, usize)` for models that track multiple substitution classes.
    type Acc;

    /// Returns a zeroed accumulator.
    fn init() -> Self::Acc;

    /// Updates `acc` for one aligned column `(a, b)` and returns the updated value.
    fn accumulate(acc: &mut Self::Acc, a: A::Symbol, b: A::Symbol) -> Self::Acc;

    /// Converts the final accumulator to an evolutionary distance.
    ///
    /// `aln_len` is the total number of alignment columns; `n_comparable` is the
    /// number of columns where neither sequence has a gap (pairwise deletion).
    /// Models should use `n_comparable` as the denominator so that gapped
    /// positions are excluded from the distance calculation.
    ///
    /// Returns [`f64::INFINITY`] when the model's formula is undefined for the
    /// observed substitution frequencies (e.g. saturation). The NJ algorithm
    /// handles infinite distances gracefully. Returns `0.0` when `n_comparable`
    /// is zero (no overlapping non-gap columns).
    fn finalize(acc: &Self::Acc, aln_len: usize, n_comparable: usize) -> f64;
}

/// Computes the pairwise distance between two aligned sequences using model `M`.
///
/// `s1` and `s2` must have equal length (i.e. already be aligned). Columns
/// where either sequence carries a gap are excluded from the comparable-site
/// count (`n_comparable`) passed to [`ModelCalculation::finalize`], implementing
/// pairwise deletion.
#[inline(always)]
pub fn pairwise_distance<M, A>(s1: &[A::Symbol], s2: &[A::Symbol]) -> f64
where
    M: ModelCalculation<A>,
    A: AlphabetEncoding,
{
    let aln_len = s1.len();
    let mut acc = M::init();
    let mut n_comparable: usize = 0;

    for k in 0..aln_len {
        if !A::is_gap(s1[k]) && !A::is_gap(s2[k]) {
            n_comparable += 1;
        }
        acc = M::accumulate(&mut acc, s1[k], s2[k]);
    }

    M::finalize(&acc, aln_len, n_comparable)
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
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, _aln_len: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        *acc as f64 / n_comparable as f64
    }
}

impl ModelCalculation<Protein> for PDiff {
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: ProteinSymbol, b: ProteinSymbol) -> Self::Acc {
        if a != b && a != ProteinSymbol::Gap && b != ProteinSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, _aln_len: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        *acc as f64 / n_comparable as f64
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
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, _aln_len: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        let p = *acc as f64 / n_comparable as f64;
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
    type Acc = (usize, usize); // (transitions, transversions)

    fn init() -> Self::Acc {
        (0, 0)
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            match (a, b) {
                (DnaSymbol::A, DnaSymbol::G)
                | (DnaSymbol::G, DnaSymbol::A)
                | (DnaSymbol::C, DnaSymbol::T)
                | (DnaSymbol::T, DnaSymbol::C) => acc.0 += 1, // transition
                _ => acc.1 += 1, // transversion
            }
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, _aln_len: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        let (ti, tv) = *acc;
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
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: ProteinSymbol, b: ProteinSymbol) -> Self::Acc {
        if a != b && a != ProteinSymbol::Gap && b != ProteinSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, _aln_len: usize, n_comparable: usize) -> f64 {
        if n_comparable == 0 {
            return 0.0;
        }
        let p = *acc as f64 / n_comparable as f64;
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
