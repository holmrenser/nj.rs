//! Alphabet definitions and encoding for DNA and protein sequences.
//!
//! Each alphabet is represented by a marker struct ([`DNA`], [`Protein`]) that
//! implements [`AlphabetEncoding`]. The associated `Symbol` types are enums
//! whose variants cover the standard characters plus a gap and an ambiguity
//! symbol (`N` for DNA, `X` for protein). Unknown input bytes are silently
//! mapped to the ambiguity symbol rather than panicking.

use serde::{Deserialize, Serialize};

/// Trait that defines a biological sequence alphabet and its byte encoding.
///
/// Implementations map raw bytes (as they appear in a FASTA file) to typed
/// symbol values. Unknown or ambiguous characters are mapped to the alphabet's
/// "wildcard" variant rather than producing an error.
///
/// The trait is generic so that [`crate::models::ModelCalculation`] and
/// [`crate::msa::MSA`] can be parameterised over any alphabet without
/// duplicating logic.
pub trait AlphabetEncoding {
    /// Encoded symbol type. Implementors must be `#[repr(u8)]` (size 1, align 1)
    /// so that a `&[Self::Symbol]` can be reinterpreted as `&[u8]` by
    /// [`as_bytes`](AlphabetEncoding::as_bytes) for the byte-level distance kernel.
    type Symbol: Copy;

    /// Byte value of the gap symbol (`Self::Symbol::Gap as u8`). Used by the
    /// branchless / SIMD distance kernel for gap detection without enum matching.
    const GAP_BYTE: u8;

    /// Encodes a byte into the corresponding symbol in the alphabet.
    fn encode(symbol: u8) -> Self::Symbol;

    /// Total number of symbols in the alphabet (including gap and ambiguity).
    fn n_symbols() -> usize;

    /// Returns `true` if `symbol` represents a gap (`-`).
    fn is_gap(symbol: Self::Symbol) -> bool;

    /// Reinterprets an encoded sequence as a raw byte slice.
    ///
    /// The default implementation relies on the trait invariant that
    /// `Self::Symbol` is `#[repr(u8)]`, so a slice of symbols has the exact same
    /// layout as a `&[u8]` of equal length. All built-in alphabets uphold this.
    fn as_bytes(seq: &[Self::Symbol]) -> &[u8] {
        debug_assert_eq!(core::mem::size_of::<Self::Symbol>(), 1);
        // SAFETY: `Self::Symbol` is `#[repr(u8)]` (1 byte, align 1), so the slice
        // has identical layout to `&[u8]` of the same length.
        unsafe { core::slice::from_raw_parts(seq.as_ptr() as *const u8, seq.len()) }
    }
}

/// A single nucleotide in the DNA alphabet.
///
/// `N` is used for unknown/ambiguous bases. `Gap` corresponds to the `-`
/// character in aligned sequences. Any byte not matching `A/C/G/T/N/-`
/// (case-insensitive) is mapped to `N`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum DnaSymbol {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
    Gap = 5,
}

/// Marker struct for the DNA alphabet.
///
/// Pass this as the type parameter `A` when working with nucleotide sequences,
/// e.g. `MSA::<DNA>` or `pairwise_distance::<JukesCantor, DNA>`.
pub struct DNA;

impl AlphabetEncoding for DNA {
    type Symbol = DnaSymbol;

    const GAP_BYTE: u8 = DnaSymbol::Gap as u8;

    fn encode(symbol: u8) -> Self::Symbol {
        match symbol {
            b'A' | b'a' => DnaSymbol::A,
            b'C' | b'c' => DnaSymbol::C,
            b'G' | b'g' => DnaSymbol::G,
            b'T' | b't' => DnaSymbol::T,
            // U (uridine) is functionally equivalent to T in distance calculations.
            b'U' | b'u' => DnaSymbol::T,
            b'N' | b'n' => DnaSymbol::N,
            b'-' => DnaSymbol::Gap,
            // IUPAC ambiguity codes: R Y S W K M B D H V — treat as N.
            b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k'
            | b'M' | b'm' | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => {
                DnaSymbol::N
            }
            _ => DnaSymbol::N, // Treat unknowns as N
        }
    }

    fn n_symbols() -> usize {
        6 // A, C, G, T, N, Gap
    }

    fn is_gap(symbol: DnaSymbol) -> bool {
        symbol == DnaSymbol::Gap
    }
}

/// A single amino acid in the protein alphabet.
///
/// Covers the 20 standard amino acids plus `X` (unknown/ambiguous) and `Gap`
/// (the `-` alignment character). Encoding is case-insensitive; any byte that
/// does not match a known amino acid letter is mapped to `X`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum ProteinSymbol {
    A = 0,
    R = 1,
    N = 2,
    D = 3,
    C = 4,
    Q = 5,
    E = 6,
    G = 7,
    H = 8,
    I = 9,
    L = 10,
    K = 11,
    M = 12,
    F = 13,
    P = 14,
    S = 15,
    T = 16,
    W = 17,
    Y = 18,
    V = 19,
    X = 20,
    Gap = 21,
}

/// Marker struct for the protein alphabet.
///
/// Pass this as the type parameter `A` when working with amino acid sequences,
/// e.g. `MSA::<Protein>` or `pairwise_distance::<Poisson, Protein>`.
pub struct Protein;

impl AlphabetEncoding for Protein {
    type Symbol = ProteinSymbol;

    const GAP_BYTE: u8 = ProteinSymbol::Gap as u8;

    fn encode(symbol: u8) -> Self::Symbol {
        // Match on the uppercased byte so lowercase residues encode identically
        // (consistent with `DNA::encode` and `detect_alphabet`, which both
        // normalise case). The gap character `-` is unaffected by uppercasing.
        match symbol.to_ascii_uppercase() {
            b'A' => ProteinSymbol::A,
            b'R' => ProteinSymbol::R,
            b'N' => ProteinSymbol::N,
            b'D' => ProteinSymbol::D,
            b'C' => ProteinSymbol::C,
            b'Q' => ProteinSymbol::Q,
            b'E' => ProteinSymbol::E,
            b'G' => ProteinSymbol::G,
            b'H' => ProteinSymbol::H,
            b'I' => ProteinSymbol::I,
            b'L' => ProteinSymbol::L,
            b'K' => ProteinSymbol::K,
            b'M' => ProteinSymbol::M,
            b'F' => ProteinSymbol::F,
            b'P' => ProteinSymbol::P,
            b'S' => ProteinSymbol::S,
            b'T' => ProteinSymbol::T,
            b'W' => ProteinSymbol::W,
            b'Y' => ProteinSymbol::Y,
            b'V' => ProteinSymbol::V,
            b'X' => ProteinSymbol::X,
            b'-' => ProteinSymbol::Gap,
            _ => ProteinSymbol::X, // Treat unknowns as X
        }
    }

    fn n_symbols() -> usize {
        22 // 20 amino acids + X + Gap
    }

    fn is_gap(symbol: ProteinSymbol) -> bool {
        symbol == ProteinSymbol::Gap
    }
}

/// Runtime alphabet discriminant used by the auto-detection heuristic and
/// the `alphabet` override field in [`crate::config::NJConfig`] /
/// [`crate::config::DistConfig`].
///
/// [`crate::detect_alphabet`] inspects the raw sequence bytes and returns
/// `DNA` unless any character outside the DNA set is found, in which case it
/// returns `Protein`. The variant is then used to select the appropriate typed
/// code path.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, ts_rs::TS)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum Alphabet {
    DNA,
    Protein,
}

/// Encodes a string literal into a `Vec<DnaSymbol>` using [`DNA::encode`].
///
/// Intended for use in tests and benchmarks where typed symbol slices are
/// needed without verbose boilerplate. Accepts both uppercase and lowercase
/// characters; unknown bytes map to `DnaSymbol::N`.
///
/// ```ignore
/// let seq = dna!("ACGT-N");
/// assert_eq!(seq[4], DnaSymbol::Gap);
/// ```
#[cfg(test)]
macro_rules! dna {
    ($s:expr) => {
        $s.bytes()
            .map(|b| $crate::alphabet::DNA::encode(b))
            .collect::<Vec<$crate::alphabet::DnaSymbol>>()
    };
}
#[cfg(test)]
pub(crate) use dna;

/// Encodes a string literal into a `Vec<ProteinSymbol>` using [`Protein::encode`].
///
/// Intended for use in tests and benchmarks. Unknown bytes map to
/// `ProteinSymbol::X`.
///
/// ```ignore
/// let seq = protein!("ARND-");
/// assert_eq!(seq[4], ProteinSymbol::Gap);
/// ```
#[cfg(test)]
macro_rules! protein {
    ($s:expr) => {
        $s.bytes()
            .map(|b| $crate::alphabet::Protein::encode(b))
            .collect::<Vec<$crate::alphabet::ProteinSymbol>>()
    };
}
#[cfg(test)]
pub(crate) use protein;
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_encoding() {
        assert_eq!(DNA::encode(b'A'), DnaSymbol::A);
        assert_eq!(DNA::encode(b'c'), DnaSymbol::C);
        assert_eq!(DNA::encode(b'G'), DnaSymbol::G);
        assert_eq!(DNA::encode(b't'), DnaSymbol::T);
        assert_eq!(DNA::encode(b'N'), DnaSymbol::N);
        assert_eq!(DNA::encode(b'-'), DnaSymbol::Gap);
        assert_eq!(DNA::encode(b'X'), DnaSymbol::N); // Unknown
    }

    #[test]
    fn test_protein_encoding() {
        assert_eq!(Protein::encode(b'A'), ProteinSymbol::A);
        assert_eq!(Protein::encode(b'R'), ProteinSymbol::R);
        assert_eq!(Protein::encode(b'N'), ProteinSymbol::N);
        assert_eq!(Protein::encode(b'D'), ProteinSymbol::D);
        assert_eq!(Protein::encode(b'C'), ProteinSymbol::C);
        assert_eq!(Protein::encode(b'X'), ProteinSymbol::X);
        assert_eq!(Protein::encode(b'-'), ProteinSymbol::Gap);
        assert_eq!(Protein::encode(b'Z'), ProteinSymbol::X); // Unknown
    }

    #[test]
    fn test_protein_encoding_lowercase() {
        // Lowercase residues must encode identically to uppercase (consistent
        // with DNA). Regression test for silent data loss where lowercase
        // protein sequences were previously all mapped to X.
        assert_eq!(Protein::encode(b'a'), ProteinSymbol::A);
        assert_eq!(Protein::encode(b'k'), ProteinSymbol::K);
        assert_eq!(Protein::encode(b'y'), ProteinSymbol::Y);
        assert_eq!(Protein::encode(b'w'), ProteinSymbol::W);
        assert_eq!(Protein::encode(b'x'), ProteinSymbol::X);
    }

    // --- dna! macro ---

    #[test]
    fn test_dna_macro_produces_correct_symbols() {
        assert_eq!(
            dna!("ACGT"),
            vec![DnaSymbol::A, DnaSymbol::C, DnaSymbol::G, DnaSymbol::T]
        );
    }

    #[test]
    fn test_dna_macro_length() {
        assert_eq!(dna!("ACGTN-").len(), 6);
    }

    #[test]
    fn test_dna_macro_gap_and_n() {
        assert_eq!(dna!("-N"), vec![DnaSymbol::Gap, DnaSymbol::N]);
    }

    #[test]
    fn test_dna_macro_lowercase() {
        assert_eq!(
            dna!("acgt"),
            vec![DnaSymbol::A, DnaSymbol::C, DnaSymbol::G, DnaSymbol::T]
        );
    }

    #[test]
    fn test_dna_macro_empty() {
        assert_eq!(dna!(""), Vec::<DnaSymbol>::new());
    }

    // --- protein! macro ---

    #[test]
    fn test_protein_macro_produces_correct_symbols() {
        assert_eq!(
            protein!("ARND"),
            vec![
                ProteinSymbol::A,
                ProteinSymbol::R,
                ProteinSymbol::N,
                ProteinSymbol::D
            ]
        );
    }

    #[test]
    fn test_protein_macro_length() {
        assert_eq!(protein!("ACDEFGHIKLMNPQRSTVWY").len(), 20);
    }

    #[test]
    fn test_protein_macro_gap_and_unknown() {
        assert_eq!(protein!("-Z"), vec![ProteinSymbol::Gap, ProteinSymbol::X]);
    }

    #[test]
    fn test_protein_macro_lowercase() {
        assert_eq!(
            protein!("arnd"),
            vec![
                ProteinSymbol::A,
                ProteinSymbol::R,
                ProteinSymbol::N,
                ProteinSymbol::D
            ]
        );
    }

    #[test]
    fn test_protein_macro_empty() {
        assert_eq!(protein!(""), Vec::<ProteinSymbol>::new());
    }
}
