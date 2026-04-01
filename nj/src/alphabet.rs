//! Alphabet definitions and encoding for DNA and protein sequences.
//!
//! Each alphabet is represented by a marker struct ([`DNA`], [`Protein`]) that
//! implements [`AlphabetEncoding`]. The associated `Symbol` types are enums
//! whose variants cover the standard characters plus a gap and an ambiguity
//! symbol (`N` for DNA, `X` for protein). Unknown input bytes are silently
//! mapped to the ambiguity symbol rather than panicking.

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
    type Symbol: Copy;

    /// Encodes a byte into the corresponding symbol in the alphabet.
    fn encode(symbol: u8) -> Self::Symbol;

    /// Total number of symbols in the alphabet (including gap and ambiguity).
    fn n_symbols() -> usize;
}

/// A single nucleotide in the DNA alphabet.
///
/// `N` is used for unknown/ambiguous bases. `Gap` corresponds to the `-`
/// character in aligned sequences. Any byte not matching `A/C/G/T/N/-`
/// (case-insensitive) is mapped to `N`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DnaSymbol {
    A,
    C,
    G,
    T,
    N,
    Gap,
}

/// Marker struct for the DNA alphabet.
///
/// Pass this as the type parameter `A` when working with nucleotide sequences,
/// e.g. `MSA::<DNA>` or `pairwise_distance::<JukesCantor, DNA>`.
pub struct DNA;

impl AlphabetEncoding for DNA {
    type Symbol = DnaSymbol;

    fn encode(symbol: u8) -> Self::Symbol {
        match symbol {
            b'A' | b'a' => DnaSymbol::A,
            b'C' | b'c' => DnaSymbol::C,
            b'G' | b'g' => DnaSymbol::G,
            b'T' | b't' => DnaSymbol::T,
            b'N' | b'n' => DnaSymbol::N,
            b'-' => DnaSymbol::Gap,
            _ => DnaSymbol::N, // Treat unknowns as N
        }
    }

    fn n_symbols() -> usize {
        6 // A, C, G, T, N, Gap
    }
}

/// A single amino acid in the protein alphabet.
///
/// Covers the 20 standard amino acids plus `X` (unknown/ambiguous) and `Gap`
/// (the `-` alignment character). Any byte that does not match a known amino
/// acid letter is mapped to `X`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProteinSymbol {
    A,
    R,
    N,
    D,
    C,
    Q,
    E,
    G,
    H,
    I,
    L,
    K,
    M,
    F,
    P,
    S,
    T,
    W,
    Y,
    V,
    X,
    Gap,
}

/// Marker struct for the protein alphabet.
///
/// Pass this as the type parameter `A` when working with amino acid sequences,
/// e.g. `MSA::<Protein>` or `pairwise_distance::<Poisson, Protein>`.
pub struct Protein;

impl AlphabetEncoding for Protein {
    type Symbol = ProteinSymbol;

    fn encode(symbol: u8) -> Self::Symbol {
        match symbol {
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
}

/// Runtime alphabet discriminant used by the auto-detection heuristic.
///
/// [`crate::detect_alphabet`] inspects the raw sequence bytes and returns
/// `DNA` unless any character outside `A/C/G/T/U/N/-` is found, in which
/// case it returns `Protein`. The variant is then used to select the
/// appropriate typed code path.
#[derive(Debug, PartialEq)]
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
macro_rules! dna {
    ($s:expr) => {
        $s.bytes()
            .map(|b| $crate::alphabet::DNA::encode(b))
            .collect::<Vec<$crate::alphabet::DnaSymbol>>()
    };
}
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
macro_rules! protein {
    ($s:expr) => {
        $s.bytes()
            .map(|b| $crate::alphabet::Protein::encode(b))
            .collect::<Vec<$crate::alphabet::ProteinSymbol>>()
    };
}
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
    fn test_protein_macro_empty() {
        assert_eq!(protein!(""), Vec::<ProteinSymbol>::new());
    }
}
