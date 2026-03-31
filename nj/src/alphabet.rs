pub trait AlphabetEncoding {
    type Symbol: Copy;

    /// Encodes a byte into the corresponding symbol in the alphabet.
    fn encode(symbol: u8) -> Self::Symbol;

    /// Total number of symbols in the alphabet.
    fn n_symbols() -> usize;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DnaSymbol {
    A,
    C,
    G,
    T,
    N,
    Gap,
}

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

#[derive(Debug, PartialEq)]
pub enum Alphabet {
    DNA,
    Protein,
}
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
}
