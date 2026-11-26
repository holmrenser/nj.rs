//! Pair-count tables: the sufficient statistic for the DNA distance models.
//!
//! For two aligned sequences, `f[a][b]` is the number of columns with base `a` in the
//! first and base `b` in the second sequence (`A=0, C=1, G=2, T=3`). Columns where
//! either sequence has anything else (a gap, `N`, an ambiguity code) are not counted:
//! pairwise deletion. Every DNA model is a function of this table, so the table can be
//! produced in different ways (a histogram over symbols here, bit-planes in
//! [`crate::prepared`]) without touching the model formulas.

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol};

/// A 4×4 table of aligned base pairs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DnaCounts {
    /// `f[a][b]`: columns with base `a` in the first and `b` in the second sequence.
    pub f: [[u32; 4]; 4],
}

impl DnaCounts {
    /// Reference construction: one histogram increment per column. Symbols with a byte
    /// value of 4 or more (`N`, gap) are skipped.
    pub fn from_bytes(a: &[u8], b: &[u8]) -> Self {
        let mut f = [[0u32; 4]; 4];
        for (&x, &y) in a.iter().zip(b) {
            if x < 4 && y < 4 {
                f[x as usize][y as usize] += 1;
            }
        }
        Self { f }
    }

    /// Reference construction from encoded DNA symbols.
    pub fn from_symbols(s1: &[DnaSymbol], s2: &[DnaSymbol]) -> Self {
        Self::from_bytes(DNA::as_bytes(s1), DNA::as_bytes(s2))
    }

    /// Number of columns where both sequences have a base.
    pub fn n_comparable(&self) -> usize {
        self.f.iter().flatten().map(|&c| c as usize).sum()
    }

    /// Number of comparable columns where the bases are equal.
    pub fn n_identical(&self) -> usize {
        (0..4).map(|k| self.f[k][k] as usize).sum()
    }

    /// Number of comparable columns where the bases differ.
    pub fn n_diff(&self) -> usize {
        self.n_comparable() - self.n_identical()
    }

    /// Transitions (A↔G, C↔T), counted in both orientations.
    pub fn transitions(&self) -> usize {
        let f = &self.f;
        (f[0][2] + f[2][0] + f[1][3] + f[3][1]) as usize
    }

    /// Transversions: every other difference.
    pub fn transversions(&self) -> usize {
        self.n_diff() - self.transitions()
    }

    /// Occurrences of each base over both sequences (so the four values sum to
    /// `2 · n_comparable`).
    pub fn base_counts(&self) -> [usize; 4] {
        let mut c = [0usize; 4];
        for a in 0..4 {
            for b in 0..4 {
                c[a] += self.f[a][b] as usize;
                c[b] += self.f[a][b] as usize;
            }
        }
        c
    }

    /// Columns where the two sequences have the (different) bases `a` and `b`, in either
    /// orientation.
    pub fn unordered(&self, a: usize, b: usize) -> usize {
        (self.f[a][b] + self.f[b][a]) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counts_skip_gaps_and_unknowns() {
        // A A C G | N - T T   vs   A G C T | A A N -
        let a = [0u8, 0, 1, 2, 4, 5, 3, 3];
        let b = [0u8, 2, 1, 3, 0, 0, 4, 5];
        let c = DnaCounts::from_bytes(&a, &b);
        assert_eq!(c.n_comparable(), 4);
        assert_eq!(c.n_identical(), 2); // A/A and C/C
        assert_eq!(c.n_diff(), 2); // A/G and G/T
        assert_eq!(c.transitions(), 1); // A/G
        assert_eq!(c.transversions(), 1); // G/T
        assert_eq!(c.base_counts(), [3, 2, 2, 1]);
        assert_eq!(c.unordered(0, 2), 1);
    }
}
