//! Lower-triangular pairwise distance matrix and entry point for NJ.
//!
//! [`DistMat`] stores only the `n*(n-1)/2` strictly lower-triangular entries
//! in a flat `Vec<f64>`, keeping memory at O(n²/2). Access and mutation via
//! [`get`](DistMat::get) / [`set`](DistMat::set) are symmetric, so callers
//! never need to worry about index ordering. Build from an [`MSA`] via
//! [`DistMat::from_msa`] (or the [`MSA::into_dist`] convenience wrapper), then
//! call [`neighbor_joining`](DistMat::neighbor_joining) to run the algorithm.

use crate::MSA;
use crate::alphabet::AlphabetEncoding;
use crate::models::{ModelCalculation, pairwise_distance};
use crate::nj::NJState;
use crate::tree::TreeNode;

/// Lower-triangular pairwise distance matrix.
///
/// For `n` taxa, only the `n*(n-1)/2` strictly lower-triangular entries
/// (`i > j`) are stored in a flat `Vec<f64>`, in row-major order. Diagonal
/// entries (`d(i,i)`) are always treated as `0.0` and are never stored.
/// Symmetry is maintained by [`get`](DistMat::get) and
/// [`set`](DistMat::set): accessing or writing `(i, j)` and `(j, i)` always
/// refer to the same stored element.
#[derive(Clone, Debug)]
pub struct DistMat {
    /// Flat storage for the strictly lower-triangular entries, packed in
    /// row-major order: index `i*(i-1)/2 + j` for `i > j`.
    pub data: Vec<f64>,
    /// Sequence/taxa names, one per row/column.
    pub names: Vec<String>,
}

impl DistMat {
    /// Creates a zeroed distance matrix for `names.len()` taxa.
    pub fn empty_with_names(names: Vec<String>) -> Self {
        let n = names.len();
        Self {
            data: vec![0.0; n.saturating_mul(n.saturating_sub(1)) / 2],
            names,
        }
    }

    /// Returns the dimension (number of nodes) of the distance matrix.
    pub fn dim(&self) -> usize {
        self.names.len()
    }

    /// Returns the flat-vector index for the `(i, j)` pair.
    ///
    /// Normalises so that the larger index is used as the row:
    /// `index = max(i,j) * (max(i,j) - 1) / 2 + min(i,j)`.
    /// The caller is responsible for ensuring `i != j`.
    fn idx(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i > j { (i, j) } else { (j, i) };
        i * (i - 1) / 2 + j
    }

    /// Gets the distance between nodes i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else {
            self.data[self.idx(i, j)]
        }
    }

    /// Sets the distance between nodes i and j.
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        if i != j {
            let index = self.idx(i, j);
            self.data[index] = val;
        }
    }

    /// Computes a pairwise distance matrix from an MSA using substitution model `M`.
    ///
    /// Iterates over all `n*(n-1)/2` sequence pairs and fills the lower triangle
    /// using [`pairwise_distance`]. Prefer [`MSA::into_dist`] over calling this
    /// directly.
    pub fn from_msa<M, A>(msa: &MSA<A>) -> DistMat
    where
        M: ModelCalculation<A>,
        A: AlphabetEncoding,
    {
        let n = msa.n_sequences;
        let mut dist = DistMat::empty_with_names(msa.identifiers.clone());

        for i in 0..n {
            let s1 = &msa.sequences[i];
            for j in 0..i {
                let s2 = &msa.sequences[j];
                let d = pairwise_distance::<M, A>(s1, s2);
                dist.set(i, j, d);
            }
        }
        dist
    }

    /// Runs the Neighbor-Joining algorithm and returns the inferred tree.
    ///
    /// Consumes `self` because the NJ state machine mutates the matrix in
    /// place during the join iterations. Returns `Err` for an empty matrix
    /// or if an internal consistency check fails.
    pub fn neighbor_joining(mut self: DistMat) -> Result<TreeNode, String> {
        NJState::new(&mut self).run()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::DNA;
    use crate::models::PDiff;
    use crate::msa::MSA;
    use crate::tree::NameOrSupport;

    /// Recursively collects the names of all leaf nodes in the given NJNode tree.
    /// Used for testing correctness of tree structure.
    fn collect_leaf_names(node: &TreeNode) -> Vec<String> {
        match &node.children {
            Some([left, right]) => {
                let mut l = collect_leaf_names(left);
                let mut r = collect_leaf_names(right);
                l.append(&mut r);
                l
            }
            None => match node.label {
                Some(NameOrSupport::Name(ref name)) => vec![name.clone()],
                _ => vec![node.identifier.to_string()],
            },
        }
    }

    #[test]
    fn test_dist_set_get_and_dim() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into(), "C".into()]);
        assert_eq!(m.dim(), 3);
        // initially zero
        assert_eq!(m.get(0, 1), 0.0);
        // set symmetric entries via set
        m.set(0, 1, 0.25);
        m.set(2, 1, 0.75);
        assert!((m.get(0, 1) - 0.25).abs() < 1e-12);
        assert!((m.get(1, 0) - 0.25).abs() < 1e-12);
        assert!((m.get(2, 1) - 0.75).abs() < 1e-12);
        assert!((m.get(1, 2) - 0.75).abs() < 1e-12);
        // diagonal always zero
        assert_eq!(m.get(1, 1), 0.0);
    }

    #[test]
    fn test_dist_empty_with_names_zero_dim() {
        let m = DistMat::empty_with_names(vec![]);
        assert_eq!(m.dim(), 0);
        assert!(m.data.is_empty());
    }

    #[test]
    fn test_dist_from_msa_basic() {
        // seq0 vs seq1 differ at middle position only
        let seqs: Vec<String> = vec!["ACG".into(), "ATG".into(), "A-G".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(seqs).unwrap();
        let mat = msa.into_dist::<PDiff>();
        // names default Seq0, Seq1, Seq2
        assert_eq!(mat.names, vec!["Seq0", "Seq1", "Seq2"]);
        // seq0 vs seq1: one mismatch across three aligned positions -> 1/3
        assert!((mat.get(0, 1) - (1.0 / 3.0)).abs() < 1e-12);
        // seq0 vs seq2: the gap position is ignored for mismatches and the full alignment length is still the denominator
        assert!((mat.get(0, 2) - 0.0).abs() < 1e-12);
        // seq1 vs seq2: same behavior, no non-gap mismatches so the distance is 0.0
        assert!((mat.get(1, 2) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dist_from_msa_gap_positions_do_not_change_denominator() {
        let seqs: Vec<String> = vec!["AT-".into(), "ACG".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(seqs).unwrap();
        let mat = msa.into_dist::<PDiff>();

        assert!((mat.get(0, 1) - (1.0 / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn test_dist_from_msa_no_overlap() {
        // no non-gap mismatches are counted, so the distance remains 0.0
        let seqs = vec![("X".into(), "A--".into()), ("Y".into(), "--A".into())];
        let msa = MSA::<DNA>::from_iter(seqs.into_iter());
        let mat = msa.into_dist::<PDiff>();
        assert_eq!(mat.names, vec!["X", "Y"]);
        assert!((mat.get(0, 1) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dist_symmetry_invariant() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into(), "C".into()]);
        m.set(2, 0, 0.25);
        assert_eq!(m.get(0, 2), 0.25);
        assert_eq!(m.get(2, 0), 0.25);
    }

    #[test]
    fn test_neighbor_joining_one_taxon() {
        let m = DistMat::empty_with_names(vec!["A".into()]);
        let tree = m
            .neighbor_joining()
            .expect("NJ should succeed for one taxon");
        assert!(tree.children.is_none());
        assert!(matches!(tree.label, Some(NameOrSupport::Name(ref s)) if s == "A"));
    }

    #[test]
    fn test_neighbor_joining_empty_matrix_returns_error() {
        let m = DistMat::empty_with_names(vec![]);
        let err = m.neighbor_joining().unwrap_err();
        assert_eq!(err, "Empty distance matrix");
    }

    #[test]
    fn test_neighbor_joining_two_taxa() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 0.6);
        let tree = m
            .neighbor_joining()
            .expect("NJ should succeed for two taxa");
        // root should be internal with two leaves A and B
        let leaves = collect_leaf_names(&tree);
        let mut leaves_sorted = leaves.clone();
        leaves_sorted.sort();
        assert_eq!(leaves_sorted, vec!["A".to_string(), "B".to_string()]);
        // branch lengths stored on child nodes
        assert!(tree.children.is_some());
        let children = tree.children.as_ref().unwrap();
        assert!((children[0].len.unwrap() - 0.3).abs() < 1e-12);
        assert!((children[1].len.unwrap() - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_neighbor_joining_three_taxa_preserves_leaves() {
        // equidistant star with small distances
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into(), "C".into()]);
        m.set(0, 1, 0.2);
        m.set(0, 2, 0.2);
        m.set(1, 2, 0.2);
        let tree = m
            .neighbor_joining()
            .expect("NJ should succeed for three taxa");
        let mut leaves = collect_leaf_names(&tree);
        leaves.sort();
        assert_eq!(
            leaves,
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
        // ensure non-negative branch lengths in the produced tree (clamping is applied in implementation)
        fn check_nonneg(node: &TreeNode) {
            match &node.children {
                None => return,
                Some(children) => {
                    assert!(children[0].len.unwrap() >= -1e-12, "left_len negative");
                    assert!(children[1].len.unwrap() >= -1e-12, "right_len negative");
                    check_nonneg(children[0].as_ref());
                    check_nonneg(children[1].as_ref());
                }
            }
        }
        check_nonneg(&tree);
    }

    #[test]
    fn test_to_newick_produces_valid_format() {
        // Build simple tree manually: root with two leaves
        let left = TreeNode::leaf(0, "L".into(), Some(0.1234));
        let right = TreeNode::leaf(1, "R".into(), Some(0.5678));
        let root = TreeNode::internal(
            2,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            Some(85),
        );
        let newick = root.to_newick();
        // Expect parentheses and two branch lengths formatted to 3 decimals (per implementation)
        assert!(newick.contains(":0.123"));
        assert!(newick.contains(":0.568"));
        assert!(newick.starts_with("("));
        assert!(newick.ends_with(")85;"));
    }

    #[test]
    fn test_two_taxa_branch_length_sum() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 1.0);
        let tree = m.neighbor_joining().unwrap();
        let children = tree.children.unwrap();
        assert!((children[0].len.unwrap() + children[1].len.unwrap() - 1.0).abs() < 1e-12);
    }
}
