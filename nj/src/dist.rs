use crate::MSA;
use crate::alphabet::AlphabetEncoding;
use crate::models::{ModelCalculation, pairwise_distance};
use crate::nj::NJState;
use crate::tree::TreeNode;

/// Triangular distance matrix with node names
#[derive(Clone, Debug)]
pub struct DistMat {
    pub data: Vec<f64>,
    pub names: Vec<String>,
}

/// (Lower) Triangular distance matrix implementation
/// Stores only the i > j entries in a flat Vec for memory efficiency.
impl DistMat {
    /// Creates an empty distance matrix with the given node names.
    pub fn empty_with_names(names: Vec<String>) -> Self {
        let n = names.len();
        Self {
            data: vec![0.0; n * (n - 1) / 2],
            names,
        }
    }

    /// Returns the dimension (number of nodes) of the distance matrix.
    pub fn dim(&self) -> usize {
        self.names.len()
    }

    /// Computes the internal index in the flat Vec for the (i, j) entry.
    /// Assumes i != j. Returns the index for the stored entry.
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

    /// Performs the Neighbor-Joining algorithm on the distance matrix.
    /// Returns the resulting tree as a (recursive) TreeNode struct.
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
    fn test_dist_from_msa_basic() {
        // seq0 vs seq1 differ at middle position only
        let seqs: Vec<String> = vec!["ACG".into(), "ATG".into(), "A-G".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(seqs);
        let mat = msa.into_dist::<PDiff>();
        // names default Seq0, Seq1, Seq2
        assert_eq!(mat.names, vec!["Seq0", "Seq1", "Seq2"]);
        // seq0 vs seq1: positions considered all three -> diffs 1/3
        assert!((mat.get(0, 1) - (1.0 / 3.0)).abs() < 1e-12);
        // seq0 vs seq2: position considered 0 and 2 -> diffs 0/2
        assert!((mat.get(0, 2) - 0.0).abs() < 1e-12);
        // seq1 vs seq2: compare positions 0 and 2 -> both are identical ('A' vs 'A', 'G' vs 'G') -> 0.0
        assert!((mat.get(1, 2) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dist_from_msa_no_overlap() {
        // no overlapping non-gap positions between the two sequences -> distance 0.0
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
