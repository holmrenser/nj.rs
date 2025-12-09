use crate::tree::TreeNode;
use bitvec::prelude::{BitVec, Lsb0};

/// Triangular distance matrix with node names
#[derive(Clone, Debug)]
pub struct DistMat {
    data: Vec<f64>,
    names: Vec<String>,
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

    /// Performs the Neighbor-Joining algorithm on the distance matrix.
    /// Returns the resulting tree as a (recursive) TreeNode struct.
    pub fn neighbor_joining(mut self: DistMat) -> Result<TreeNode, String> {
        NJState::new(&mut self).run()
    }
}

/// Computes branch lengths for the new node joining i and j.
/// Clamps lengths to be non-negative.
fn compute_branch_lengths(
    d_ij: f64,
    row_sums: &[f64],
    i: usize,
    j: usize,
    active_count: usize,
) -> (f64, f64) {
    let n = active_count as f64;

    let li = (0.5 * d_ij + (row_sums[i] - row_sums[j]) / (2.0 * (n - 2.0))).max(0.0);
    let lj = (d_ij - li).max(0.0);

    (li, lj)
}
/// Internal state for the NJ algorithm.
/// Holds the distance matrix, active nodes, row sums, and current tree nodes.
struct NJState<'a> {
    dist: &'a mut DistMat,
    active: BitVec<u8, Lsb0>,
    row_sums: Vec<f64>,
    nodes: Vec<Option<TreeNode>>,
    next_internal: usize,
}

/// Implementation of the NJ algorithm state and methods.
impl<'a> NJState<'a> {
    pub fn new(dist: &'a mut DistMat) -> Self {
        let n = dist.dim();

        let nodes = dist
            .names
            .iter()
            .map(|name| Some(TreeNode::leaf(name.clone(), None)))
            .collect();

        let row_sums = (0..n)
            .map(|i| (0..n).map(|j| dist.get(i, j)).sum())
            .collect();

        NJState {
            dist,
            active: BitVec::repeat(true, n),
            row_sums,
            nodes,
            next_internal: n,
        }
    }

    /// Run NJ and return a TreeNode or Err.
    pub fn run(mut self) -> Result<TreeNode, String> {
        let n = self.dist.dim();
        if n == 0 {
            return Err("Empty distance matrix".into());
        }
        if n == 1 {
            return Ok(self.nodes[0].take().unwrap());
        }

        for _ in 0..(n - 2) {
            let (i, j, d_ij) = self.select_min_q_pair().ok_or("No pair found")?;

            let (li, lj) =
                compute_branch_lengths(d_ij, &self.row_sums, i, j, self.active.count_ones());

            self.join_nodes(i, j, li, lj);
            self.active.set(j, false);

            self.update_distances(i, j);
        }

        if self.active.count_ones() != 2 {
            return Err(format!(
                "Expected 2 active nodes but found {}",
                self.active.count_ones()
            ));
        }

        Ok(self.final_join())
    }

    /// Select pair i,j minimizing the Q-metric.
    fn select_min_q_pair(&self) -> Option<(usize, usize, f64)> {
        let n_active = self.active.count_ones() as f64;

        (0..self.dist.dim())
            .filter(|&i| self.active[i])
            .flat_map(|i| {
                (0..i).filter(move |&j| self.active[j]).map(move |j| {
                    let d_ij = self.dist.get(i, j);
                    let q = (n_active - 2.0) * d_ij - self.row_sums[i] - self.row_sums[j];
                    (i, j, q, d_ij)
                })
            })
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
            .map(|(i, j, _, d)| (i, j, d))
    }

    // Update distances after joining i (the new node) and removing j.
    fn update_distances(&mut self, i: usize, j: usize) {
        let d_ij = self.dist.get(i, j);

        // iterate ones returns indices of active bits — efficient
        for k in self.active.iter_ones() {
            if k == i {
                continue;
            }
            let d_ik = self.dist.get(i, k);
            let d_jk = self.dist.get(j, k);
            let d_new = 0.5 * (d_ik + d_jk - d_ij);

            // maintain row sums incrementally
            self.row_sums[i] += d_new - d_ik - d_jk;
            self.row_sums[k] += d_new - d_ik - d_jk;

            self.dist.set(i, k, d_new);
        }

        self.row_sums[j] = 0.0;
    }

    /// Join nodes i and j into a new internal node at index i.
    fn join_nodes(&mut self, i: usize, j: usize, li: f64, lj: f64) {
        let mut left = self.nodes[i].take().expect("node i exists");
        left.len = li;

        let mut right = self.nodes[j].take().expect("node j exists");
        right.len = lj;

        let name = format!("Node{}", self.next_internal);
        self.next_internal += 1;

        self.nodes[i] = Some(TreeNode::internal(
            name,
            Some([Box::new(left), Box::new(right)]),
            0.0,
        ));
    }

    /// Final join of the last two remaining active nodes.
    fn final_join(mut self) -> TreeNode {
        let mut it = self.active.iter_ones();
        let i = it.next().unwrap();
        let j = it.next().unwrap();

        let d_ij = self.dist.get(i, j);

        let mut left = self.nodes[i].take().unwrap();
        let mut right = self.nodes[j].take().unwrap();

        left.len = d_ij / 2.0;
        right.len = d_ij / 2.0;

        let name = format!("Node{}", self.next_internal);
        TreeNode::internal(name, Some([Box::new(left), Box::new(right)]), 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::FastaSequence;
    use crate::msa::MSA;

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
            None => vec![node.name.clone()],
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
        let msa = MSA::from_unnamed_sequences(seqs);
        let mat = msa.into_dist();
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
        let seqs = vec![
            FastaSequence {
                identifier: "X".into(),
                sequence: "A--".into(),
            },
            FastaSequence {
                identifier: "Y".into(),
                sequence: "--A".into(),
            },
        ];
        let msa = MSA::from_iter(seqs.into_iter());
        //let seqs: Ve c<String> = vec!["ACG".into(), "ATG".into(), "A-G".into()];
        //let msa = MSA::from_unnamed_sequences(seqs);
        let mat = msa.into_dist();
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
        assert_eq!(tree.name, "A");
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
        assert!((children[0].len - 0.3).abs() < 1e-12);
        assert!((children[1].len - 0.3).abs() < 1e-12);
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
                    assert!(children[0].len >= -1e-12, "left_len negative");
                    assert!(children[1].len >= -1e-12, "right_len negative");
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
        let left = TreeNode::leaf("L".into(), Some(0.1234));
        let right = TreeNode::leaf("R".into(), Some(0.5678));
        let root = TreeNode::internal("root".into(), Some([Box::new(left), Box::new(right)]), 0.0);
        // to_newick hides internal names when requested
        let s_hidden = root.to_newick(true);
        // Expect parentheses and two branch lengths formatted to 3 decimals (per implementation)
        assert!(s_hidden.contains(":0.123"));
        assert!(s_hidden.contains(":0.568"));
        // with internal name shown, name appears at end
        let s_named = root.to_newick(false);
        assert!(s_named.ends_with("root;"));
    }

    #[test]
    fn test_two_taxa_branch_length_sum() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 1.0);
        let tree = m.neighbor_joining().unwrap();
        let children = tree.children.unwrap();
        assert!((children[0].len + children[1].len - 1.0).abs() < 1e-12);
    }
}
