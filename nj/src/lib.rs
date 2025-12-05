//! Neighbor-Joining tree builder and triangular distance matrix utilities.
//!
//! This module provides a compact representation for pairwise distances (TriMat),
//! a simple tree node structure used to represent Neighbor-Joining (NJ) trees
//! (NJNode), a robust, borrow-safe implementation of the Neighbor-Joining
//! algorithm (neighbor_joining), and a helper to build a TriMat from a multiple
//! sequence alignment (tri_from_msa).
//!
//! # Types
//!
//! - TriMat
//!   - Compact triangular (upper-triangular stored as a flat Vec) distance
//!     matrix with associated node names.
//!   - Access via `get(i, j)` and `set(i, j, value)`. Diagonal distances are
//!     always treated as 0. The internal index mapping follows a strict
//!     triangular packing so only i != j is stored.
//!   - Construct with `TriMat::with_names(names)`; `dim()` returns the number
//!     of taxa.
//!
//! - TreeNode
//!   - Binary tree node representing either a leaf (no children) or an
//!     internal node with exactly two children.
//!   - Stores child branch lengths on the parent as `left_len` and `right_len`.
//!   - `to_newick(hide_internal)` emits a Newick string with branch lengths
//!     formatted to 3 decimal places. When `hide_internal` is true the internal
//!     node labels are omitted.
//!
//! # Functions
//!
//! - neighbor_joining(dist: TriMat) -> Result<NJNode, String>
//!   - Performs the Neighbor-Joining algorithm on the provided triangular
//!     distance matrix and returns the inferred tree as an `NJNode`.
//!   - The implementation is iterator-oriented and avoids borrowing issues by
//!     moving and reorganizing node ownership via a HashMap of indices -> nodes.
//!   - Uses a BitVec to track active taxa and maintains a `row_sums` vector
//!     for efficient Q-matrix computation.
//!   - Branch lengths computed for joins are clamped to >= 0.0 to avoid
//!     negative lengths caused by rounding or non-ultrametric inputs.
//!   - Returns an `Err` if the input matrix is empty or if an internal
//!     consistency check fails (e.g., no pair to join or unexpected number of
//!     remaining active nodes).
//!   - Complexity: O(n^2) memory for the distance representation and O(n^3)
//!     worst-case time if implemented naively; this implementation reduces
//!     repeated work via row sums but still iterates over active pairs while
//!     choosing minimal Q entries.
//!
//! - tri_from_msa(msa: &[String], names: Option<&[String]>) -> TriMat
//!   - Builds a `TriMat` from a multiple sequence alignment (MSA).
//!   - Pairwise distance = (number of differing positions) / (number of
//!     compared positions), ignoring any column where either sequence has a
//!     gap (`'-'`). If two sequences share no compared positions the distance
//!     is set to 0.0.
//!   - If `names` is omitted, default names are generated as `Seq0`, `Seq1`, ....
//!
//! # Usage notes and invariants
//!
//! - TriMat stores only i != j entries; `get(i,i)` is defined to be 0.0.
//! - Indices used with TriMat must be valid (0..dim()). `set` silently ignores
//!   attempts to set diagonal entries.
//! - neighbor_joining consumes the distance matrix (it is taken by value) and
//!   returns a tree whose leaves preserve the original node name strings.
//! - The NJ implementation expects a symmetric distance matrix; asymmetries
//!   will lead to undefined/incorrect results.
//! - Branch lengths are clamped to non-negative values to keep the tree valid
//!   in downstream tools; if you require different handling of negative lengths
//!   adjust the implementation accordingly.
//!
//! # Examples
//!
//! Construct a TriMat from an MSA and produce a Newick tree:
//!
//! ```ignore
//! let msa = vec!["ACG".into(), "ATG".into(), "A-G".into()];
//! let tm = tri_from_msa(&msa, None);
//! let tree = neighbor_joining(tm).expect("NJ failed");
//! let newick = tree.to_newick(true);
//! ```
//!
//! # Errors & Panics
//!
//! - Functions in this module return `Result` where appropriate to indicate
//!   recoverable errors (e.g., empty input). Internal indexing errors return
//!   descriptive strings rather than panicking.
//! - The code assumes reasonably small `n` such that triangular storage in a
//!   Vec is practical; very large numbers of taxa may cause memory pressure.
//!
//! # Testing
//!
//! The accompanying unit tests exercise TriMat indexing, distance-building
//! behavior with gaps, two- and three-taxon NJ behavior, and Newick output
//! formatting. They also validate that branch lengths are non-negative after
//! clamping.
mod config;
mod tree;
pub use crate::config::{FastaSequence, MSA, NJConfig};
use crate::tree::TreeNode;
use bitvec::prelude::*;

/// Triangular distance matrix with node names
#[derive(Clone, Debug)]
pub struct DistMat {
    pub data: Vec<f64>,
    pub node_names: Vec<String>,
}

/// (Lower) Triangular distance matrix implementation
/// Stores only the i > j entries in a flat Vec for memory efficiency.
impl DistMat {
    pub fn empty_with_names(names: Vec<String>) -> Self {
        let n = names.len();
        Self {
            data: vec![0.0; n * (n - 1) / 2],
            node_names: names,
        }
    }

    pub fn dim(&self) -> usize {
        self.node_names.len()
    }

    fn idx(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i > j { (i, j) } else { (j, i) };
        i * (i - 1) / 2 + j
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else {
            self.data[self.idx(i, j)]
        }
    }

    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        if i != j {
            let index = self.idx(i, j);
            self.data[index] = val;
        }
    }
}
/// Selects the pair (i, j) with the minimum Q(i, j) value from the active set.
/// Returns (i, j, d_ij) where d_ij is the distance between i and j.
fn select_min_q_pair(
    dist: &DistMat,
    active: &BitVec<u8, Lsb0>,
    row_sums: &[f64],
) -> Option<(usize, usize, f64)> {
    let n_active = active.count_ones() as f64;

    (0..dist.dim())
        .filter(|&i| active[i])
        .flat_map(|i| {
            (0..i).filter(move |&j| active[j]).map(move |j| {
                let d_ij = dist.get(i, j);
                let q_ij = (n_active - 2.0) * d_ij - row_sums[i] - row_sums[j];
                (i, j, q_ij, d_ij)
            })
        })
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        .map(|(i, j, _, d)| (i, j, d))
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

/// Joins nodes i and j into a new internal node at index i.
/// Sets the branch lengths for the children.
/// Removes node j by setting it to None.
/// The new internal node has the given name.
fn join_nodes(nodes: &mut [Option<TreeNode>], i: usize, j: usize, li: f64, lj: f64, name: String) {
    let mut left = nodes[i].take().unwrap();
    left.len = li;

    let mut right = nodes[j].take().unwrap();
    right.len = lj;

    nodes[i] = Some(TreeNode::internal(
        name,
        Some([Box::new(left), Box::new(right)]),
        0.0,
    ));
}

/// Updates the distance matrix and row sums after joining nodes i and j.
/// The new node is at index i; index j is removed.
/// The distances are updated according to the NJ formula.
/// Row sums are adjusted accordingly.
fn update_distances(
    dist: &mut DistMat,
    row_sums: &mut [f64],
    active: &BitVec<u8, Lsb0>,
    i: usize,
    j: usize,
) {
    let d_ij = dist.get(i, j);

    for k in active.iter_ones() {
        if k == i {
            continue;
        }
        let d_ik = dist.get(i, k);
        let d_jk = dist.get(j, k);
        let d_new = 0.5 * (d_ik + d_jk - d_ij);

        row_sums[i] += d_new - d_ik - d_jk;
        row_sums[k] += d_new - d_ik - d_jk;

        dist.set(i, k, d_new);
    }

    row_sums[j] = 0.0;
}

/// Final join of the last two active nodes into a root node.
/// Sets branch lengths accordingly.
/// The new root node has the given name.
fn final_join(
    nodes: &mut [Option<TreeNode>],
    dist: &DistMat,
    active: &BitVec<u8, Lsb0>,
    name: String,
) -> TreeNode {
    let mut iter = active.iter_ones();
    let i = iter.next().unwrap();
    let j = iter.next().unwrap();

    let d_ij = dist.get(i, j);

    let mut left = nodes[i].take().unwrap();
    let mut right = nodes[j].take().unwrap();

    left.len = d_ij / 2.0;
    right.len = d_ij / 2.0;

    TreeNode::internal(name, Some([Box::new(left), Box::new(right)]), 0.0)
}

/// Neighbor Joining, iterator-based and borrow-safe
/// Performs the Neighbor-Joining algorithm on the provided triangular distance matrix
/// and returns the inferred tree as a TreeNode.
/// /// # Arguments
/// * `dist` - A triangular distance matrix representing pairwise distances between taxa.
/// # Returns
/// A Result containing the root TreeNode of the inferred NJ tree, or an error
/// string if the input is invalid or the algorithm fails.
/// # Errors
/// Returns an error if the input distance matrix is empty or if internal
/// consistency checks fail (e.g., no pair to join, unexpected number of
/// remaining active nodes).
fn neighbor_joining(mut dist: DistMat) -> Result<TreeNode, String> {
    // Number of genes/proteins in the distance matrix
    let n = dist.dim();
    if n == 0 {
        return Err("Empty distance matrix".into());
    }
    if n == 1 {
        return Ok(TreeNode::leaf(dist.node_names[0].clone(), None));
    }

    // Bit mask indicating which elements are currently active.
    let mut active: BitVec<u8, Lsb0> = BitVec::repeat(true, n);

    let mut nodes: Vec<Option<TreeNode>> = dist
        .node_names
        .iter()
        .map(|name| Some(TreeNode::leaf(name.clone(), None)))
        .collect();

    let mut row_sums: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| dist.get(i, j)).sum())
        .collect();

    let mut next_internal = n;

    for _ in 0..(n - 2) {
        let (i, j, d_ij) = select_min_q_pair(&dist, &active, &row_sums).ok_or("No pair found")?;

        let (l_i, l_j) =
            compute_branch_lengths(d_ij, &row_sums, i, j, active.count_ones() as usize);

        let name = format!("Node{}", next_internal);
        next_internal += 1;

        join_nodes(&mut nodes, i, j, l_i, l_j, name.clone());
        active.set(j, false);

        update_distances(&mut dist, &mut row_sums, &active, i, j);
    }

    if active.count_ones() != 2 {
        return Err(format!(
            "Expected 2 remaining active nodes, found {}. Input may be malformed.",
            active.count_ones()
        ));
    }

    Ok(final_join(
        &mut nodes,
        &dist,
        &active,
        format!("Node{}", next_internal),
    ))
}

/// Build TriMat from MSA in iterator style
/// Constructs a triangular distance matrix (`TriMat`) from a multiple sequence alignment (MSA).
///
/// # Arguments
/// * `msa` - A slice of strings, each representing a sequence in the alignment. All sequences should be of equal length.
/// * `names` - An optional slice of names for the sequences. If not provided, default names "Seq0", "Seq1", ... are used.
///
/// # Behavior
/// Pairwise distances are calculated as the proportion of differing, non-gap (`'-'`) positions between sequences.
/// Positions where either sequence has a gap are ignored in the calculation.
/// If no valid positions exist between a pair, the distance is set to 0.0.
fn dist_from_msa(msa: &[FastaSequence]) -> DistMat {
    let n = msa.len();
    let names: Vec<String> = msa.into_iter().map(|fs| fs.identifier.clone()).collect();
    let mut dist = DistMat::empty_with_names(names);

    (0..n).for_each(|i| {
        (0..i).for_each(|j| {
            let (diffs, valid) = msa[i]
                .sequence
                .chars()
                .zip(msa[j].sequence.chars())
                .filter(|(a, b)| *a != '-' && *b != '-')
                .fold((0, 0), |(d, v), (a, b)| {
                    (d + if a != b { 1 } else { 0 }, v + 1)
                });
            let d = if valid > 0 {
                diffs as f64 / valid as f64
            } else {
                0.0
            };
            dist.set(i, j, d);
        })
    });
    dist
}

pub fn nj(conf: NJConfig) -> Result<String, String> {
    let dist = dist_from_msa(&conf.msa);
    let tree = neighbor_joining(dist).map_err(|e| format!("neighbor-joining failed: {e}"))?;
    let newick = format!("{};", tree.to_newick(conf.hide_internal));
    Ok(newick)
}

#[cfg(test)]
mod tests {
    use super::*;

    pub trait FastaReader {
        fn from_unnamed_sequences(sequences: Vec<String>) -> Self;
    }

    impl FastaReader for MSA {
        fn from_unnamed_sequences(sequences: Vec<String>) -> Self {
            sequences
                .into_iter()
                .enumerate()
                .map(|(i, s)| FastaSequence {
                    identifier: format!("Seq{}", i),
                    sequence: s,
                })
                .collect()
        }
    }

    /// Recursively collects the names of all leaf nodes in the given NJNode tree.
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
    fn test_distmat_set_get_and_dim() {
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
        let mat = dist_from_msa(&msa);
        // names default Seq0, Seq1, Seq2
        assert_eq!(mat.node_names, vec!["Seq0", "Seq1", "Seq2"]);
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
        let msa: MSA = vec![
            FastaSequence {
                identifier: "X".into(),
                sequence: "A--".into(),
            },
            FastaSequence {
                identifier: "Y".into(),
                sequence: "--A".into(),
            },
        ];
        //let seqs: Ve c<String> = vec!["ACG".into(), "ATG".into(), "A-G".into()];
        //let msa = MSA::from_unnamed_sequences(seqs);
        let mat = dist_from_msa(&msa);
        assert_eq!(mat.node_names, vec!["X", "Y"]);
        assert!((mat.get(0, 1) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_distmat_symmetry_invariant() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into(), "C".into()]);
        m.set(2, 0, 0.25);
        assert_eq!(m.get(0, 2), 0.25);
        assert_eq!(m.get(2, 0), 0.25);
    }

    #[test]
    fn test_neighbor_joining_one_taxon() {
        let m = DistMat::empty_with_names(vec!["A".into()]);
        let tree = neighbor_joining(m).expect("NJ should succeed for one taxon");
        assert!(tree.children.is_none());
        assert_eq!(tree.name, "A");
    }

    #[test]
    fn test_neighbor_joining_two_taxa() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 0.6);
        let tree = neighbor_joining(m).expect("NJ should succeed for two taxa");
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
        let tree = neighbor_joining(m).expect("NJ should succeed for three taxa");
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
        assert!(s_named.ends_with("root"));
    }

    #[test]
    fn test_nj_wrapper_adds_semicolon() {
        let msa = MSA::from_unnamed_sequences(vec!["A".into(), "A".into()]);
        let conf = NJConfig {
            msa,
            hide_internal: true,
        };
        let out = nj(conf).unwrap();
        assert!(out.ends_with(';'));
    }

    #[test]
    fn test_nj_deterministic_order() {
        let msa = MSA::from_unnamed_sequences(vec!["ACG".into(), "ATG".into(), "AGG".into()]);
        let conf = NJConfig {
            msa,
            hide_internal: true,
        };

        let t1 = nj(conf.clone()).unwrap();
        let t2 = nj(conf).unwrap();
        assert_eq!(t1, t2);
    }

    #[test]
    fn test_two_taxa_branch_length_sum() {
        let mut m = DistMat::empty_with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 1.0);
        let tree = neighbor_joining(m).unwrap();
        let children = tree.children.unwrap();
        assert!((children[0].len + children[1].len - 1.0).abs() < 1e-12);
    }
}
