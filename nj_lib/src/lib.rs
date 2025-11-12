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
//! - NJNode
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
pub use crate::config::{FastaSequence, NJConfig, MSA};
use crate::tree::NJNode;
use bitvec::prelude::*;
use std::collections::HashMap;

/// Triangular distance matrix with node names
#[derive(Clone, Debug)]
pub struct TriMat {
    pub data: Vec<f64>,
    pub node_names: Vec<String>,
}

impl TriMat {
    pub fn with_names(names: Vec<String>) -> Self {
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

/// Neighbor Joining, iterator-based and borrow-safe
fn neighbor_joining(mut dist: TriMat) -> Result<NJNode, String> {
    let n = dist.dim();
    if n == 0 {
        return Err("Empty distance matrix".to_string());
    }
    if n == 1 {
        return Ok(NJNode::leaf(dist.node_names[0].clone()));
    }

    let mut active: BitVec<u8, Lsb0> = BitVec::repeat(true, n);
    let mut nodes: HashMap<usize, NJNode> = (0..n)
        .map(|i| (i, NJNode::leaf(dist.node_names[i].clone())))
        .collect();
    let mut row_sums: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| dist.get(i, j)).sum())
        .collect();
    let mut next_internal = n;

    for _ in 0..(n - 2) {
        let active_count = active.count_ones() as f64;
        let active_ref = &active;
        let row_sums_ref = &row_sums;
        let dist_ref = &dist;

        // Find minimal Q entry
        let pair_opt = (0..n)
            .filter(|&i| active_ref[i])
            .flat_map(|i| {
                (0..i).filter(move |&j| active_ref[j]).map(move |j| {
                    (
                        i,
                        j,
                        (active_count - 2.0) * dist_ref.get(i, j)
                            - row_sums_ref[i]
                            - row_sums_ref[j],
                        dist_ref.get(i, j),
                    )
                })
            })
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
            .map(|(i, j, _, d)| (i, j, d));

        let (i_min, j_min, dij) = match pair_opt {
            Some(t) => t,
            None => return Err("Failed to find a pair to join: no active pair found".to_string()),
        };

        // Clamp branch lengths to zero to avoid negative values due to non-ultrametric input or rounding errors
        let li = (0.5 * dij + (row_sums[i_min] - row_sums[j_min]) / (2.0 * (active_count - 2.0)))
            .max(0.0);
        let lj = (dij - li).max(0.0);
        let internal_name = format!("Node{}", next_internal);
        next_internal += 1;

        // Remove nodes safely
        let left_node = nodes
            .remove(&i_min)
            .ok_or_else(|| format!("Internal error: node {} missing during join", i_min))?;
        let right_node = nodes
            .remove(&j_min)
            .ok_or_else(|| format!("Internal error: node {} missing during join", j_min))?;
        let new_node = NJNode::internal(internal_name, left_node, right_node, li, lj);
        nodes.insert(i_min, new_node);
        active.set(j_min, false);

        // Update distances
        (0..n).filter(|&k| active[k] && k != i_min).for_each(|k| {
            let dik = dist.get(i_min, k);
            let djk = dist.get(j_min, k);
            let d_new = 0.5 * (dik + djk - dij);
            row_sums[i_min] = row_sums[i_min] - dik - djk + d_new;
            row_sums[k] = row_sums[k] - dik - djk + d_new;
            dist.set(i_min, k, d_new);
        });
        // Reset row sum for deactivated node to avoid affecting future calculations
        row_sums[j_min] = 0.0;
    }

    let remaining: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
    if remaining.len() != 2 {
        return Err(format!(
            "Expected 2 remaining active nodes, found {}. Input may be malformed.",
            remaining.len()
        ));
    }
    let i = remaining[0];
    let j = remaining[1];
    let dij = dist.get(i, j);
    let root_name = format!("Node{}", next_internal);
    let left = nodes
        .remove(&i)
        .ok_or_else(|| format!("Internal error: remaining node {} not found", i))?;
    let right = nodes
        .remove(&j)
        .ok_or_else(|| format!("Internal error: remaining node {} not found", j))?;
    Ok(NJNode::internal(
        root_name,
        left,
        right,
        dij / 2.0,
        dij / 2.0,
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
fn dist_from_msa(msa: &[FastaSequence]) -> TriMat {
    let n = msa.len();
    let names: Vec<String> = msa.into_iter().map(|fs| fs.header.clone()).collect();
    let mut dist = TriMat::with_names(names);

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

pub trait FastaReader {
    fn from_unnamed_sequences(sequences: Vec<String>) -> Self;
}

impl FastaReader for MSA {
    fn from_unnamed_sequences(sequences: Vec<String>) -> Self {
        sequences
            .into_iter()
            .enumerate()
            .map(|(i, s)| FastaSequence {
                header: format!("Seq{}", i),
                sequence: s,
            })
            .collect()
    }
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

    /// Recursively collects the names of all leaf nodes in the given NJNode tree.
    fn collect_leaf_names(node: &NJNode) -> Vec<String> {
        match (&node.left, &node.right) {
            (Some(left), Some(right)) => {
                let mut l = collect_leaf_names(left);
                let mut r = collect_leaf_names(right);
                l.append(&mut r);
                l
            }
            (None, None) => vec![node.name.clone()],
            _ => unreachable!("panic: NJNode should be either leaf or internal"),
        }
    }

    #[test]
    fn test_trimatrix_set_get_and_dim() {
        let mut m = TriMat::with_names(vec!["A".into(), "B".into(), "C".into()]);
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
    fn test_tri_from_msa_basic() {
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
    fn test_tri_from_msa_no_overlap() {
        // no overlapping non-gap positions between the two sequences -> distance 0.0
        let msa: MSA = vec![
            FastaSequence {
                header: "X".into(),
                sequence: "A--".into(),
            },
            FastaSequence {
                header: "Y".into(),
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
    fn test_neighbor_joining_two_taxa() {
        let mut m = TriMat::with_names(vec!["A".into(), "B".into()]);
        m.set(0, 1, 0.6);
        let tree = neighbor_joining(m).expect("NJ should succeed for two taxa");
        // root should be internal with two leaves A and B
        let leaves = collect_leaf_names(&tree);
        let mut leaves_sorted = leaves.clone();
        leaves_sorted.sort();
        assert_eq!(leaves_sorted, vec!["A".to_string(), "B".to_string()]);
        // branch lengths stored on the root node (left_len/right_len)
        assert!((tree.left_len - 0.3).abs() < 1e-12);
        assert!((tree.right_len - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_neighbor_joining_three_taxa_preserves_leaves() {
        // equidistant star with small distances
        let mut m = TriMat::with_names(vec!["A".into(), "B".into(), "C".into()]);
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
        fn check_nonneg(node: &NJNode) {
            assert!(node.left_len >= -1e-12, "left_len negative");
            assert!(node.right_len >= -1e-12, "right_len negative");
            if let (Some(l), Some(r)) = (&node.left, &node.right) {
                check_nonneg(l);
                check_nonneg(r);
            }
        }
        check_nonneg(&tree);
    }

    #[test]
    fn test_to_newick_produces_valid_format() {
        // Build simple tree manually: root with two leaves
        let left = NJNode::leaf("L".into());
        let right = NJNode::leaf("R".into());
        let root = NJNode::internal("root".into(), left, right, 0.1234, 0.5678);
        // to_newick hides internal names when requested
        let s_hidden = root.to_newick(true);
        // Expect parentheses and two branch lengths formatted to 3 decimals (per implementation)
        assert!(s_hidden.contains(":0.123"));
        assert!(s_hidden.contains(":0.568"));
        // with internal name shown, name appears at end
        let s_named = root.to_newick(false);
        assert!(s_named.ends_with("root"));
    }
}
