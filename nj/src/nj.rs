//! Core Neighbor-Joining algorithm implementation.
//!
//! The algorithm runs in three phases:
//! 1. **Initialisation** ([`NJState::new`]) — build leaf nodes, precompute
//!    row sums of the distance matrix, and mark all nodes active.
//! 2. **Iteration** ([`NJState::run`]) — repeat `n-2` times: select the pair
//!    `(i, j)` minimising the Q-criterion, join them into a new internal node,
//!    update the distance matrix, and deactivate `j`.
//! 3. **Final join** ([`NJState::final_join`]) — connect the two remaining
//!    active nodes to form the root, splitting the remaining distance equally.
//!
//! The distance matrix is mutated in-place during iteration, so [`NJState`]
//! holds a mutable reference. All active node tracking is done via a
//! [`bitvec::prelude::BitVec`] for compact storage.

use bitvec::prelude::{BitVec, Lsb0};

use crate::{distance_matrix::DistMat, tree::TreeNode};

/// Computes the branch lengths from the new internal node to taxa `i` and `j`.
///
/// Uses the standard NJ formula:
/// ```text
/// l_i = 0.5 * d(i,j) + (S_i - S_j) / (2 * (n - 2))
/// l_j = d(i,j) - l_i
/// ```
/// where `S_i` and `S_j` are the row sums of the Q-reduced distance matrix
/// and `n` is the number of currently active taxa. Both lengths are clamped
/// to `0.0` to prevent negative branch lengths from rounding or
/// non-ultrametric inputs.
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
/// Mutable state for a single NJ run.
///
/// Owns the partial tree nodes, tracks which taxa are still active, and
/// maintains incremental row sums so the Q-criterion can be evaluated in
/// O(n²) per iteration rather than recomputing from scratch.
pub struct NJState<'a> {
    dist: &'a mut DistMat,
    active: BitVec<u8, Lsb0>,
    row_sums: Vec<f64>,
    nodes: Vec<Option<TreeNode>>,
    next_internal: usize,
}

impl<'a> NJState<'a> {
    /// Initialises the NJ state from a distance matrix.
    ///
    /// Creates one leaf [`TreeNode`] per taxon, precomputes the full row sums
    /// `S_i = Σ_j d(i,j)`, and marks all `n` taxa as active.
    pub fn new(dist: &'a mut DistMat) -> Self {
        let n = dist.dim();

        let nodes = dist
            .names
            .iter()
            .enumerate()
            .map(|(i, name)| Some(TreeNode::leaf(i, name.clone(), None)))
            .collect();

        let row_sums = (0..n)
            .map(|i| (0..n).map(|j| dist.get(i, j)).sum())
            .collect();

        let active = BitVec::repeat(true, n);

        NJState {
            dist,
            active,
            row_sums,
            nodes,
            next_internal: n, // internal nodes get identifiers starting after the leaves
        }
    }

    /// Runs the NJ algorithm to completion and returns the root [`TreeNode`].
    ///
    /// Returns `Err` if the distance matrix is empty, if no pair can be found
    /// during an iteration (should be unreachable with a valid matrix), or if
    /// the number of active nodes after the main loop is not exactly 2.
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
            // Must precede marking j inactive: update_distances iterates active nodes.
            self.update_distances(i, j);
            self.active.set(j, false);
        }

        if self.active.count_ones() != 2 {
            return Err(format!(
                "Expected 2 active nodes but found {}",
                self.active.count_ones()
            ));
        }

        self.final_join()
    }

    /// Selects the active pair `(i, j)` with the smallest Q value.
    ///
    /// Q(i,j) = (n-2)·d(i,j) - S_i - S_j, where `n` is the current number
    /// of active taxa and `S_i` is the precomputed row sum. Returns the pair
    /// indices and the raw distance `d(i,j)` (not the Q value). NaN Q values
    /// (which can arise when distances are infinite) are treated as the worst
    /// case via `unwrap_or(Ordering::Greater)`.
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
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Greater))
            .map(|(i, j, _, d)| (i, j, d))
    }

    /// Merges nodes `i` and `j` into a new internal node stored at index `i`.
    ///
    /// Assigns branch lengths `li` and `lj` to the child nodes, then replaces
    /// `nodes[i]` with a new internal [`TreeNode`]. Node `j` is consumed here
    /// and its slot is left empty; it will be marked inactive by the caller.
    fn join_nodes(&mut self, i: usize, j: usize, li: f64, lj: f64) {
        let mut left = self.nodes[i].take().expect("node i exists");
        left.len = Some(li);

        let mut right = self.nodes[j].take().expect("node j exists");
        right.len = Some(lj);

        self.nodes[i] = Some(TreeNode::internal(
            self.next_internal,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            None,
        ));
        self.next_internal += 1;
    }

    /// Updates the distance matrix and row sums after joining `i` and `j`.
    ///
    /// For every remaining active taxon `k` (excluding `i` and `j`), applies
    /// the classic NJ distance update formula:
    /// ```text
    /// d(new, k) = 0.5 * (d(i,k) + d(j,k) - d(i,j))
    /// ```
    /// Row sums are maintained incrementally so they need not be recomputed
    /// from scratch each iteration. Must be called **before** marking `j`
    /// inactive, since the loop iterates over currently active nodes.
    fn update_distances(&mut self, i: usize, j: usize) {
        let d_ij = self.dist.get(i, j);
        self.row_sums[i] -= d_ij;

        for k in self.active.iter_ones() {
            if k == i || k == j {
                continue;
            }
            let d_ik = self.dist.get(i, k);
            let d_jk = self.dist.get(j, k);
            let d_new = 0.5 * (d_ik + d_jk - d_ij);
            self.row_sums[i] = self.row_sums[i] - d_ik + d_new;
            self.row_sums[k] = self.row_sums[k] - d_ik - d_jk + d_new;
            self.dist.set(i, k, d_new);
        }
        self.row_sums[j] = 0.0;
    }

    /// Connects the two remaining active nodes to form the root.
    ///
    /// The remaining distance `d(i,j)` is split equally between the two
    /// branches (`d/2` each), producing an unrooted-style midpoint root.
    fn final_join(mut self) -> Result<TreeNode, String> {
        let mut it = self.active.iter_ones();
        let i = it.next().expect("at least one active node");
        let j = it.next().expect("exactly two active nodes");

        let d_ij = self.dist.get(i, j);

        let mut left = self.nodes[i].take().expect("node i exists");
        let mut right = self.nodes[j].take().expect("node j exists");

        left.len = Some(d_ij / 2.0);
        right.len = Some(d_ij / 2.0);

        Ok(TreeNode::internal(
            self.next_internal,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            None,
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::distance_matrix::DistMat;

    #[test]
    fn test_nj_simple() {
        let names = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let dist = DistMat {
            names,
            data: vec![5.0, 9.0, 10.0, 9.0, 10.0, 8.0],
        };
        let tree = dist.neighbor_joining().unwrap();

        if let Some([left, right]) = &tree.children {
            assert_eq!(right.identifier, 3);
            assert_eq!(right.len.unwrap(), 2.0);

            if let Some([c, ab]) = &left.children {
                assert_eq!(c.identifier, 2);
                assert_eq!(c.len.unwrap(), 4.0);
                assert_eq!(ab.len.unwrap(), 3.0);

                if let Some([b, a]) = &ab.children {
                    assert_eq!(b.identifier, 1);
                    assert_eq!(a.identifier, 0);
                    assert_eq!(b.len.unwrap(), 3.0);
                    assert_eq!(a.len.unwrap(), 2.0);
                } else {
                    panic!("Left internal child should contain the A/B cherry");
                }
            } else {
                panic!("Left child of root should have two children");
            }
        } else {
            panic!("Root should have two children");
        }
    }
}
