use bitvec::prelude::{BitVec, Lsb0};

use crate::{dist::DistMat, tree::TreeNode};

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
pub struct NJState<'a> {
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

        // Initialize tree nodes as leaves with identifiers corresponding to their index.
        let nodes = dist
            .names
            .iter()
            .enumerate()
            .map(|(i, name)| Some(TreeNode::leaf(i, name.clone(), None)))
            .collect();

        // Precompute row sums for the initial distance matrix.
        let row_sums = (0..n)
            .map(|i| (0..n).map(|j| dist.get(i, j)).sum())
            .collect();

        // All nodes start as active (not yet joined).
        let active = BitVec::repeat(true, n);

        NJState {
            dist,
            active,
            row_sums,
            nodes,
            next_internal: n, // internal nodes get identifiers starting after the leaves
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
            // Select pair i,j minimizing Q-metric.
            let (i, j, d_ij) = self.select_min_q_pair().ok_or("No pair found")?;
            // Compute branch lengths for the new node to i and j.
            let (li, lj) =
                compute_branch_lengths(d_ij, &self.row_sums, i, j, self.active.count_ones());

            // New node goes in position i, j is removed
            self.join_nodes(i, j, li, lj);
            // Update distances and row sums after joining.
            // It is important to do this before marking j as inactive, since we need j's distances.
            self.update_distances(i, j);
            // Mark j as inactive.
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
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Greater))
            .map(|(i, j, _, d)| (i, j, d))
    }

    /// Join nodes i and j into a new internal node at index i.
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

    // Update distances after joining i (the new node) and removing j.
    fn update_distances(&mut self, i: usize, j: usize) {
        // Get the distance between i and j before they are joined, which is needed for the update.
        let d_ij = self.dist.get(i, j);

        self.row_sums[i] -= d_ij;

        // Update distances to the new node i and maintain row sums incrementally.
        for k in self.active.iter_ones() {
            // Skip j since it's being removed and skip i since it's the new node.
            if k == i || k == j {
                continue;
            }
            let d_ik = self.dist.get(i, k);
            let d_jk = self.dist.get(j, k);
            // Classic NJ update formula
            let d_new = 0.5 * (d_ik + d_jk - d_ij);

            // Update row sums for i and k. We can do this incrementally by subtracting the old distances and adding the new one.
            self.row_sums[i] = self.row_sums[i] - d_ik + d_new;
            self.row_sums[k] = self.row_sums[k] - d_ik - d_jk + d_new;
            // Update the distance matrix with the new distance to i.
            self.dist.set(i, k, d_new);
        }
        self.row_sums[j] = 0.0;
    }

    /// Final join of the last two remaining active nodes.
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
    use crate::dist::DistMat;

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
