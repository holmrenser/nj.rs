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

        let nodes = dist
            .names
            .iter()
            .enumerate()
            .map(|(i, name)| Some(TreeNode::leaf(i, name.clone(), None)))
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
        left.len = Some(li);

        let mut right = self.nodes[j].take().expect("node j exists");
        right.len = Some(lj);

        // let name = format!("Node{}", self.next_internal);
        self.next_internal += 1;

        self.nodes[i] = Some(TreeNode::internal(
            self.next_internal,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            None,
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

        left.len = Some(d_ij / 2.0);
        right.len = Some(d_ij / 2.0);

        // let name = format!("Node{}", self.next_internal);
        TreeNode::internal(
            self.next_internal + 1,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            None,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::DistMat;
    use crate::tree::NameOrSupport;

    #[test]
    fn test_nj_simple() {
        let names = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let data = vec![
            0.0, 5.0, 9.0, 9.0, //
            5.0, 0.0, 10.0, 10.0, //
            9.0, 10.0, 0.0, 8.0, //
            9.0, 10.0, 8.0, 0.0, //
        ];
        let dist = DistMat { names, data }; // 4x4 matrix   
        let tree = dist.neighbor_joining().unwrap();
        // The expected tree structure is ((A:2.0,B:3.0):5.0,(C:4.0,D:4.0):3.0);
        if let Some([left, right]) = &tree.children {
            // Check left subtree (A,B)
            if let Some([a, b]) = &left.children {
                assert_eq!(a.identifier, 1);
                assert_eq!(b.identifier, 0);
                assert_eq!(a.len.unwrap(), 2.0);
                assert_eq!(b.len.unwrap(), 3.0);
            } else {
                panic!("Left child of root should have two children");
            }
            // Check right subtree (C,D)
            if let Some([c, d]) = &right.children {
                assert_eq!(c.identifier, 2);
                assert_eq!(d.identifier, 3);
                assert_eq!(c.len.unwrap(), 4.0);
                assert_eq!(d.len.unwrap(), 4.0);
            } else {
                panic!("Right child of root should have two children");
            }
        } else {
            panic!("Root should have two children");
        }
    }
}
