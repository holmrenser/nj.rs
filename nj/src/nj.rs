//! Core Neighbor-Joining algorithm implementation.
//!
//! The algorithm runs in three phases:
//! 1. **Initialisation** ([`NJState::new`]) — build leaf nodes and precompute the row
//!    sums of the distance matrix.
//! 2. **Iteration** ([`NJState::run`]) — repeat `n-2` times: select the pair `(i, j)`
//!    minimising the Q-criterion, join them into a new internal node, update the
//!    distance matrix, and retire one slot.
//! 3. **Final join** ([`NJState::final_join`]) — connect the two remaining nodes to
//!    form the root, splitting the remaining distance equally.
//!
//! # Layout
//!
//! The matrix is updated in place and kept **compact**: active nodes occupy slots
//! `0..m` of the packed lower triangle, so row `r` is one contiguous slice
//! `data[r(r-1)/2 .. r(r-1)/2 + r]` and the Q search is a plain reduction over slices
//! with no per-element active-set test or index arithmetic. When two nodes join, the
//! retired slot is refilled by moving the last slot into it (swap-remove, O(m)).
//!
//! Slot order therefore differs from the original taxon order, so every place where
//! the order of operations can change a result is pinned to *original* indices:
//!
//! * `Q_ij = (m−2)·d_ij − (S_i + S_j)` is symmetric in `i` and `j`, so its value does
//!   not depend on which slot holds which node;
//! * ties are resolved by the smallest original `(i, j)` pair (`i > j`), as in a scan
//!   over original indices;
//! * the incremental row-sum updates visit nodes in ascending original index, so the
//!   floating-point rounding is the same as for a scan over original indices.
//!
//! The straightforward implementation is kept under `cfg(test)` in
//! `nj_reference.rs`; the tests require bit-identical trees.

use crate::{distance_matrix::DistMat, tree::TreeNode};

/// Index of the pair `(a, b)`, `a != b`, in the packed lower triangle.
#[inline(always)]
fn tri(a: usize, b: usize) -> usize {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi * (hi - 1) / 2 + lo
}

/// Computes the branch lengths from the new internal node to taxa `i` and `j`.
///
/// Uses the standard NJ formula:
/// ```text
/// l_i = 0.5 * d(i,j) + (S_i - S_j) / (2 * (m - 2))
/// l_j = d(i,j) - l_i
/// ```
/// where `S_i` and `S_j` are the row sums and `m` is the number of currently
/// active taxa. Both lengths are clamped to `0.0` to prevent negative branch
/// lengths from rounding or non-ultrametric inputs.
fn compute_branch_lengths(d_ij: f64, s_i: f64, s_j: f64, active_count: usize) -> (f64, f64) {
    let n = active_count as f64;

    let li = (0.5 * d_ij + (s_i - s_j) / (2.0 * (n - 2.0))).max(0.0);
    let lj = (d_ij - li).max(0.0);

    (li, lj)
}

/// Minimum over `t` of `scale·d[t] − (s_i + s[t])`: the Q criterion for one row.
///
/// Every element is computed with the same two separately rounded operations as the
/// scalar reference (no fused multiply-add), and `min` is exact, so the result does
/// not depend on how the loop is vectorised.
#[inline]
fn row_min_q(d: &[f64], s: &[f64], scale: f64, s_i: f64) -> f64 {
    debug_assert_eq!(d.len(), s.len());
    let mut m = [f64::INFINITY; 4];
    let (cd, cs) = (d.chunks_exact(4), s.chunks_exact(4));
    let (rd, rs) = (cd.remainder(), cs.remainder());
    for (x, y) in cd.zip(cs) {
        for l in 0..4 {
            m[l] = m[l].min(scale * x[l] - (s_i + y[l]));
        }
    }
    let mut r = m[0].min(m[1]).min(m[2].min(m[3]));
    for (x, y) in rd.iter().zip(rs) {
        r = r.min(scale * x - (s_i + y));
    }
    r
}

/// Mutable state for a single NJ run.
///
/// Owns the partial tree nodes and works on the distance matrix in place. Slots
/// `0..m` hold the active nodes; `orig` maps a slot to the index its node would have
/// in a scan over original indices (the larger index of a joined pair survives).
pub struct NJState<'a> {
    dist: &'a mut DistMat,
    /// Number of active nodes (slots `0..m`).
    m: usize,
    /// Slot → original index.
    orig: Vec<usize>,
    /// Original index → slot (valid for active nodes).
    slot_of: Vec<usize>,
    /// Active original indices in ascending order.
    active_sorted: Vec<usize>,
    /// Row sums, by slot.
    row_sums: Vec<f64>,
    /// Partial trees, by slot.
    nodes: Vec<Option<TreeNode>>,
    next_internal: usize,
}

impl<'a> NJState<'a> {
    /// Initialises the NJ state from a distance matrix.
    ///
    /// Creates one leaf [`TreeNode`] per taxon, precomputes the full row sums
    /// `S_i = Σ_j d(i,j)`, and places taxon `i` in slot `i`.
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
            m: n,
            orig: (0..n).collect(),
            slot_of: (0..n).collect(),
            active_sorted: (0..n).collect(),
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
            let (si, sj) = self.select_min_q_pair().ok_or("No pair found")?;
            let d_ij = self.dist.data[tri(si, sj)];
            let (li, lj) = compute_branch_lengths(
                d_ij,
                self.row_sums[si],
                self.row_sums[sj],
                self.m,
            );
            self.join_nodes(si, sj, li, lj);
            self.update_distances(si, sj, d_ij);
            self.retire(sj);
        }

        if self.m != 2 {
            return Err(format!("Expected 2 active nodes but found {}", self.m));
        }

        self.final_join()
    }

    /// Selects the active pair with the smallest Q value and returns its slots
    /// `(si, sj)`, where `si` holds the node with the larger original index.
    ///
    /// `Q(i,j) = (m-2)·d(i,j) - (S_i + S_j)`. Each row's minimum is found with one
    /// pass over its slice; only rows whose minimum equals the global one are
    /// rescanned, and among equal Q values the smallest original `(i, j)` wins.
    fn select_min_q_pair(&self) -> Option<(usize, usize)> {
        if self.m < 2 {
            return None;
        }
        let scale = (self.m - 2) as f64;
        let (data, s) = (&self.dist.data, &self.row_sums);

        let mut best = f64::INFINITY;
        let mut row_min = vec![f64::INFINITY; self.m];
        for r in 1..self.m {
            let start = r * (r - 1) / 2;
            let q = row_min_q(&data[start..start + r], &s[..r], scale, s[r]);
            row_min[r] = q;
            if q < best {
                best = q;
            }
        }

        let mut pick: Option<(usize, usize, usize, usize)> = None; // (hi_orig, lo_orig, slot_hi, slot_lo)
        for r in 1..self.m {
            if row_min[r] != best {
                continue;
            }
            let start = r * (r - 1) / 2;
            for t in 0..r {
                if scale * data[start + t] - (s[r] + s[t]) == best {
                    let (a, b) = (self.orig[r], self.orig[t]);
                    let cand = if a > b { (a, b, r, t) } else { (b, a, t, r) };
                    if pick.is_none_or(|p| (cand.0, cand.1) < (p.0, p.1)) {
                        pick = Some(cand);
                    }
                }
            }
        }
        pick.map(|(_, _, hi_slot, lo_slot)| (hi_slot, lo_slot))
    }

    /// Merges the nodes in slots `si` (larger original index) and `sj` into a new
    /// internal node stored in slot `si`.
    ///
    /// Assigns branch lengths `li` and `lj` to the child nodes. The node in slot
    /// `sj` is consumed here; the caller retires the slot.
    fn join_nodes(&mut self, si: usize, sj: usize, li: f64, lj: f64) {
        let mut left = self.nodes[si].take().expect("node i exists");
        left.len = Some(li);

        let mut right = self.nodes[sj].take().expect("node j exists");
        right.len = Some(lj);

        self.nodes[si] = Some(TreeNode::internal(
            self.next_internal,
            Some([Box::new(left), Box::new(right)]),
            Some(0.0),
            None,
        ));
        self.next_internal += 1;
    }

    /// Updates the distance matrix and row sums after joining the nodes in slots
    /// `si` and `sj`.
    ///
    /// For every remaining node `k` (excluding the joined pair), applies the classic
    /// NJ distance update `d(new, k) = 0.5 * (d(i,k) + d(j,k) - d(i,j))`. Row sums are
    /// maintained incrementally, visiting nodes in ascending **original** index so the
    /// floating-point rounding matches a scan over original indices.
    fn update_distances(&mut self, si: usize, sj: usize, d_ij: f64) {
        let (oi, oj) = (self.orig[si], self.orig[sj]);
        self.row_sums[si] -= d_ij;

        for &ok in &self.active_sorted {
            if ok == oi || ok == oj {
                continue;
            }
            let sk = self.slot_of[ok];
            let d_ik = self.dist.data[tri(si, sk)];
            let d_jk = self.dist.data[tri(sj, sk)];
            let d_new = 0.5 * (d_ik + d_jk - d_ij);
            self.row_sums[si] = self.row_sums[si] - d_ik + d_new;
            self.row_sums[sk] = self.row_sums[sk] - d_ik - d_jk + d_new;
            self.dist.data[tri(si, sk)] = d_new;
        }
    }

    /// Removes slot `sj` by moving the last active slot into it, keeping the active
    /// rows contiguous.
    fn retire(&mut self, sj: usize) {
        let oj = self.orig[sj];
        let pos = self
            .active_sorted
            .binary_search(&oj)
            .expect("retired node is active");
        self.active_sorted.remove(pos);

        let last = self.m - 1;
        if sj != last {
            for t in 0..last {
                if t != sj {
                    self.dist.data[tri(sj, t)] = self.dist.data[tri(last, t)];
                }
            }
            self.row_sums[sj] = self.row_sums[last];
            self.orig[sj] = self.orig[last];
            self.slot_of[self.orig[sj]] = sj;
            self.nodes[sj] = self.nodes[last].take();
        }
        self.m -= 1;
    }

    /// Connects the two remaining nodes to form the root.
    ///
    /// The remaining distance `d(i,j)` is split equally between the two
    /// branches (`d/2` each), producing an unrooted-style midpoint root. The node
    /// with the smaller original index is the left child.
    fn final_join(mut self) -> Result<TreeNode, String> {
        let (first, second) = if self.orig[0] < self.orig[1] { (0, 1) } else { (1, 0) };
        let d_ij = self.dist.data[tri(0, 1)];

        let mut left = self.nodes[first].take().expect("node i exists");
        let mut right = self.nodes[second].take().expect("node j exists");

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

    // --- the compacted implementation must reproduce the reference bit for bit ---

    use crate::nj_reference;
    use crate::tree::{NameOrSupport, TreeNode};

    fn same_node(a: &TreeNode, b: &TreeNode) -> bool {
        let same_label = match (&a.label, &b.label) {
            (None, None) => true,
            (Some(NameOrSupport::Name(x)), Some(NameOrSupport::Name(y))) => x == y,
            (Some(NameOrSupport::Support(x)), Some(NameOrSupport::Support(y))) => x == y,
            _ => false,
        };
        let same_len = match (a.len, b.len) {
            (None, None) => true,
            (Some(x), Some(y)) => x.to_bits() == y.to_bits(),
            _ => false,
        };
        let same_children = match (&a.children, &b.children) {
            (None, None) => true,
            (Some([al, ar]), Some([bl, br])) => same_node(al, bl) && same_node(ar, br),
            _ => false,
        };
        a.identifier == b.identifier && same_label && same_len && same_children
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn matrix(n: usize, mut value: impl FnMut(usize, usize) -> f64) -> DistMat {
        let names = (0..n).map(|i| format!("t{i}")).collect();
        let mut d = DistMat::empty_with_names(names);
        for i in 1..n {
            for j in 0..i {
                d.set(i, j, value(i, j));
            }
        }
        d
    }

    fn assert_matches_reference(d: DistMat, what: &str) {
        let expected = nj_reference::NJState::new(&mut d.clone()).run().unwrap();
        let got = d.neighbor_joining().unwrap();
        assert!(
            same_node(&expected, &got),
            "{what}: trees differ\n  reference {}\n  new       {}",
            expected.to_newick(),
            got.to_newick()
        );
    }

    #[test]
    fn test_matches_reference_on_random_matrices() {
        let mut r = Rng(0xC0FF_EE11);
        for n in 1..=45 {
            for rep in 0..6 {
                assert_matches_reference(
                    matrix(n, |_, _| 0.05 + r.unit()),
                    &format!("random n={n} rep={rep}"),
                );
            }
        }
        assert_matches_reference(matrix(150, |_, _| 0.05 + r.unit()), "random n=150");
    }

    #[test]
    fn test_matches_reference_with_many_ties() {
        // Distances drawn from a handful of integers: most Q values tie exactly.
        let mut r = Rng(0xBEEF_0001);
        for n in 3..=40 {
            for rep in 0..6 {
                assert_matches_reference(
                    matrix(n, |_, _| (1 + r.next() % 3) as f64),
                    &format!("integer ties n={n} rep={rep}"),
                );
            }
        }
    }

    #[test]
    fn test_matches_reference_with_duplicate_taxa() {
        // Taxa that are copies of each other (zero distance, identical rows).
        let mut r = Rng(0xD0D0_0002);
        for n in 4..=40 {
            let base: Vec<usize> = (0..n).map(|i| if r.next() % 3 == 0 && i > 0 { i - 1 } else { i }).collect();
            let pos: Vec<f64> = (0..n).map(|_| r.unit() * 4.0).collect();
            assert_matches_reference(
                matrix(n, |i, j| {
                    let (a, b) = (base[i], base[j]);
                    if a == b { 0.0 } else { (pos[a] - pos[b]).abs() + 0.25 }
                }),
                &format!("duplicates n={n}"),
            );
        }
    }

    #[test]
    fn test_matches_reference_on_degenerate_matrices() {
        assert_matches_reference(matrix(2, |_, _| 0.7), "two taxa");
        assert_matches_reference(matrix(3, |_, _| 1.0), "three equal");
        assert_matches_reference(matrix(12, |_, _| 1.0), "all equal");
        assert_matches_reference(matrix(12, |_, _| 0.0), "all zero");
    }
}
