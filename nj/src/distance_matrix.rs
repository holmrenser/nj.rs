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
use crate::models::{ModelCalculation, RateHet};
use crate::prepared::Prepared;
use crate::nj::NJState;
use crate::tree::TreeNode;
use serde::{Deserialize, Serialize};

/// Rows per block in [`DistMat::from_prepared`].
const ROW_BLOCK: usize = 16;

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
        let n_entries = n
            .checked_mul(n.saturating_sub(1))
            .expect("distance matrix too large: n * (n-1) overflows usize")
            / 2;
        Self {
            data: vec![0.0; n_entries],
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

    /// Replaces every non-finite entry (`inf` from a saturated pair, or `NaN`) by
    /// `max_distance` and returns the affected `(i, j)` pairs with `i > j`.
    ///
    /// Neighbor-joining cannot use infinite distances (infinite row sums turn into
    /// `NaN` in the update step), so callers apply this before
    /// [`neighbor_joining`](DistMat::neighbor_joining). A capped value stands for
    /// "at least this far". Finite entries, however large, are left untouched.
    pub fn cap_non_finite(&mut self, max_distance: f64) -> Vec<(usize, usize)> {
        let mut capped = Vec::new();
        for i in 1..self.dim() {
            for j in 0..i {
                let k = i * (i - 1) / 2 + j;
                if !self.data[k].is_finite() {
                    self.data[k] = max_distance;
                    capped.push((i, j));
                }
            }
        }
        capped
    }

    /// Computes a pairwise distance matrix from an MSA using substitution model `M`.
    ///
    /// Iterates over all `n*(n-1)/2` sequence pairs and fills the lower triangle
    /// using [`pairwise_distance`]. When the `parallel` feature is enabled, pairs
    /// are computed in parallel on the Rayon thread pool. Prefer [`MSA::into_dist`]
    /// over calling this directly.
    ///
    /// The `Send + Sync` bounds are required for the `parallel` feature and are
    /// trivially satisfied by all built-in alphabet and model types.
    pub fn from_msa<M, A>(msa: &MSA<A>, rates: RateHet) -> DistMat
    where
        M: ModelCalculation<A> + Send + Sync,
        A: AlphabetEncoding + Sync,
        A::Symbol: Send + Sync,
    {
        let prepared = Prepared::new(msa);
        Self::from_prepared::<M, A>(&prepared, msa.identifiers.clone(), rates)
    }

    /// Computes the matrix from a [`Prepared`] alignment.
    ///
    /// Rows are processed in blocks of [`ROW_BLOCK`]; each block owns a contiguous slice of
    /// the packed triangle and is filled column by column, so one sequence `j` is reused
    /// by every row of the block and the results are written in place (no per-pair
    /// intermediate vectors). With the `parallel` feature the blocks run on the Rayon
    /// thread pool.
    pub fn from_prepared<M, A>(data: &Prepared<'_, A>, names: Vec<String>, rates: RateHet) -> DistMat
    where
        M: ModelCalculation<A> + Send + Sync,
        A: AlphabetEncoding + Sync,
        A::Symbol: Send + Sync,
    {
        let mut dist = DistMat::empty_with_names(names);
        let n = dist.dim();
        debug_assert_eq!(n, data.len());

        // Entries stored before row `r` (row `r` has `r` entries).
        let before = |r: usize| if r == 0 { 0 } else { r * (r - 1) / 2 };

        let mut blocks: Vec<(usize, usize, &mut [f64])> = Vec::with_capacity(n / ROW_BLOCK + 1);
        let mut rest: &mut [f64] = &mut dist.data;
        let mut r0 = 0;
        while r0 < n {
            let r1 = (r0 + ROW_BLOCK).min(n);
            let (head, tail) = rest.split_at_mut(before(r1) - before(r0));
            blocks.push((r0, r1, head));
            rest = tail;
            r0 = r1;
        }

        let fill = |(r0, r1, out): (usize, usize, &mut [f64])| {
            let base = before(r0);
            for j in 0..r1 {
                for i in r0.max(j + 1)..r1 {
                    out[before(i) - base + j] = M::distance_prepared(data, i, j, rates);
                }
            }
        };

        #[cfg(feature = "parallel")]
        {
            use rayon::iter::{IntoParallelIterator, ParallelIterator};
            blocks.into_par_iter().for_each(fill);
        }
        #[cfg(not(feature = "parallel"))]
        blocks.into_iter().for_each(fill);

        dist
    }

    /// Converts to a [`DistanceResult`] by expanding to a full symmetric n×n matrix.
    pub fn into_result(self) -> DistanceResult {
        let n = self.dim();
        let matrix = (0..n)
            .map(|i| (0..n).map(|j| self.get(i, j)).collect())
            .collect();
        DistanceResult { names: self.names, matrix }
    }

    /// Borrows `self` to produce a [`DistanceResult`] without consuming the matrix.
    ///
    /// Useful when the distance matrix is needed both as output and for running NJ.
    pub fn to_result(&self) -> DistanceResult {
        let n = self.dim();
        let matrix = (0..n)
            .map(|i| (0..n).map(|j| self.get(i, j)).collect())
            .collect();
        DistanceResult { names: self.names.clone(), matrix }
    }

    /// Returns the mean of all `n*(n-1)/2` unique pairwise distances.
    ///
    /// Returns `0.0` for matrices with fewer than 2 taxa (no pairs exist).
    pub fn average(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().sum::<f64>() / self.data.len() as f64
    }

    /// Runs the Neighbor-Joining algorithm and returns the inferred tree.
    ///
    /// Consumes `self` because the NJ state machine mutates the matrix in
    /// place during the join iterations. Returns `Err` for an empty matrix,
    /// if any distance is not finite (call [`cap_non_finite`](DistMat::cap_non_finite)
    /// first), or if an internal consistency check fails.
    pub fn neighbor_joining(mut self: DistMat) -> Result<TreeNode, String> {
        for i in 1..self.dim() {
            for j in 0..i {
                if !self.get(i, j).is_finite() {
                    return Err(format!(
                        "distance between '{}' and '{}' is not finite; cap it with \
                         DistMat::cap_non_finite before running neighbor-joining",
                        self.names[i], self.names[j]
                    ));
                }
            }
        }
        NJState::new(&mut self).run()
    }
}

/// Serializable pairwise distance matrix returned by [`crate::distance_matrix`].
///
/// Expands the internal lower-triangular storage to a full symmetric n×n matrix
/// for ease of consumption in JSON, Python, and TypeScript.
#[derive(Clone, Debug, PartialEq, ts_rs::TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct DistanceResult {
    /// Sequence/taxa names, one per row/column.
    pub names: Vec<String>,
    /// Full symmetric n×n distance matrix. `matrix[i][j]` is the distance
    /// between taxa `i` and `j`; diagonal entries are `0.0`.
    pub matrix: Vec<Vec<f64>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::DNA;
    use crate::models::PDiff;
    use crate::msa::MSA;
    use crate::tree::NameOrSupport;

    fn named_matrix(n: usize) -> DistMat {
        let names = (0..n).map(|i| format!("t{i}")).collect();
        let mut d = DistMat::empty_with_names(names);
        for i in 1..n {
            for j in 0..i {
                d.set(i, j, 0.1 * (i + j) as f64 + 0.1);
            }
        }
        d
    }

    #[test]
    fn test_cap_non_finite_replaces_inf_and_nan_only() {
        let mut d = named_matrix(4);
        d.set(2, 0, f64::INFINITY);
        d.set(3, 1, f64::NAN);
        d.set(3, 2, 123.0); // large but finite: untouched
        let before = d.get(1, 0);
        let capped = d.cap_non_finite(10.0);
        assert_eq!(capped, vec![(2, 0), (3, 1)]);
        assert_eq!(d.get(2, 0), 10.0);
        assert_eq!(d.get(0, 2), 10.0); // symmetric access
        assert_eq!(d.get(3, 1), 10.0);
        assert_eq!(d.get(3, 2), 123.0);
        assert_eq!(d.get(1, 0), before);
    }

    #[test]
    fn test_cap_non_finite_is_noop_on_finite_matrix() {
        let mut d = named_matrix(5);
        let original = d.data.clone();
        assert!(d.cap_non_finite(10.0).is_empty());
        assert_eq!(d.data, original);
    }

    #[test]
    fn test_neighbor_joining_rejects_non_finite_distance() {
        let mut d = named_matrix(4);
        d.set(3, 1, f64::INFINITY);
        let err = d.neighbor_joining().unwrap_err();
        assert!(err.contains("'t3'") && err.contains("'t1'"), "{err}");
    }

    #[test]
    fn test_neighbor_joining_after_cap_has_finite_newick() {
        let mut d = named_matrix(5);
        d.set(4, 0, f64::INFINITY);
        d.cap_non_finite(10.0);
        let newick = d.neighbor_joining().unwrap().to_newick();
        assert!(!newick.contains("NaN") && !newick.contains("inf"), "{newick}");
    }

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
        let mat = msa.into_dist::<PDiff>(RateHet::NONE);
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
    fn test_dist_from_msa_gap_positions_excluded_from_denominator() {
        // AT- vs ACG: pos 2 is gap in first sequence → not comparable.
        // n_comparable=2, diffs=1 (T vs C at pos 1), so 1/2.
        let seqs: Vec<String> = vec!["AT-".into(), "ACG".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(seqs).unwrap();
        let mat = msa.into_dist::<PDiff>(RateHet::NONE);

        assert!((mat.get(0, 1) - (1.0 / 2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_dist_from_msa_no_overlap() {
        // no non-gap mismatches are counted, so the distance remains 0.0
        let seqs = vec![("X".into(), "A--".into()), ("Y".into(), "--A".into())];
        let msa = MSA::<DNA>::from_iter(seqs);
        let mat = msa.into_dist::<PDiff>(RateHet::NONE);
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
                None => (),
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

    // --- the prepared, tiled path must reproduce a plain per-pair loop bit for bit ---

    use crate::alphabet::{DNA as DnaAlphabet, Protein as ProteinAlphabet};
    use crate::models::{
        JukesCantor, Kimura2P, KimuraProtein, Poisson, TajimaNei, Tamura, pairwise_distance_with,
    };

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
    }

    fn random_msa<A: AlphabetEncoding>(rng: &mut Rng, alphabet: &[u8], n: usize, len: usize) -> MSA<A> {
        let pairs = (0..n).map(|i| {
            let seq: String = (0..len).map(|_| alphabet[(rng.next() % alphabet.len() as u64) as usize] as char).collect();
            (format!("s{i}"), seq)
        });
        MSA::<A>::try_from_pairs(pairs).unwrap()
    }

    fn assert_matches_reference<M, A>(msa: &MSA<A>, rates: RateHet, what: &str)
    where
        M: ModelCalculation<A> + Send + Sync,
        A: AlphabetEncoding + Sync,
        A::Symbol: Send + Sync,
    {
        let got = msa.into_dist::<M>(rates);
        for i in 0..msa.n_sequences {
            for j in 0..i {
                let want = pairwise_distance_with::<M, A>(&msa.sequences[i], &msa.sequences[j], rates);
                assert_eq!(got.get(i, j).to_bits(), want.to_bits(), "{what}: pair ({i},{j})");
            }
        }
    }

    #[test]
    fn test_prepared_matrix_matches_reference_dna() {
        let mut rng = Rng(0xABCD_EF01_2345_6789);
        // Sizes around the row-block boundary; lengths around the 64-column word boundary.
        for n in [1usize, 2, 3, 15, 16, 17, 31, 33, 50] {
            for len in [1usize, 63, 64, 65, 200] {
                let msa = random_msa::<DnaAlphabet>(&mut rng, b"ACGTACGTACGTACGTNR-", n, len);
                let what = format!("n={n} len={len}");
                for rates in [RateHet::NONE, RateHet { gamma_shape: Some(0.5), p_invar: 0.2 }] {
                    assert_matches_reference::<PDiff, _>(&msa, rates, &what);
                    assert_matches_reference::<JukesCantor, _>(&msa, rates, &what);
                    assert_matches_reference::<Kimura2P, _>(&msa, rates, &what);
                    assert_matches_reference::<TajimaNei, _>(&msa, rates, &what);
                    assert_matches_reference::<Tamura, _>(&msa, rates, &what);
                }
            }
        }
    }

    #[test]
    fn test_prepared_matrix_matches_reference_protein() {
        let mut rng = Rng(0x0F0F_1234_5678_9ABC);
        for n in [1usize, 2, 16, 17, 40] {
            for len in [1usize, 65, 150] {
                let msa = random_msa::<ProteinAlphabet>(&mut rng, b"ARNDCQEGHILKMFPSTWYVXX-", n, len);
                let what = format!("n={n} len={len}");
                assert_matches_reference::<PDiff, _>(&msa, RateHet::NONE, &what);
                assert_matches_reference::<Poisson, _>(&msa, RateHet::NONE, &what);
                assert_matches_reference::<KimuraProtein, _>(&msa, RateHet::NONE, &what);
            }
        }
    }
}
