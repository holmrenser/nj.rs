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
mod dist;
mod msa;
mod tree;
use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::HashMap;

pub use crate::config::{FastaSequence, MSA, NJConfig};
use crate::dist::DistMat;
use crate::tree::TreeNode;

/// Collect leaves into a BitVec according to index map.
fn bitset_of(node: &TreeNode, idx: &HashMap<String, usize>, out: &mut BitVec<u8, Lsb0>) {
    match &node.children {
        None => {
            let i = idx[&node.name];
            out.set(i, true);
        }
        Some([l, r]) => {
            bitset_of(l, idx, out);
            bitset_of(r, idx, out);
        }
    }
}

pub fn count_clades_bitvec_noset(
    tree: &TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counter: &mut HashMap<Vec<u8>, usize>,
) {
    if let Some([l, r]) = &tree.children {
        // compute clade bitvec
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(tree, idx, &mut bv);

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            // unique by structure; no HashSet needed
            counter
                .entry(bv.as_raw_slice().to_vec())
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }

        // recursion
        count_clades_bitvec_noset(l, idx, n_taxa, counter);
        count_clades_bitvec_noset(r, idx, n_taxa, counter);
    }
}

pub fn bootstrap_clade_counts(
    trees: &[TreeNode],
    idx: &HashMap<String, usize>,
    n_taxa: usize,
) -> HashMap<Vec<u8>, usize> {
    let mut counter = HashMap::new();

    for t in trees {
        count_clades_bitvec_noset(t, idx, n_taxa, &mut counter);
    }

    counter
}

/// Main entry point for Neighbor-Joining tree construction.
/// Takes an NJConfig with MSA sequences, internal node name hiding option,
/// and number of bootstrap samples.
pub fn nj(conf: NJConfig) -> Result<String, String> {
    let msa = MSA::from_iter(conf.msa.into_iter());
    let clade_counts = if conf.n_bootstrap_samples > 0 {
        let mut idx_map = HashMap::new();
        for (i, fs) in msa.sequences.iter().enumerate() {
            idx_map.insert(fs.identifier.clone(), i);
        }
        let bootstap_trees: Vec<TreeNode> = (0..conf.n_bootstrap_samples)
            .map(|i| {
                println!(
                    "Generating bootstrap tree {}/{}",
                    i + 1,
                    conf.n_bootstrap_samples
                );
                return msa.bootstrap().into_dist().neighbor_joining().unwrap();
            })
            .collect();
        Some(bootstrap_clade_counts(
            &bootstap_trees,
            &idx_map,
            msa.n_sequences,
        ))
    } else {
        None
    };

    let main_tree = msa.into_dist().neighbor_joining()?;
    let newick = match clade_counts {
        Some(counts) => {
            // annotate internal nodes with bootstrap support
            fn annotate_bootstrap(
                node: &mut TreeNode,
                idx: &HashMap<String, usize>,
                n_taxa: usize,
                counts: &HashMap<Vec<u8>, usize>,
            ) {
                if node.children.is_some() {
                    // compute clade bitvec
                    let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
                    bitset_of(node, idx, &mut bv);

                    let n = bv.count_ones();
                    if n > 1 && n < n_taxa {
                        if let Some(c) = counts.get(&bv.as_raw_slice().to_vec()) {
                            node.name = format!("{}:{}", node.name, c);
                        }
                    }

                    // recursion
                    if let Some([l, r]) = &mut node.children {
                        annotate_bootstrap(l, idx, n_taxa, counts);
                        annotate_bootstrap(r, idx, n_taxa, counts);
                    }
                }
            }

            let mut idx_map = HashMap::new();
            for (i, fs) in msa.sequences.iter().enumerate() {
                idx_map.insert(fs.identifier.clone(), i);
            }
            let mut main_tree_mut = main_tree;
            annotate_bootstrap(&mut main_tree_mut, &idx_map, msa.n_sequences, &counts);
            main_tree_mut.to_newick(conf.show_internal)
        }
        None => main_tree.to_newick(conf.show_internal),
    };
    //let newick = main_tree.to_newick(conf.hide_internal);
    Ok(newick)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nj_wrapper_adds_semicolon() {
        let msa = MSA::from_unnamed_sequences(vec!["A".into(), "A".into()]);
        let conf = NJConfig {
            msa: msa.sequences,
            show_internal: false,
            n_bootstrap_samples: 1,
        };
        let out = nj(conf).unwrap();
        assert!(out.ends_with(';'));
    }

    #[test]
    fn test_nj_deterministic_order() {
        let msa = MSA::from_unnamed_sequences(vec!["ACG".into(), "ATG".into(), "AGG".into()]);
        let conf = NJConfig {
            msa: msa.sequences,
            show_internal: false,
            n_bootstrap_samples: 1,
        };

        let t1 = nj(conf.clone()).unwrap();
        let t2 = nj(conf).unwrap();
        assert_eq!(t1, t2);
    }
}
