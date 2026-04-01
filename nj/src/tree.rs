//! Binary tree representation for Neighbor-Joining output.
//!
//! Trees are built bottom-up by the NJ algorithm: the algorithm starts with
//! leaf nodes and repeatedly joins pairs into new internal nodes. The final
//! result is a single [`TreeNode`] at the root.

use std::fmt::{Display, Formatter, Result as FmtResult};

/// A node in a strictly bifurcating (binary) phylogenetic tree.
///
/// Every node is either a **leaf** (`children` is `None`, `label` is a
/// [`NameOrSupport::Name`]) or an **internal node** (`children` holds exactly
/// two boxed children, `label` is optionally a [`NameOrSupport::Support`]).
///
/// Branch lengths are stored on the *child* node rather than on the edge —
/// `len` is the length of the branch connecting this node to its parent.
/// The root has `len = Some(0.0)` by convention.
#[derive(Clone, Debug)]
pub struct TreeNode {
    /// Numeric identifier assigned during NJ construction. Leaves receive the
    /// index of the original sequence (0-based); internal nodes receive
    /// sequential identifiers starting after the last leaf index.
    pub identifier: usize,
    /// Node label: a sequence name for leaves, or a bootstrap count for
    /// internal nodes (assigned after bootstrap sampling). `None` for internal
    /// nodes that have not been assigned a support value.
    pub label: Option<NameOrSupport>,
    /// Child nodes; `None` for leaf nodes, `Some([left, right])` for internal.
    pub children: Option<[Box<TreeNode>; 2]>,
    /// Length of the branch from this node to its parent. `None` before the
    /// branch length is computed; `Some(0.0)` for the root.
    pub len: Option<f64>,
}

/// Label carried by a [`TreeNode`]: either a sequence name or a bootstrap count.
///
/// Leaves always carry a [`Name`](NameOrSupport::Name). Internal nodes start
/// with no label and may have a [`Support`](NameOrSupport::Support) count
/// attached after bootstrap analysis via [`crate::add_bootstrap_to_tree`].
#[derive(Clone, Debug)]
pub enum NameOrSupport {
    /// Sequence identifier (leaf nodes only).
    Name(String),
    /// Number of bootstrap replicates in which this clade was recovered.
    Support(usize),
}

impl Display for TreeNode {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", self.to_newick())
    }
}

impl TreeNode {
    /// Creates a leaf node with the given sequence name and optional branch length.
    pub fn leaf(identifier: usize, name: String, len: Option<f64>) -> Self {
        Self {
            identifier,
            label: Some(NameOrSupport::Name(name)),
            children: None,
            len,
        }
    }
    /// Creates an internal node with children and an optional bootstrap support value.
    ///
    /// # Panics
    /// Panics if `children` is `None` (the NJ algorithm always provides children
    /// for internal nodes, so this is a programming error).
    pub fn internal(
        identifier: usize,
        children: Option<[Box<TreeNode>; 2]>,
        len: Option<f64>,
        support: Option<usize>,
    ) -> Self {
        Self {
            identifier,
            children: Some(children.unwrap()),
            len,
            label: match support {
                Some(s) => Some(NameOrSupport::Support(s)),
                None => None,
            },
        }
    }

    /// Recursively serialises the subtree rooted at `self` to Newick notation.
    ///
    /// Branch lengths are formatted with three decimal places (`:.3`). Internal
    /// node labels (bootstrap support values) are appended after the closing
    /// parenthesis. The trailing `;` is **not** included — call [`to_newick`]
    /// for the complete, semicolon-terminated string.
    ///
    /// [`to_newick`]: TreeNode::to_newick
    fn to_newick_recursion(&self) -> String {
        match &self.children {
            Some([left, right]) => {
                let left_part = match left.len {
                    Some(len) => {
                        let left_str = left.to_newick_recursion();
                        format!("{}:{:.3}", left_str, len)
                    }
                    None => left.to_newick_recursion(),
                };
                let right_part = match right.len {
                    Some(len) => {
                        let right_str = right.to_newick_recursion();
                        format!("{}:{:.3}", right_str, len)
                    }
                    None => right.to_newick_recursion(),
                };

                let label_str = match self.label {
                    Some(NameOrSupport::Support(s)) => format!("{}", s),
                    Some(NameOrSupport::Name(ref n)) => n.clone(),
                    None => "".to_string(),
                };

                format!("({},{}){}", left_part, right_part, label_str)
            }
            None => match self.label {
                Some(NameOrSupport::Support(s)) => format!("{}", s),
                Some(NameOrSupport::Name(ref n)) => n.clone(),
                None => "".to_string(),
            },
        }
    }

    /// Returns the full Newick representation of this tree, including the
    /// trailing semicolon required by the Newick standard.
    ///
    /// Branch lengths are rounded to three decimal places. Bootstrap support
    /// values (if present) appear as integer labels on internal nodes.
    pub fn to_newick(&self) -> String {
        let mut newick = self.to_newick_recursion();
        newick.push(';');
        newick
    }
}
