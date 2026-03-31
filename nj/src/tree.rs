use std::fmt::{Display, Formatter, Result as FmtResult};

/// A node in a (bifurcating) phylogenetic tree.
/// Can be either a leaf node (no children) or an internal node (two children).
/// Each branch has an associated length.
/// Leaf nodes represent sequences, while internal nodes represent common ancestors.
/// Building a Tree means recursively combining nodes until a single root node remains.
#[derive(Clone, Debug)]
pub struct TreeNode {
    /// Unique identifier for the node.
    pub identifier: usize,
    /// Label for the node: either a name (for leaves) or bootstrap support (for internal nodes).
    pub label: Option<NameOrSupport>,
    /// Child nodes (None for leaf nodes).
    pub children: Option<[Box<TreeNode>; 2]>,
    /// Branch length leading to this node.
    pub len: Option<f64>,
}

/// Label for a tree node: either a name (String) or bootstrap support (usize).
#[derive(Clone, Debug)]
pub enum NameOrSupport {
    /// Name label for leaf nodes.
    Name(String),
    /// Bootstrap support for internal nodes.
    Support(usize),
}

impl Display for TreeNode {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", self.to_newick())
    }
}

impl TreeNode {
    /// Creates a leaf node with the given name.
    pub fn leaf(identifier: usize, name: String, len: Option<f64>) -> Self {
        Self {
            identifier,
            label: Some(NameOrSupport::Name(name)),
            children: None,
            len,
        }
    }
    /// Creates an internal node with the given name, left and right children, and branch lengths.
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

    /// Recursively converts the tree to Newick format.
    /// Returns the Newick string representation of the tree (without the trailing semicolon).
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

    /// Converts the tree to Newick format (with trailing semicolon).
    pub fn to_newick(&self) -> String {
        let mut newick = self.to_newick_recursion();
        newick.push(';');
        newick
    }
}
