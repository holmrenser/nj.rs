/// A node in a (bifurcating) phylogenetic tree.
/// Can be either a leaf node (no children) or an internal node (two children).
/// Each branch has an associated length.
/// Leaf nodes represent sequences, while internal nodes represent common ancestors.
/// Building a Tree means recursively combining nodes until a single root node remains.
#[derive(Clone, Debug)]
pub struct TreeNode {
    pub name: String,
    pub children: Option<[Box<TreeNode>; 2]>,
    pub len: f64,
}

impl TreeNode {
    /// Creates a leaf node with the given name.
    pub fn leaf(name: String, len: Option<f64>) -> Self {
        Self {
            name,
            children: None,
            len: match len {
                Some(l) => l,
                None => 0.0,
            },
        }
    }
    /// Creates an internal node with the given name, left and right children, and branch lengths.
    pub fn internal(name: String, children: Option<[Box<TreeNode>; 2]>, len: f64) -> Self {
        Self {
            name,
            children: Some(children.unwrap()),
            len,
        }
    }

    /// Recursively converts the tree to Newick format.
    /// If hide_internal is true, internal node names are omitted.
    /// Returns the Newick string representation of the tree.
    fn to_newick_recursion(&self, hide_internal: bool) -> String {
        match &self.children {
            Some([left, right]) => {
                let left_str = left.to_newick(hide_internal);
                let right_str = right.to_newick(hide_internal);
                let name_str = if hide_internal {
                    "".to_string()
                } else {
                    self.name.clone()
                };
                format!(
                    "({}:{:.3},{}:{:.3}){}",
                    left_str, left.len, right_str, right.len, name_str
                )
            }
            None => self.name.clone(),
        }
    }

    /// Converts the tree to Newick format (with trailing semicolon).
    pub fn to_newick(&self, hide_internal: bool) -> String {
        let mut newick = self.to_newick_recursion(hide_internal);
        newick.push(';');
        newick
    }
}
