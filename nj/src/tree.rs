/// NJ tree node
#[derive(Clone, Debug)]
pub struct NJNode {
    pub name: String,
    pub left: Option<Box<NJNode>>,
    pub right: Option<Box<NJNode>>,
    pub left_len: f64,
    pub right_len: f64,
}

impl NJNode {
    pub fn leaf(name: String) -> Self {
        Self {
            name,
            left: None,
            right: None,
            left_len: 0.0,
            right_len: 0.0,
        }
    }

    pub fn internal(
        name: String,
        left: NJNode,
        right: NJNode,
        left_len: f64,
        right_len: f64,
    ) -> Self {
        Self {
            name,
            left: Some(Box::new(left)),
            right: Some(Box::new(right)),
            left_len,
            right_len,
        }
    }

    /// Exports the tree in Newick format.
    ///
    /// Branch lengths are formatted to three decimal places.
    /// If `hide_internal` is true, internal node names are omitted from the output.
    ///
    /// # Panics
    /// This function will panic if called on a node with only one child (invalid NJNode).
    pub fn to_newick(&self, hide_internal: bool) -> String {
        match (&self.left, &self.right) {
            (Some(left), Some(right)) => {
                let left_str = left.to_newick(hide_internal);
                let right_str = right.to_newick(hide_internal);
                let name_str = if hide_internal {
                    "".to_string()
                } else {
                    self.name.clone()
                };
                format!(
                    "({}:{:.3},{}:{:.3}){}",
                    left_str, self.left_len, right_str, self.right_len, name_str
                )
            }
            (None, None) => self.name.clone(),
            _ => panic!(
                "NJNode must have either zero or two children; found a node with only one child, which is invalid."
            ),
        }
    }
}
