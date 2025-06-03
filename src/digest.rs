use crate::b_field_element::BFieldElement;

/// The result of hashing a sequence of elements, for [Tip5].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Digest(pub [BFieldElement; Digest::LEN]);

impl Digest {
    /// The number of [elements](BFieldElement) in a digest.
    pub const LEN: usize = 5;

    /// Creates a new digest from an array of elements.
    pub const fn new(elements: [BFieldElement; Self::LEN]) -> Self {
        Self(elements)
    }

    pub const fn values(self) -> [BFieldElement; Self::LEN] {
        self.0
    }
}
