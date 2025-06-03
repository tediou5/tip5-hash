use std::iter;

use crate::b_field_element::BFieldElement;
use itertools::Itertools;
use num_traits::ConstOne;
use num_traits::ConstZero;

pub const RATE: usize = 10;

/// The hasher [Domain] differentiates between the modes of hashing.
///
/// The main purpose of declaring the domain is to prevent collisions between different types of
/// hashing by introducing defining differences in the way the hash function's internal state
/// (e.g. a sponge state's capacity) is initialized.
#[derive(Debug, PartialEq, Eq)]
pub enum Domain {
    /// The `VariableLength` domain is used for hashing objects that potentially serialize to more
    /// than [`RATE`] number of field elements.
    VariableLength,

    /// The `FixedLength` domain is used for hashing objects that always fit within [RATE] number
    /// of fields elements, e.g. a pair of [Digest](crate::math::digest::Digest)s.
    FixedLength,
}

/// A [cryptographic sponge][sponge]. Should only be based on a cryptographic permutation, e.g.,
/// [`Tip5`][tip5].
///
/// [sponge]: https://keccak.team/files/CSF-0.1.pdf
/// [tip5]: crate::prelude::Tip5
pub trait Sponge: Send + Sync {
    const RATE: usize;

    fn init() -> Self;

    fn absorb(&mut self, input: [BFieldElement; RATE]);

    fn squeeze(&mut self) -> [BFieldElement; RATE];

    fn pad_and_absorb_all(&mut self, input: &[BFieldElement]) {
        // pad input with [1, 0, 0, …] – padding is at least one element
        let padded_length = (input.len() + 1).next_multiple_of(RATE);
        let padding_iter =
            iter::once(&BFieldElement::ONE).chain(iter::repeat(&BFieldElement::ZERO));
        let padded_input = input.iter().chain(padding_iter).take(padded_length);

        for chunk in padded_input.chunks(RATE).into_iter() {
            // the padded input has length some multiple of `RATE`
            let absorb_elems = chunk.cloned().collect_vec().try_into().unwrap();
            self.absorb(absorb_elems);
        }
    }
}
