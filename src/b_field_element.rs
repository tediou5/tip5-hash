use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{ConstOne, ConstZero, One, Zero};

/// Base field element ∈ ℤ_{2^64 - 2^32 + 1}.
///
/// In Montgomery representation. This implementation follows <https://eprint.iacr.org/2022/274.pdf>
/// and <https://github.com/novifinancial/winterfell/pull/101/files>.
#[derive(Debug, Copy, Clone, Default, Hash, PartialEq, Eq)]
pub struct BFieldElement(u64);

impl BFieldElement {
    pub const BYTES: usize = 8;

    /// The base field's prime, _i.e._, 2^64 - 2^32 + 1.
    pub const P: u64 = 0xffff_ffff_0000_0001;
    pub const MAX: u64 = Self::P - 1;

    /// 2^128 mod P; this is used for conversion of elements into Montgomery representation.
    const R2: u64 = 0xffff_fffe_0000_0001;

    #[inline]
    pub const fn new(value: u64) -> Self {
        Self(Self::montyred((value as u128) * (Self::R2 as u128)))
    }

    /// Montgomery reduction
    #[inline(always)]
    pub const fn montyred(x: u128) -> u64 {
        // See reference above for a description of the following implementation.
        let xl = x as u64;
        let xh = (x >> 64) as u64;
        let (a, e) = xl.overflowing_add(xl << 32);

        let b = a.wrapping_sub(a >> 32).wrapping_sub(e as u64);

        let (r, c) = xh.overflowing_sub(b);

        // See https://github.com/Neptune-Crypto/twenty-first/pull/70 for various ways
        // of expressing this.
        r.wrapping_sub((1 + !Self::P) * c as u64)
    }

    #[inline]
    pub const fn value(&self) -> u64 {
        self.canonical_representation()
    }

    #[inline]
    const fn canonical_representation(&self) -> u64 {
        Self::montyred(self.0 as u128)
    }
}

impl BFieldElement {
    #[must_use]
    #[inline]
    pub fn inverse(&self) -> Self {
        #[inline(always)]
        const fn exp(base: BFieldElement, exponent: u64) -> BFieldElement {
            let mut res = base;
            let mut i = 0;
            while i < exponent {
                res = BFieldElement(BFieldElement::montyred(res.0 as u128 * res.0 as u128));
                i += 1;
            }
            res
        }

        let x = *self;
        assert_ne!(
            x,
            Self::zero(),
            "Attempted to find the multiplicative inverse of zero."
        );

        let bin_2_ones = x.square() * x;
        let bin_3_ones = bin_2_ones.square() * x;
        let bin_6_ones = exp(bin_3_ones, 3) * bin_3_ones;
        let bin_12_ones = exp(bin_6_ones, 6) * bin_6_ones;
        let bin_24_ones = exp(bin_12_ones, 12) * bin_12_ones;
        let bin_30_ones = exp(bin_24_ones, 6) * bin_6_ones;
        let bin_31_ones = bin_30_ones.square() * x;
        let bin_31_ones_1_zero = bin_31_ones.square();
        let bin_32_ones = bin_31_ones.square() * x;

        exp(bin_31_ones_1_zero, 32) * bin_32_ones
    }

    #[inline(always)]
    fn square(self) -> Self {
        self * self
    }

    /// Return the raw bytes or 8-bit chunks of the Montgomery
    /// representation, in little-endian byte order
    pub const fn raw_bytes(&self) -> [u8; 8] {
        self.0.to_le_bytes()
    }

    #[inline]
    pub const fn raw_u64(&self) -> u64 {
        self.0
    }

    /// Take a slice of 8 bytes and interpret it as an integer in
    /// little-endian byte order, and cast it to a BFieldElement
    /// in Montgomery representation
    pub const fn from_raw_bytes(bytes: &[u8; 8]) -> Self {
        Self(u64::from_le_bytes(*bytes))
    }

    #[inline]
    pub const fn from_raw_u64(e: u64) -> BFieldElement {
        BFieldElement(e)
    }
}

impl Zero for BFieldElement {
    #[inline]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self == &Self::ZERO
    }
}

impl ConstZero for BFieldElement {
    const ZERO: Self = Self::new(0);
}

impl One for BFieldElement {
    #[inline]
    fn one() -> Self {
        Self::ONE
    }

    #[inline]
    fn is_one(&self) -> bool {
        self == &Self::ONE
    }
}

impl ConstOne for BFieldElement {
    const ONE: Self = Self::new(1);
}

impl Add for BFieldElement {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        // Compute a + b = a - (p - b).
        let (x1, c1) = self.0.overflowing_sub(Self::P - rhs.0);

        // The following if/else is equivalent to the commented-out code below but
        // the if/else was found to be faster.
        // let adj = 0u32.wrapping_sub(c1 as u32);
        // Self(x1.wrapping_sub(adj as u64))
        // See
        // https://github.com/Neptune-Crypto/twenty-first/pull/70
        if c1 {
            Self(x1.wrapping_add(Self::P))
        } else {
            Self(x1)
        }
    }
}

impl AddAssign for BFieldElement {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl SubAssign for BFieldElement {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl MulAssign for BFieldElement {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul for BFieldElement {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self(Self::montyred((self.0 as u128) * (rhs.0 as u128)))
    }
}

impl Neg for BFieldElement {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::zero() - self
    }
}

impl Sub for BFieldElement {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let (x1, c1) = self.0.overflowing_sub(rhs.0);

        // The following code is equivalent to the commented-out code below
        // but they were determined to have near-equiavalent running times. Maybe because
        // subtraction is not used very often.
        // See: https://github.com/Neptune-Crypto/twenty-first/pull/70
        // 1st alternative:
        // if c1 {
        //     Self(x1.wrapping_add(Self::P))
        // } else {
        //     Self(x1)
        // }
        // 2nd alternative:
        // let adj = 0u32.wrapping_sub(c1 as u32);
        // Self(x1.wrapping_sub(adj as u64))
        Self(x1.wrapping_sub((1 + !Self::P) * c1 as u64))
    }
}

impl Div for BFieldElement {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, other: Self) -> Self {
        other.inverse() * self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_montgomery_reduction() {
        let input = 2_609_026_890_597_981_882u128;
        let expected = 11_259_563_268_822_605_859u64;

        let result = BFieldElement::montyred(input);

        assert_eq!(result, expected, "Montgomery reduction failed");
    }

    #[test]
    fn test_montify() {
        let value = 12_045_832_659_793_544_965;
        let bfe = BFieldElement::new(value);
        let expected = 9_712_864_734_344_745_984u64;

        let result = BFieldElement::new(value).0;

        assert_eq!(result, expected, "Montification failed");

        let red = bfe.value();
        assert_eq!(red, value, "Canonical representation failed");

    }
}