//! Actually only provides the [ZMod2] struct and corresponding methods to 
//! represent integers modulo 2 and to perform modulo 2 arithmetic on them.

pub use std::fmt;
pub use num_traits::{Zero, One};
pub use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};


/// Used to represent congruence classes of the ring of integers modulo 2 with
/// the corresponding arithmetic.
/// Note: only implements the methods and traits that are necessary for this
///       assignment.
#[derive(Clone,PartialEq, Copy)]
pub struct ZMod2 {
    #[doc(hidden)]
    val: bool
}

impl From<u32> for ZMod2 {
    /// Maps an unsigned integer `value` to its congruence class of `value` modulo 2.
    fn from(value: u32) -> Self {
        ZMod2 {
            val: value % 2 != 0,
        }
    }
}

impl From<ZMod2> for u32 {
    /// Maps a congruence class `value` modulo 2 to the representative integer of
    /// that class.
    fn from(value: ZMod2) -> Self {
        value.val as u32
    }
}

impl fmt::Display for ZMod2 {
    /// Formats the congruence class `self` modulo 2 by the formatting of its
    /// representative integer.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val as i32)
    }
}

impl fmt::Debug for ZMod2 {
    /// Formats the congruence class `self` modulo 2 by the formatting of its
    /// representative integer.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val as i32)
    }
}

impl Zero for ZMod2 {
    fn zero() -> Self{
        ZMod2{val: false}
    }

    fn is_zero(&self) -> bool {
        !self.val
    }

    fn set_zero(&mut self) {
        self.val = false
    }
}

impl One for ZMod2 {
    fn one() -> Self{
        ZMod2{val: true}
    }

    fn is_one(&self) -> bool {
        self.val
    }

    fn set_one(&mut self) {
        self.val = true
    }
}

impl Add for ZMod2 {
    type Output = ZMod2;

    fn add(self, rhs: Self) -> Self::Output {
        ZMod2{val: self.val ^ rhs.val}
    }
}

impl AddAssign for ZMod2 {
    fn add_assign(&mut self, rhs: Self) {
        self.val ^= rhs.val;
    }
}

impl Sub for ZMod2 {
    type Output = ZMod2;

    fn sub(self, rhs: Self) -> Self::Output {
        ZMod2{val: self.val ^ rhs.val}
    }
}

impl SubAssign for ZMod2 {
    fn sub_assign(&mut self, rhs: Self) {
        self.val ^= rhs.val;
    }
}

impl Mul for ZMod2 {
    type Output = ZMod2;

    fn mul(self, rhs: Self) -> Self::Output {
        ZMod2{val: self.val && rhs.val}
    }
}

impl MulAssign for ZMod2 {
    fn mul_assign(&mut self, rhs: Self) {
        self.val = self.val && rhs.val;
    }
}
