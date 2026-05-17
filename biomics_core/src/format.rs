//! Zero-allocation integer and floating-point formatting utilities.
//!
//! These wrappers expose the `itoa` and `ryu` fast formatting crates through a
//! single, consistent API. Callers stack-allocate the `Buffer` once and reuse
//! it across iterations to avoid any heap traffic.

/// Format a `u64` into a caller-supplied `itoa::Buffer` with zero heap allocation.
///
/// The returned `&str` is borrowed from `buf` and valid until the next call to
/// `fmt_u64` (or any other `itoa::Buffer` method) on the same buffer.
#[inline]
pub fn fmt_u64(n: u64, buf: &mut itoa::Buffer) -> &str {
    buf.format(n)
}

/// Format a `f64` into a caller-supplied `ryu::Buffer` with zero heap allocation.
///
/// Uses Ryū's shortest-round-trip algorithm, which produces the fewest digits
/// needed to uniquely identify the value.
#[inline]
pub fn fmt_f64(n: f64, buf: &mut ryu::Buffer) -> &str {
    buf.format(n)
}
