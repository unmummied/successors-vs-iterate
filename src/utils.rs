use std::iter::{from_fn, successors};

pub const ERROR_TOO_LARGE_P: &str = "Error: `p` should be `p` < 65...";
pub const ERROR_BASE_IS_LESS_THAN_2: &str = "Error: `a` should be 2 <= `a`...";
pub const UNREACHABLE_DIVERGENCE: &str =
    "Unreachable: The convergence is mathematically guaranteed.";
pub const UNREACHABLE_EMPTY: &str = "Unreachable: At the very least, `init` is `Some`.";

/// The 'dual' of [`std::iter::Iterator::fold`]
///
/// Returns an iterator that yields the first elements of the pairs returned by `f`, updating the internal state with the second elements while `f` returns `Some`.
fn unfold<A, B, F>(init: B, f: F) -> impl Iterator<Item = A>
where
    F: FnMut(B) -> Option<(A, B)>,
{
    let mut state = Some((init, f));
    from_fn(move || {
        let (curr, mut f) = state.take()?;
        let (res, next) = f(curr)?;
        state = Some((next, f));
        Some(res)
    })
}

/// Pure [`unfold`]
///
/// Similar to [Data.List.unfoldr](https://hackage-content.haskell.org/package/base-4.22.0.0/docs/Data-List.html#v:unfoldr) in Haskell
pub fn unfold_ref<A, B, F>(init: &B, f: &F) -> impl Iterator<Item = A>
where
    F: Fn(&B) -> Option<(A, B)>,
{
    successors(f(init), |(_, b)| f(b)).map(|(a, _)| a)
}

/// Modular exponentiation
///
/// Returns the remainder when an integer `base` is raised to the power `exp`, and then divided by a positive integer `modulo`.
fn pow_mod(base: u32, exp: u32, modulo: u32) -> Option<u32> {
    if modulo == 0 {
        return None;
    }

    let (base, modulo) = (u64::from(base), u64::from(modulo));

    let (res, _) = unfold_ref(&exp, &|&bits| {
        (bits != 0).then_some((bits & 1 == 1, bits >> 1))
    })
    .fold((1 % modulo, base % modulo), |(acc, b), p| {
        let acc = if p { acc * b % modulo } else { acc };
        (acc, (b * b) % modulo)
    });
    res.try_into().ok()
}

/// 'Naive' primality test
///
/// Returns whether `n` is prime.
pub const fn is_prime(n: u32) -> bool {
    if n <= 2 {
        return n == 2;
    }
    if n.is_multiple_of(2) {
        return false;
    }
    let mut d = 3;
    while d <= n.isqrt() {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_pow_mod() {
        fn naive_pow_mod(mut base: u64, mut exp: u64, modulo: u64) -> u64 {
            let mut res = 1 % modulo;
            base %= modulo;

            while exp > 0 {
                if exp & 1 == 1 {
                    res = (res * base) % modulo;
                }
                base = (base * base) % modulo;
                exp >>= 1;
            }

            res
        }

        assert_eq!(pow_mod(0, 0, 0), None);
        assert_eq!(pow_mod(0, 0, 1), Some(0));
        assert_eq!(pow_mod(0, 0, 7), Some(1)); // convention check
        assert_eq!(pow_mod(5, 1, 7), Some(5));
        assert_eq!(pow_mod(5, 2, 1), Some(0));
        assert_eq!(pow_mod(2, 10, 1000), Some(1024 % 1000));
        assert_eq!(pow_mod(3, 0, 5), Some(1));
        assert_eq!(pow_mod(0, 5, 7), Some(0));
        assert_eq!(pow_mod(1, 999, 13), Some(1));
        assert_eq!(pow_mod(2, 3, 5), Some(3));
        assert_eq!(pow_mod(0, 0, 1), Some(0));
        assert_eq!(pow_mod(3, 4, 7), Some(4));
        assert_eq!(pow_mod(5, 3, 13), Some(8));
        assert_eq!(pow_mod(7, 0, 11), Some(1));
        assert_eq!(pow_mod(0, 5, 10), Some(0));
        assert_eq!(pow_mod(0, 0, 10), Some(1));
        assert_eq!(pow_mod(5, 3, 1), Some(0));
        assert_eq!(pow_mod(10, 2, 3), Some(1));
        assert_eq!(pow_mod(14, 3, 5), Some(4));
        assert_eq!(pow_mod(2, 1_000_000, 13), Some(3));
        assert_eq!(pow_mod(3, 1_000_000_000, 17), Some(1));
        assert_eq!(
            pow_mod(123_456_789, 987_654_321, 1_000_000_007),
            Some(652_541_198)
        );
        assert_eq!(pow_mod(4_294_967_295, 2, 4_294_967_291), Some(16));
        assert_eq!(pow_mod(5, 10, 0), None);
        assert_eq!(pow_mod(10, 10, 2), Some(0));
        assert_eq!(pow_mod(10, 10, 3), Some(1));
        assert_eq!(
            pow_mod(123, 456, 97),
            Some(naive_pow_mod(123, 456, 97) as _)
        );
        assert_eq!(
            pow_mod(2, 20, 1_000),
            Some(naive_pow_mod(2, 20, 1_000) as _)
        );
        assert_eq!(
            pow_mod(7, 31, 1_000_000_007),
            Some(naive_pow_mod(7, 31, 1_000_000_007) as _)
        );
        assert_eq!(
            pow_mod(12345, 6, 97),
            Some(naive_pow_mod(12345, 6, 97) as _)
        );

        assert_eq!(
            pow_mod(99991, 12345, 1_000_000_007),
            Some(naive_pow_mod(99991, 12345, 1_000_000_007) as _)
        );
    }

    #[test]
    fn test_is_prime() {
        let primes = [
            // A000040
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
            181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
        ];
        let odd_composites = [
            // A071904
            9, 15, 21, 25, 27, 33, 35, 39, 45, 49, 51, 55, 57, 63, 65, 69, 75, 77, 81, 85, 87, 91,
            93, 95, 99, 105, 111, 115, 117, 119, 121, 123, 125, 129, 133, 135, 141, 143, 145, 147,
            153, 155, 159, 161, 165, 169, 171, 175, 177, 183, 185, 187, 189, 195, 201, 203, 205,
        ];

        for p in primes {
            assert!(is_prime(p));
        }
        for c in odd_composites {
            assert!(!is_prime(c));
        }
    }
}
