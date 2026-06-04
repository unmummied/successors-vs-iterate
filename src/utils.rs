use std::iter::{from_fn, successors};

pub const ERROR_TOO_LARGE_FOR_LUCAS_LEHMER: &str =
    "Error: `p` should be less than `65`... see https://oeis.org/A000043";
pub const UNREACHABLE_DIVERGENCE: &str =
    "Unreachable: The convergence is mathematically guaranteed...";
pub const UNREACHABLE_EMPTY: &str = "Unreachable: At the very least, `init` is `Some`...";

/// 'Naive' primality test
///
/// Returns whether `n` is prime.
pub const fn is_prime(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n.is_multiple_of(2) {
        return n == 2;
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

/// Similar function of the `Data.List.unfoldr` in Haskell
///
/// The 'dual' to fold: while fold reduces a list to a summary value, unfold builds a list from a seed value.
fn unfold_ref<A, B, F>(init: &B, f: &F) -> impl Iterator<Item = A>
where
    F: Fn(&B) -> Option<(A, B)>,
{
    successors(f(init), |(_, b)| f(b)).map(|(a, _)| a)
}

/// Mutable and ownership-consuming version of `unfold_ref`
///
/// This version takes ownership of the initial seed value and allows the generator function to maintain and mutate state directly via `FnMut`.
fn unfold<A, B, F>(init: B, mut f: F) -> impl Iterator<Item = A>
where
    F: FnMut(&mut B) -> Option<(A, B)>,
{
    let mut state = init;
    from_fn(move || {
        let (res, next) = f(&mut state)?;
        state = next;
        Some(res)
    })
}

/// Modular exponentiation
///
/// Returns the remainder when an integer `base` is raised to the power `exp`, and then divided by a positive integer `modulo`.
fn pow_mod(base: u32, exp: u32, modulo: u32) -> Option<u32> {
    if modulo == 0 {
        return None;
    }

    let (base, modulo) = (u64::from(base), u64::from(modulo));

    let (res, _) = unfold(exp, |&mut bits| {
        (bits != 0).then_some((bits & 1 == 1, bits >> 1))
    })
    .fold((1 % modulo, base % modulo), |(acc, b), p| {
        let acc = if p { (acc * b) % modulo } else { acc };
        (acc, (b * b) % modulo)
    });
    res.try_into().ok()
}

/// R. P. Brent's improved cycle detection algorithm
///
/// Returns `mu` (the length of the tail) and `lambda` (the length of the cycle).
fn detect_cycle_brent<T, F>(x0: &T, f: &F) -> (usize, usize)
where
    T: Clone + Eq,
    F: Fn(T) -> T,
{
    let mut tortoise = x0.clone();
    let mut hare = x0.clone();

    let mut power = 1;
    let mut lambda = 1;

    hare = f(hare);
    while tortoise != hare {
        if power == lambda {
            tortoise = hare.clone();
            power *= 2;
            lambda = 0;
        }
        hare = f(hare);
        lambda += 1;
    }

    tortoise = x0.clone();
    hare = x0.clone();
    for _ in 0..lambda {
        hare = f(hare);
    }

    let mut mu = 0;
    while tortoise != hare {
        tortoise = f(tortoise);
        hare = f(hare);
        mu += 1;
    }

    (mu, lambda)
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
}
