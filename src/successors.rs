use std::{collections::HashMap, iter::successors};

use crate::utils::{ERROR_TOO_LARGE_FOR_LUCAS_LEHMER, UNREACHABLE_DIVERGENCE, is_prime};

/// Greatest common divisor
///
/// Returns the largest positive integer that divides each of the integers.
pub fn gcd(m: u32, n: u32) -> u32 {
    if m == 0 && n == 0 {
        return 0; // see https://en.wikipedia.org/wiki/B%C3%A9zout's_identity
    }
    #[allow(clippy::nonminimal_bool)]
    if !(m <= n) {
        return gcd(n, m);
    }
    successors(Some((m, n)), |&(k, d)| k.checked_rem(d).map(|q| (d, q)))
        .last()
        .map_or_else(|| unreachable!("{UNREACHABLE_DIVERGENCE}"), |(res, _)| res)
}

/// Prime Generator
///
/// Returns an iterator that yields the prime numbers in ascending order up to `u32::MAX`.
pub fn primes() -> impl Iterator<Item = u32> {
    let mut map = HashMap::from([(4, 2)]);
    successors(Some((2u32, true)), move |(pred, _)| {
        let n = pred + 1;
        let is_prime = match map.remove(&n) {
            None => {
                if let Some(square) = n.checked_mul(n) {
                    map.insert(square, n);
                }
                true
            }
            Some(p) => {
                let mut skipped = n + p;
                while map.contains_key(&skipped) {
                    skipped += p;
                }
                map.insert(skipped, p);
                false
            }
        };
        Some((n, is_prime))
    })
    .filter_map(|(n, is_prime)| is_prime.then_some(n))
}

/// Collatz orbit
///
/// Returns an iterator that yields the Collatz orbit starting from `n`.
///
/// The iterator terminates once it reaches 1.
pub fn collatz(n: u32) -> impl Iterator<Item = u32> {
    successors(Some(n), |&prev| {
        (1 < prev).then(|| {
            if prev.is_multiple_of(2) {
                prev / 2
            } else {
                3 * prev + 1
            }
        })
    })
}

/// Lucas-Lehmer primality test
///
/// Returns whether the `p`-th Mersenne number is prime, where `p` is a prime number.
pub fn lucas_lehmer(p: u8) -> bool {
    assert!(p <= 64, "{ERROR_TOO_LARGE_FOR_LUCAS_LEHMER}");
    if p < 3 {
        return p == 2;
    }
    if !is_prime(p.into()) {
        // M_p is composite if p is composite.
        return false;
    }
    let m_p = (1 << p) - 1;
    successors(Some(4u128), |&prev| {
        (prev != 0).then(|| (prev * prev - 2) % m_p)
    })
    .nth(usize::from(p) - 2)
        == Some(0)
}

/// The successive convergents of the continued fraction of square root
///
/// Returns an iterator that yields the partial quotients (coefficients) of the continued fraction expansion of square root of `n`.
///
/// If `n` is a (perfect) square number, it yields `n.isqrt()` and terminates, otherwise, it loops infinitely.
pub fn conti_frac_sqrt(n: u32) -> impl Iterator<Item = u32> {
    let a0 = n.isqrt();
    successors(Some(((0, 1), a0)), move |&((prev_m, prev_d), prev_a)| {
        let next_m = prev_d * prev_a - prev_m;
        let next_d = (n - next_m * next_m).checked_div(prev_d)?;
        let next_a = (a0 + next_m).checked_div(next_d)?;
        Some(((next_m, next_d), next_a))
    })
    .map(|(_, a)| a)
}
