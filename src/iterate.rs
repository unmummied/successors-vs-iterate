use itertools::{Itertools, iterate};
use std::{collections::HashMap, iter::repeat_n};

use crate::utils::{
    ERROR_TOO_LARGE_FOR_LUCAS_LEHMER, UNREACHABLE_DIVERGENCE, UNREACHABLE_EMPTY, is_prime,
};

pub fn gcd(m: u32, n: u32) -> u32 {
    if m == 0 && n == 0 {
        return 0;
    }
    #[allow(clippy::nonminimal_bool)]
    if !(m <= n) {
        return gcd(n, m);
    }

    iterate((m, n), |&(k, d)| (d, k.checked_rem(d).unwrap_or(0)))
        .find_or_last(|&(_, q)| q == 0)
        .map_or_else(|| unreachable!("{UNREACHABLE_DIVERGENCE}"), |(res, _)| res)
}

pub fn primes() -> impl Iterator<Item = u32> {
    let mut map = HashMap::from([(4, 2)]);
    iterate((2u32, true), move |(pred, _)| {
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
        (n, is_prime)
    })
    .filter_map(|(n, is_prime)| is_prime.then_some(n))
}

pub fn collatz(n: u32) -> impl Iterator<Item = u32> {
    iterate(n, |&prev| {
        if prev.is_multiple_of(2) {
            prev / 2
        } else {
            3 * prev + 1
        }
    })
    .take_while_inclusive(|&x| 1 < x)
}

pub fn lucas_lehmer(p: u8) -> bool {
    assert!(p <= 64, "{ERROR_TOO_LARGE_FOR_LUCAS_LEHMER}");
    if p < 3 {
        return p == 2;
    }
    if !is_prime(p.into()) {
        return false;
    }

    let m_p = (1 << p) - 1;
    iterate(4u128, |&prev| {
        if prev != 0 {
            (prev * prev - 2) % m_p
        } else {
            0
        }
    })
    .nth(usize::from(p) - 2)
        == Some(0)
}

pub fn conti_frac_sqrt(n: u32) -> impl Iterator<Item = u32> {
    let a0 = n.isqrt();
    iterate(((0, 1), a0), move |&((prev_m, prev_d), prev_a)| {
        let next_m = prev_d * prev_a - prev_m;
        let next_d = (n - next_m * next_m) / prev_d;
        let next_a = (a0 + next_m) / next_d;
        ((next_m, next_d), next_a)
    })
    .take_while(|&((_, d), _)| d != 0)
    .map(|(_, a)| a)
}

pub fn detect_cycle_floyd<T, F>(x0: &T, f: &F) -> (usize, usize)
where
    T: Clone + Eq,
    F: Fn(&T) -> T,
{
    let orbit = |init| iterate(init, |x| f(x));

    let tortoise = orbit(x0.clone());
    let hare = tortoise.clone().step_by(2);
    let reunion = tortoise
        .clone()
        .zip(hare)
        .skip(1)
        .find_map(|(t, h)| (t == h).then_some(t))
        .unwrap_or_else(|| unreachable!("{UNREACHABLE_EMPTY}"));

    let cycle = orbit(f(&reunion));
    let lambda = 1 + cycle
        .enumerate()
        .find_map(|(lambda, x)| (x == reunion).then_some(lambda))
        .unwrap_or_else(|| unreachable!("{UNREACHABLE_EMPTY}"));

    let cycle = orbit(reunion);
    let mu = cycle
        .zip(tortoise)
        .enumerate()
        .find_map(|(mu, (x, t))| (x == t).then_some(mu))
        .unwrap_or_else(|| unreachable!("{UNREACHABLE_EMPTY}"));

    (mu, lambda)
}

pub fn detect_cycle_brent<T, F>(x0: &T, f: &F) -> (usize, usize)
where
    T: Clone + Eq,
    F: Fn(&T) -> T,
{
    let orbit = |init| iterate(init, |x| f(x));

    let lambda = 1
        + (0..)
            .flat_map(|exp| 0..(1 << exp))
            .zip(orbit(x0.clone()))
            .filter_map(|(offset, t)| (offset == 0).then_some(t))
            .enumerate()
            .flat_map(|(exp, t)| repeat_n(t, 1 << exp).enumerate())
            .zip(orbit(f(x0)))
            .find_map(|((lambda, t), h)| (t == h).then_some(lambda))
            .unwrap_or_else(|| unreachable!("{UNREACHABLE_EMPTY}"));

    let mu = orbit(x0.clone())
        .skip(lambda)
        .zip(orbit(x0.clone()))
        .enumerate()
        .find_map(|(mu, (x, t))| (x == t).then_some(mu))
        .unwrap_or_else(|| unreachable!("{UNREACHABLE_EMPTY}"));

    (mu, lambda)
}
