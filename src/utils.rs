use anyhow::anyhow;
use std::iter::{from_fn, successors};

pub const ERROR_TOO_LARGE_P: &str = "Error: 'p' must satisfy p < 65...";
pub const UNREACHABLE_DIVERGENCE: &str =
    "Unreachable: The convergence is mathematically guaranteed.";
pub const UNREACHABLE_EMPTY: &str = "Unreachable: At the very least, 'init' is a Some.";

/// Miller-Rabin primality test
///
/// Returns whether `n` is prime.
///
/// # Deterministic test against tiny n
///
/// Generally, trying all bases such that a < 2 ln^2 n is sufficient under the generalized Riemann hypothesis.
/// However, it is not necessary to try every base up to the bound for small `n`, as much smaller sets of potential witnesses are known to suffice.
///
/// For example:
/// - If `n` is an [`u32`], the set of witnesses {2, 7, 61} is sufficient.
/// - If `n` is an [`u64`], the set of witnesses {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} is sufficient.
pub fn is_prime(n: u32) -> bool {
    if n < 2 {
        return false;
    }

    [2, 7, 61]
        .into_iter()
        .take_while(|&a| a < n)
        .all(|a| is_strong_probable_prime(n, a).unwrap_or(false))
}

/// Strong probable primality test
///
/// Returns whether `n` is a strong probable prime on base `a`.
///
/// # Errors
///
/// Returns an error unless:
/// - `3 <= n`
/// - `n` is odd
/// - `2 <= a < n`
fn is_strong_probable_prime(n: u32, a: u32) -> anyhow::Result<bool> {
    const ERROR_INVALID_N: &str = "Error: 'n' must be odd and greater than 2...";
    const ERROR_INVALID_A: &str = "Error: 'a' must satisfy 2 <= a < n...";

    #[allow(clippy::nonminimal_bool)]
    if !(3 <= n) || n.is_multiple_of(2) {
        return Err(anyhow!(ERROR_INVALID_N));
    }
    #[allow(clippy::nonminimal_bool)]
    if !(2 <= a) || !(a < n) {
        return Err(anyhow!(ERROR_INVALID_A));
    }

    // m = d * 2^s
    let m = n - 1;
    let s = m.trailing_zeros();
    let d = m >> s;

    let mut x = mod_pow(a, d, n)?;

    if x == 1 || x == m {
        return Ok(true);
    }

    for _ in 1..s {
        x = mod_pow(x, 2, n)?;
        if x == m {
            return Ok(true);
        }
        if x == 1 {
            return Ok(false);
        }
    }

    Ok(false)
}

/// Modular exponentiation
///
/// Returns the remainder when an integer `base` is raised to the power `exp`, and then divided by a positive integer `modulo`.
///
/// # Hylomorphism
///
/// Conceptually, this implementation is a hylomorphism composed of an anamorphism ([`unfold_ref`]) and a catamorphism ([`fold`]).
///
/// Since the intermediate sequence is represented as a lazy [`Iterator`], it is consumed as it is produced, achiving deforestation.
fn mod_pow(base: u32, exp: u32, modulo: u32) -> anyhow::Result<u32> {
    const ERROR_ZERO_MODULO: &str = "Error: modulo must be non-zero...";
    if modulo == 0 {
        return Err(anyhow!(ERROR_ZERO_MODULO));
    }

    let (base, modulo) = (u64::from(base), u64::from(modulo));

    let (res, _) = unfold_ref(&exp, &|&bits| {
        (bits != 0).then_some((bits & 1 == 1, bits >> 1))
    })
    .fold((1 % modulo, base % modulo), |(acc, b), p| {
        let acc = if p { acc * b % modulo } else { acc };
        (acc, (b * b) % modulo)
    });
    Ok(res.try_into()?)
}

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
fn unfold_ref<A, B, F>(init: &B, f: &F) -> impl Iterator<Item = A>
where
    F: Fn(&B) -> Option<(A, B)>,
{
    successors(f(init), |(_, b)| f(b)).map(|(a, _)| a)
}

#[cfg(test)]
mod tests {
    use std::convert::identity;

    use super::*;

    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_pow_mod() -> anyhow::Result<()> {
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

        assert!(mod_pow(0, 0, 0).is_err());
        assert_eq!(mod_pow(0, 0, 1)?, 0);
        assert_eq!(mod_pow(0, 0, 7)?, 1);
        assert_eq!(mod_pow(5, 1, 7)?, 5);
        assert_eq!(mod_pow(5, 2, 1)?, 0);
        assert_eq!(mod_pow(2, 10, 1000)?, 1024 % 1000);
        assert_eq!(mod_pow(3, 0, 5)?, 1);
        assert_eq!(mod_pow(0, 5, 7)?, 0);
        assert_eq!(mod_pow(1, 999, 13)?, 1);
        assert_eq!(mod_pow(2, 3, 5)?, 3);
        assert_eq!(mod_pow(0, 0, 1)?, 0);
        assert_eq!(mod_pow(3, 4, 7)?, 4);
        assert_eq!(mod_pow(5, 3, 13)?, 8);
        assert_eq!(mod_pow(7, 0, 11)?, 1);
        assert_eq!(mod_pow(0, 5, 10)?, 0);
        assert_eq!(mod_pow(0, 0, 10)?, 1);
        assert_eq!(mod_pow(5, 3, 1)?, 0);
        assert_eq!(mod_pow(10, 2, 3)?, 1);
        assert_eq!(mod_pow(14, 3, 5)?, 4);
        assert_eq!(mod_pow(2, 1_000_000, 13)?, 3);
        assert_eq!(mod_pow(3, 1_000_000_000, 17)?, 1);
        assert_eq!(
            mod_pow(123_456_789, 987_654_321, 1_000_000_007)?,
            652_541_198
        );
        assert_eq!(mod_pow(4_294_967_295, 2, 4_294_967_291)?, 16);
        assert!(mod_pow(5, 10, 0).is_err());
        assert_eq!(mod_pow(10, 10, 2)?, 0);
        assert_eq!(mod_pow(10, 10, 3)?, 1);
        assert_eq!(mod_pow(123, 456, 97)?, naive_pow_mod(123, 456, 97) as _);
        assert_eq!(mod_pow(2, 20, 1_000)?, naive_pow_mod(2, 20, 1_000) as _);
        assert_eq!(
            mod_pow(7, 31, 1_000_000_007)?,
            naive_pow_mod(7, 31, 1_000_000_007) as _
        );
        assert_eq!(mod_pow(12345, 6, 97)?, naive_pow_mod(12345, 6, 97) as _);
        assert_eq!(
            mod_pow(99991, 12345, 1_000_000_007)?,
            naive_pow_mod(99991, 12345, 1_000_000_007) as _
        );

        Ok(())
    }

    #[test]
    fn test_is_prime() {
        assert!(!is_prime(0));

        // A000040
        let primes = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
            181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
        ];
        // A018252
        let odd_composites = [
            1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34,
            35, 36, 38, 39, 40, 42, 44, 45, 46, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 60, 62, 63,
            64, 65, 66, 68, 69, 70, 72, 74, 75, 76, 77, 78, 80, 81, 82, 84, 85, 86, 87, 88,
        ];

        for p in primes {
            assert!(is_prime(p));
        }
        for c in odd_composites {
            assert!(!is_prime(c));
        }
    }

    #[test]
    fn test_is_strong_probable_prime() {
        assert!(is_strong_probable_prime(3, 2).is_ok_and(identity));
        assert!(is_strong_probable_prime(5, 2).is_ok_and(identity));

        assert!(is_strong_probable_prime(3, 3).is_err());
        assert!(is_strong_probable_prime(5, 5).is_err());

        assert!(is_strong_probable_prime(3, 1).is_err());
        assert!(is_strong_probable_prime(1, 2).is_err());
        assert!(is_strong_probable_prime(2, 2).is_err());

        // A000040
        let primes = [
            // 2, 3, 5, 7, 11, 13, 17,
            19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
            109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
            199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
        ];
        // A071904
        let odd_composites = [
            // 9,
            15, 21, 25, 27, 33, 35, 39, 45, 49, 51, 55, 57, 63, 65, 69, 75, 77, 81, 85, 87,
            // 91,
            93, 95, 99, 105, 111, 115, 117, 119, // 121,
            123, 125, 129, 133, 135, 141, 143, // 145,
            147, 153, 155, 159, 161, 165, 169, 171, 175, 177, 183, 185, 187, 189, 195, 201, 203,
            205,
        ];

        // A001262
        let sp2 = [
            2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665, 80581,
            85489, 88357, 90751, 104_653, 130_561, 196_093, 220_729, 233_017, 252_601, 253_241,
            256_999, 271_951, 280_601, 314_821, 357_761, 390_937, 458_989, 476_971, 486_737,
        ];
        // A020229
        let sp3 = [
            // 121,
            703, 1891, 3281, 8401, 8911, 10585, 12403, 16531, 18721, 19345, 23521, 31621, 44287,
            47197, 55969, 63139, 74593, 79003, 82513, 87913, 88573, 97567, 105_163, 111_361,
            112_141, 148_417, 152_551, 182_527, 188_191, 211_411, 218_791, 221_761, 226_801,
        ];
        // A020230
        let sp4 = [
            341, 1387, 2047, 3277, 4033, 4371, 4681, 5461, 8321, 8911, 10261, 13747, 14491, 15709,
            15841, 19951, 29341, 31621, 42799, 49141, 49981, 52633, 60787, 65077, 65281, 74665,
            80581, 83333, 85489, 88357, 90751, 104_653, 123_251, 129_921, 130_561, 137_149,
        ];
        // A020231
        let sp5 = [
            781, 1541, 5461, 5611, 7813, 13021, 14981, 15751, 24211, 25351, 29539, 38081, 40501,
            44801, 53971, 79381, 100_651, 102_311, 104_721, 112_141, 121_463, 133_141, 141_361,
            146_611, 195_313, 211_951, 216_457, 222_301, 251_521, 289_081, 290_629, 298_271,
            315_121,
        ];
        // A020232
        let sp6 = [
            217, 481, 1111, 1261, 2701, 3589, 5713, 6533, 11041, 14701, 20017, 29341, 34441, 39493,
            43621, 46657, 46873, 49141, 49661, 58969, 74023, 74563, 76921, 83333, 87061, 92053,
            94657, 94697, 97751, 97921, 109_061, 115_921, 125_563, 128_627, 151_387, 173_377,
        ];
        // A020243
        let sp17 = [
            // 9, 91, 145,
            781, 1111, 2821, 4033, 4187, 5365, 5833, 6697, 7171, 15805, 19729, 21781, 22791, 24211,
            26245, 31621, 33001, 33227, 34441, 35371, 38081, 42127, 49771, 71071, 74665, 77293,
            78881, 88831, 96433, 97921, 98671, 101_101, 102_311, 125_563, 129_493,
        ];

        for p in primes {
            assert!(is_strong_probable_prime(p, 2).is_ok_and(identity));
            assert!(is_strong_probable_prime(p, 3).is_ok_and(identity));
            assert!(is_strong_probable_prime(p, 4).is_ok_and(identity));
            assert!(is_strong_probable_prime(p, 5).is_ok_and(identity));
            assert!(is_strong_probable_prime(p, 6).is_ok_and(identity));
            assert!(is_strong_probable_prime(p, 17).is_ok_and(identity));
        }
        for c in odd_composites {
            assert!(!is_strong_probable_prime(c, 2).is_ok_and(identity));
            assert!(!is_strong_probable_prime(c, 3).is_ok_and(identity));
            assert!(!is_strong_probable_prime(c, 4).is_ok_and(identity));
            assert!(!is_strong_probable_prime(c, 5).is_ok_and(identity));
            assert!(!is_strong_probable_prime(c, 6).is_ok_and(identity));
            assert!(!is_strong_probable_prime(c, 17).is_ok_and(identity));
        }

        for (sp, a) in itertools::chain!(
            sp2.map(|sp| (sp, 2)),
            sp3.map(|sp| (sp, 3)),
            sp4.map(|sp| (sp, 4)),
            sp5.map(|sp| (sp, 5)),
            sp6.map(|sp| (sp, 6)),
            sp17.map(|sp| (sp, 17)),
        ) {
            assert!(is_strong_probable_prime(sp, a).is_ok_and(identity));
        }
    }
}
