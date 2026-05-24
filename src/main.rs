mod iterate;
mod successors;
mod utils;

fn main() {
    println!("Hello, world!");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd() {
        let test = |f: fn(_, _) -> _| {
            assert_eq!(f(48, 18), 6);
            assert_eq!(f(101, 10), 1);
            assert_eq!(f(0, 5), 5);
            assert_eq!(f(5, 0), 5);
            assert_eq!(f(0, 0), 0);
            assert_eq!(f(1, 10), 1);
            assert_eq!(f(10, 1), 1);
            assert_eq!(f(7, 7), 7);
            assert_eq!(f(13 * 100, 17 * 100), 100);
        };
        test(iterate::gcd);
        test(successors::gcd);
    }

    #[test]
    fn test_primes() {
        let naive = (2..1_000).filter(|&n| utils::is_prime(n));
        let successors = successors::primes();
        let iterate = iterate::primes();
        assert!(
            successors
                .zip(iterate)
                .zip(naive)
                .into_iter()
                .all(|((s, i), n)| s == i && i == n)
        );
    }

    #[test]
    fn test_collatz() {
        fn test<I: IntoIterator<Item = u32>>(f: fn(u32) -> I) {
            assert!(f(6).into_iter().eq([6, 3, 10, 5, 16, 8, 4, 2, 1]));
            assert!(
                f(11)
                    .into_iter()
                    .eq([11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1,])
            );
            assert!(f(27).into_iter().eq([
                27, 82, 41, 124, 62, 31, 94, 47, 142, 71, 214, 107, 322, 161, 484, 242, 121, 364,
                182, 91, 274, 137, 412, 206, 103, 310, 155, 466, 233, 700, 350, 175, 526, 263, 790,
                395, 1186, 593, 1780, 890, 445, 1336, 668, 334, 167, 502, 251, 754, 377, 1132, 566,
                283, 850, 425, 1276, 638, 319, 958, 479, 1438, 719, 2158, 1079, 3238, 1619, 4858,
                2429, 7288, 3644, 1822, 911, 2734, 1367, 4102, 2051, 6154, 3077, 9232, 4616, 2308,
                1154, 577, 1732, 866, 433, 1300, 650, 325, 976, 488, 244, 122, 61, 184, 92, 46, 23,
                70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1,
            ]));
        }
        test(iterate::collatz);
        test(successors::collatz);
    }

    #[test]
    fn test_lucas_lehmer() {
        let test = |f: fn(_) -> _| {
            assert!([2, 3, 5, 7, 13, 17, 19, 31, 61].into_iter().all(&f));
            assert!(![11, 23, 29, 37, 41, 43, 47, 53, 59].into_iter().any(&f));
        };
        test(iterate::lucas_lehmer);
        test(successors::lucas_lehmer);
    }
}
