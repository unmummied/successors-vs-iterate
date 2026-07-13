mod iterate;
mod successors;
mod utils;

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
            assert_eq!(f(60, 60), 60);
            assert_eq!(f(1071, 462), 21);
            assert_eq!(f(13 * 100, 17 * 100), 100);
        };
        test(successors::gcd);
        test(iterate::gcd);
    }

    #[test]
    fn test_primes() {
        let naive = (0..1_000).filter(|&n| utils::is_prime(n));
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
        fn test<F, I>(f: F)
        where
            F: Fn(u32) -> I,
            I: Iterator<Item = u32>,
        {
            assert!(f(6).eq([6, 3, 10, 5, 16, 8, 4, 2, 1]));
            assert!(f(11).eq([11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1,]));
            assert!(f(27).eq([
                27, 82, 41, 124, 62, 31, 94, 47, 142, 71, 214, 107, 322, 161, 484, 242, 121, 364,
                182, 91, 274, 137, 412, 206, 103, 310, 155, 466, 233, 700, 350, 175, 526, 263, 790,
                395, 1186, 593, 1780, 890, 445, 1336, 668, 334, 167, 502, 251, 754, 377, 1132, 566,
                283, 850, 425, 1276, 638, 319, 958, 479, 1438, 719, 2158, 1079, 3238, 1619, 4858,
                2429, 7288, 3644, 1822, 911, 2734, 1367, 4102, 2051, 6154, 3077, 9232, 4616, 2308,
                1154, 577, 1732, 866, 433, 1300, 650, 325, 976, 488, 244, 122, 61, 184, 92, 46, 23,
                70, 35, 106, 53, 160, 80, 40, 20, 10, 5, 16, 8, 4, 2, 1,
            ]));
        }
        test(successors::collatz);
        test(iterate::collatz);
    }

    #[test]
    fn test_is_mersenne_exp() {
        let test = |f: fn(_) -> _| {
            assert!([2, 3, 5, 7, 13, 17, 19, 31, 61].into_iter().all(&f));
            assert!(![11, 23, 29, 37, 41, 43, 47, 53, 59].into_iter().any(&f));
        };
        test(successors::is_mersenne_exp);
        test(iterate::is_mersenne_exp);
    }

    #[test]
    fn test_conti_frac_sqrt() {
        fn test<I, F>(f: F)
        where
            F: Fn(u32) -> I,
            I: Iterator<Item = u32>,
        {
            (0..10)
                .map(|n| (n, n * n)) // perfect squares
                .for_each(|(n, nn)| assert_eq!(f(nn).collect::<Vec<_>>(), [n]));

            for (n, a0, period) in [
                (0, 0, vec![]),
                (1, 1, vec![]),
                (2, 1, vec![2]),
                (3, 1, vec![1, 2]),
                (4, 2, vec![]),
                (5, 2, vec![4]),
                (6, 2, vec![2, 4]),
                (7, 2, vec![1, 1, 1, 4]),
                (8, 2, vec![1, 4]),
                (9, 3, vec![]),
                (10, 3, vec![6]),
                (11, 3, vec![3, 6]),
                (12, 3, vec![2, 6]),
                (13, 3, vec![1, 1, 1, 1, 6]),
                (14, 3, vec![1, 2, 1, 6]),
                (15, 3, vec![1, 6]),
                (16, 4, vec![]),
                (17, 4, vec![8]),
                (18, 4, vec![4, 8]),
                (19, 4, vec![2, 1, 3, 1, 2, 8]),
                (20, 4, vec![2, 8]),
                (21, 4, vec![1, 1, 2, 1, 1, 8]),
                (22, 4, vec![1, 2, 4, 2, 1, 8]),
                (23, 4, vec![1, 3, 1, 8]),
                (24, 4, vec![1, 8]),
                (25, 5, vec![]),
                (26, 5, vec![10]),
                (27, 5, vec![5, 10]),
                (28, 5, vec![3, 2, 3, 10]),
                (29, 5, vec![2, 1, 1, 2, 10]),
                (30, 5, vec![2, 10]),
                (31, 5, vec![1, 1, 3, 5, 3, 1, 1, 10]),
            ] {
                let actual = f(n);
                let expected = [a0].into_iter().chain(period.into_iter().cycle());
                assert_eq!(
                    actual.take(100).collect::<Vec<_>>(),
                    expected.take(100).collect::<Vec<_>>(),
                );
            }
        }

        test(successors::conti_frac_sqrt);
        // test(iterate::conti_frac_sqrt);
    }

    #[allow(clippy::type_complexity)]
    fn cycle_detector_tester(detector: fn(&usize, &dyn Fn(&usize) -> usize) -> (usize, usize)) {
        let tests = [
            (3, 4, 0),
            (0, 5, 0),
            (0, 1, 0),
            (5, 1, 0),
            (1, 10, 0),
            (10, 2, 0),
            (3, 4, 10),
            (0, 3, 5),
            (0, 1, 10),
            (10, 1, 10),
            (0, 100, 0),
            (1, 100, 5),
            (100, 5, 0),
        ];

        for (mu, lambda, x0) in tests {
            let f = move |x: &usize| {
                x.checked_sub(x0 + mu) // entered the cycle?
                    .map_or_else(|| x + 1, |over| x0 + mu + (over + 1) % lambda)
            };

            let (detected_mu, detected_lambda) = detector(&x0, &f);
            assert_eq!(mu, detected_mu);
            assert_eq!(lambda, detected_lambda);
        }
    }

    #[test]
    fn test_detect_cycle_floyd() {
        cycle_detector_tester(|x0, f| successors::detect_cycle_floyd(x0, &f));
        cycle_detector_tester(|x0, f| iterate::detect_cycle_floyd(x0, &f));
    }

    #[test]
    fn test_detect_cycle_brent() {
        cycle_detector_tester(|x0, f| successors::detect_cycle_brent(x0, &f));
        cycle_detector_tester(|x0, f| iterate::detect_cycle_brent(x0, &f));
    }
}
