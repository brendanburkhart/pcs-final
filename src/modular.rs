use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};

pub fn inverse(a: BigUint, n: &BigUint) -> BigUint {
    let mut t0: BigInt = BigInt::from(0i32);
    let mut r0: BigInt = BigInt::from_biguint(Sign::Plus, n.clone());

    let mut t1: BigInt = BigInt::from(1i32);
    let mut r1: BigInt = BigInt::from_biguint(Sign::Plus, a.clone());

    let zero: BigInt = Zero::zero();
    let one: BigInt = One::one();

    while r1 != zero {
        let quotient = &r0 / &r1; // non-modular division
        (t0, t1) = (t1.clone(), &t0 - (&quotient * &t1));
        (r0, r1) = (r1.clone(), &r0 - (&quotient * &r1));
    }

    assert_eq!(r0, one); // if gcd is not one, we cannot invert

    // Euclidean algorithm terminates with |t0| < n, so either 0 <= t0 < n or 0 <= (t0 + n) < n
    if t0.sign() == Sign::Minus {
        t0 = t0 + BigInt::from_biguint(Sign::Plus, n.clone());
        assert_ne!(t0.sign(), Sign::Minus);
    }

    return t0.to_biguint().expect("must be positive");
}
