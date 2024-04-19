use num_bigint::{BigUint, RandBigInt};
use num_modular::{ModularCoreOps, ModularPow};
use num_traits::{One,Zero};
use std::fmt;

use crate::{classgroup::{HASSE_INTERVAL, INV4, P, PRIMES}, modular};

// Elliptic curve in Montgomery form
//     See Costello & Smith, https://eprint.iacr.org/2017/212.pdf for details
//     Also see CSIDH paper, Castryck et al., https://csidh.isogeny.org/csidh-20181118.pdf
// Curves are defined over the finite field with P elements
// TODO: implement Montgomery arithmetic for BigUint?
pub struct MontgomeryCurve {
    pub a: BigUint,
    ap2d4: BigUint, // (a+2)/4
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: BigUint,
    pub z: BigUint,
}

impl Point {
    pub fn from_x(x: BigUint) -> Point {
        return Point {
            x, z: One::one()
        };
    }

    pub fn is_zero(&self) -> bool {
        return self.z.is_zero();
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.z.is_zero() {
            return write!(f, "1 : 0");
        } else {
            let z_inv = modular::inverse(self.z.clone(), &P);
            let x = (&self.x).mulm(&z_inv, &P);
            return write!(f, "{} : 1", x);
        }
    }
}

impl MontgomeryCurve {
    pub fn new(a: BigUint) -> MontgomeryCurve {
        let a = BigUint::addm(a, BigUint::zero(), &P); // ensure A is reduced mod p
        let two: BigUint = BigUint::from(2u32);
        let ap2d4 = ((&a).addm(&two, &P))
                         .mulm(&(*INV4), &P);

        MontgomeryCurve { a, ap2d4: ap2d4.clone() }
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    fn add3(&self, p: &Point, q: &Point, pq: &Point) -> Point {
        let square = BigUint::from(2u32);

        let v0 = (&p.x).addm(&p.z, &P);
        let v1 = (&q.x).subm(&q.z, &P);
        let v1 = &v1.mulm(&v0, &P);

        let v0 = (&p.x).subm(&p.z, &P);
        let v2 = (&q.x).addm(&q.z, &P);
        let v2 = v2.mulm(&v0, &P);

        let v3 = &v1.addm(&v2, &P);
        let v3 = &v3.powm(&square, &P);

        let v4 =  &v1.subm(&v2, &P);
        let v4 = &v4.powm(&square, &P);

        return Point {
            x: (&pq.z).mulm(v3, &P),
            z: (&pq.x).mulm(v4, &P),
        }
    }

    // Given P not equal to origin or infinity, computes P+P
    // Algorithm 2 in Costello & Smith
    fn double(&self, p: &Point) -> Point {
        let square = BigUint::from(2u32);

        let v1 = (&p.x).addm(&p.z, &P);
        let v1 = v1.powm(&square, &P);

        let v2 = (&p.x).subm(&p.z, &P);
        let v2 = v2.powm(&square, &P);

        let x = (&v1).mulm(&v2, &P);

        let v1 = (&v1).subm(&v2, &P);
        let v3 = (&self.ap2d4).mulm(&v1, &P);
        let v3 = (&v3).addm(&v2, &P);

        return Point {
            x: x.clone(),
            z: (&v1).mulm(&v3, &P),
        }
    }

    // Given P not equal to origin or infinity, computes x([k]P)
    // Algorithm 4 in Costello & Smith
    pub fn mult(&self, p: &Point, k: BigUint) -> Point {
        if k.is_zero() {
            return Point {
                x: One::one(),
                z: Zero::zero(),
            };
        }

        let ell = k.bits();
        let mut x0 = (*p).clone();
        let mut x1 = self.double(p);

        // standard double-and-add, but adapted to x-only arithmetic
        for i in (0..=(ell-2)).rev() {
            let sum = self.add3(&x0, &x1, &p);

            if k.bit(i) {
                x0 = sum;
                x1 = self.double(&x1);
            } else {
                x0 = self.double(&x0);
                x1 = sum;
            }
        }

        return Point {
            x: x0.x,
            z: x0.z
        }
    }

    // Verify curve is nonsingular, i.e. A^2 != 4
    pub fn is_nonsingular(&self) -> bool {
        if self.a == BigUint::from(2u32) {
            return false;
        } else if (&self.a).addm(BigUint::from(2u32), &P) == BigUint::zero() {
            return false
        }

        return true;
    }

    // Algorithm 1 in Castryck et al.
    // Only valid for nonsingular curves
    // TODO: speed-up by sharing computations of (p+1)/ell_i
    pub fn is_supersingular(&self) -> bool {
        // computes (p+1)/(ell_i)
        fn leave_out_one_product(i: usize) -> BigUint {
            let mut k = BigUint::from(4u32);
            for j in 0..PRIMES.len() {
                if j != i {
                    let ell = BigUint::from(PRIMES[j]);
                    k = (&k).mulm(ell, &P);
                }
            }

            return k;
        }

        print!("{}\n\n", leave_out_one_product(0));

        loop {
            let x = rand::thread_rng().gen_biguint_below(&P);
            let p = Point::from_x(x);

            let mut order_divisor = BigUint::one();
            for i in 0..PRIMES.len() {
                let k = leave_out_one_product(i);

                let q = self.mult(&p, k); // [(p+1)/ell_i] * k
                let r = self.mult(&q, BigUint::from(PRIMES[i])); // ell_i * q = [p+1] * k
                if !r.is_zero() {
                    print!("{}, {}, {}, {}, {}\n", BigUint::from(PRIMES[i]), order_divisor, p, q, r);
                    return false; // ell_i divides p+1 but not group order, so cannot be supersingular
                }

                // p+1 is multiple of group order but (p+1)/ell_i is not => ell_i divides group order
                if !q.is_zero() {
                    order_divisor = order_divisor.mulm(BigUint::from(PRIMES[i]), &P);
                }

                // once the Hasse interval [(p+1) - 2*sqrt(p), (p+1) + 2*sqrt(2)] contains only one multiple
                // of the order_divisor, and order_divisor divides p+1, we can be sure the group order is p+1
                if order_divisor > *HASSE_INTERVAL {
                    return true;
                }
            }

            // failed to confirm or refute supersingularity
            // try again with another random point
        }
    }
}
