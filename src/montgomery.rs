use num_bigint::BigUint;
use num_modular::{ModularCoreOps, ModularPow};
use num_traits::{One,Zero};
use std::fmt;

use crate::{classgroup::{INV4, P}, modular};

// Elliptic curve in Montgomery form
//     See https://eprint.iacr.org/2017/212.pdf for details
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
        let two: BigUint = BigUint::from(2u32);
        let ap2d4 = ((&a).addm(&two, &P))
                         .mulm(&(*INV4), &P);

        MontgomeryCurve { a, ap2d4: ap2d4.clone() }
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    pub fn add3(&self, p: &Point, q: &Point, pq: &Point) -> Point {
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
    pub fn double(&self, p: &Point) -> Point {
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
    pub fn ladder(&self, p: &Point, k: BigUint) -> Point {
        if k.is_zero() {
            return Point {
                x: One::one(),
                z: Zero::zero(),
            };
        }

        let ell = k.bits();
        let mut x0 = (*p).clone();
        let mut x1 = self.double(p);
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
}
