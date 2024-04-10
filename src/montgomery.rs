use num_bigint::BigUint;

// Elliptic curve in Montgomery form
//     See https://eprint.iacr.org/2017/212.pdf for details
// Curves are defined over the finite field with P elements
// TODO: implement Montgomery arithmetic for BigUint?
pub struct MontgomeryCurve {
    pub a: BigUint,
}

pub struct Point {
    pub x: BigUint,
    pub z: BigUint,
}

impl MontgomeryCurve {
    pub fn new(a: BigUint) -> MontgomeryCurve {
        MontgomeryCurve { a }
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    pub fn add3(&self, p: Point, q: Point, pq: Point) -> Point {
        let v0 = &p.x + &p.z;
        let v1 = &q.x - &q.z;
        let v1 = v1 * v0;

        let v0 = &p.x - &p.z;
        let v2 = &q.x + &q.z;
        let v2 = v2 * v0;

        let v3 = &v1 + &v2;
        let v3 = v3.pow(2);

        let v4 = &v1 - &v2;
        let v4 = v4.pow(2);

        return Point {
            x: pq.z * v3,
            z: pq.x * v4,
        }
    }

    // Given P not equal to origin or infinity, computes P+P
    // Algorithm 2 in Costello & Smith
    pub fn double(&self, p: Point) -> Point {
        let v1 = &p.x + &p.z;
        let v1 = v1.pow(2);

        let v2 = &p.x - &p.z;
        let v2 = v2.pow(2);

        let x = &v1 * &v2;

        let v1 = &v1 - &v2;
        let v3 = (&self.a + 2u32)/4u32;
        let v3 = &v3 - &v2;

        return Point {
            x: x,
            z: v1 * v3,
        }
    }

    
}
