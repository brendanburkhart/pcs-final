use classgroup::PRIMES;
use num_bigint::BigUint;

use crate::montgomery::MontgomeryCurve;

mod classgroup;
mod modular;
mod montgomery;

fn main() {
    let mut prod: BigUint = BigUint::from(1u32);

    for prime in PRIMES {
        prod = prod * BigUint::from(prime);
    }

    let result = BigUint::from(4u32) * prod - BigUint::from(1u32);
    assert_eq!(result, *classgroup::P);

    let m = MontgomeryCurve::new(BigUint::from(0u32));
    let p = montgomery::Point {
        x: BigUint::from(89279u32),
        z: BigUint::from(1u32)
    };

    let k = 157u32;

    print!("{:?}\n", p);

    let mut a = p.clone();
    let mut b = m.double(&p);

    for _i in 0..(k-1) {
        let temp = b.clone();
        b = m.add3(&b, &p, &a);
        a = temp;
    }

    print!("{}\n", m.ladder(&p, BigUint::from(k)));
    print!("{}\n", a);

    // print!("{}\n", *INV4);
    // assert_eq!(BigUint::from(1u32), (&*INV4).mulm(BigUint::from(4u32), &P));
    // print!("{}\n", (&*INV4).mulm(BigUint::from(4u32), &P));
}
