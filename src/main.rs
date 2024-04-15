use classgroup::PRIMES;
use num_bigint::BigUint;
use num_modular::ModularCoreOps;

use crate::classgroup::{INV4, P};

mod classgroup;
mod montgomery;

fn main() {
    let mut prod: BigUint = BigUint::from(1u32);

    for prime in PRIMES {
        prod = prod * BigUint::from(prime);
    }

    let result = BigUint::from(4u32) * prod - BigUint::from(1u32);
    assert_eq!(result, *classgroup::P);

    print!("{}\n", *INV4);
    assert_eq!(BigUint::from(1u32), (&*INV4).mulm(BigUint::from(4u32), &P));
    print!("{}\n", (&*INV4).mulm(BigUint::from(4u32), &P));
}
