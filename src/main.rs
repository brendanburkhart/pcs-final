use classgroup::PRIMES;
use num_bigint::BigUint;
use num_modular::ModularCoreOps;

mod classgroup;
mod montgomery;

fn main() {
    let mut prod: BigUint = BigUint::from(1u32);

    for prime in PRIMES {
        prod = prod * BigUint::from(prime);
    }

    let result = BigUint::from(4u32) * prod - BigUint::from(1u32);
    assert_eq!(result, *classgroup::P);

    let a = BigUint::from(3u32);
    let b = BigUint::from(5u32);
    let c = BigUint::from(2u32);

    let modulus = BigUint::from(7u32);

    assert_eq!(a, b.mulm(c, &modulus));
}
