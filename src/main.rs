use std::str::FromStr;

use num_bigint::BigUint;

use crate::montgomery::MontgomeryCurve;

mod classgroup;
mod modular;
mod montgomery;

fn main() {
    let m0 = MontgomeryCurve::new(BigUint::from(0u32));
    let m1 = MontgomeryCurve::new(BigUint::from(2u32));
    let minus_two: BigUint = (&*classgroup::P) - BigUint::from(2u32);
    let m2 = MontgomeryCurve::new(minus_two);
    let big = BigUint::from_str("3761345407298064040496734252078593163266383672948683311982036640586525413088935630843747460636021803468026326856059472255583832635669585412131689148241997").unwrap();
    let m3 = MontgomeryCurve::new(big);
    print!("m0: {}, {}\n", m0.is_nonsingular(), m0.is_supersingular());
    print!("m1: {}, {}\n", m1.is_nonsingular(), m1.is_supersingular());
    print!("m2: {}, {}\n", m2.is_nonsingular(), m2.is_supersingular());
    print!("m3: {}, {}\n", m3.is_nonsingular(), m3.is_supersingular());
}
