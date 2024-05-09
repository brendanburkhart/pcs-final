use std::mem;
use std::ops::{RemAssign, ShrAssign};

use crypto_bigint::{Integer, NonZero, U512, Zero};
use crypto_bigint::modular::constant_mod::ResidueParams;

use crate::constants::{ModP, ModulusP};

//TODO: Bernstein-Yang?
#[allow(dead_code)]
pub fn jacobi_vartime(pt: &ModP) -> i8 {
    let mut t = 1;
    // (a / P) = (aR / P)/(R / P)
    // (R / P) = 1 ==> (a / P) = (aR / P)
    // so we don't need to convert from montgomery form
    let mut a = *pt.as_montgomery();
    let mut n = ModulusP::MODULUS;
    while !bool::from(a.is_zero()) {
        while bool::from(a.is_even()) {
            a.shr_assign(1);
            let r = n.as_limbs()[0].0 & 7u64;
            if (r == 3) || (r == 5) {
                t *= -1;
            }
        }
        mem::swap(&mut a, &mut n);
        if ((a.as_limbs()[0].0 & 3u64) == 3) && ((n.as_limbs()[0].0 & 3u64) == 3) {
            t *= -1;
        }
        a.rem_assign(NonZero::from_uint(n));
    }
    return if n == U512::ONE { t } else { 0 };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jacobi() {
        assert_eq!(jacobi_vartime(&(ModP::new(&U512::ZERO))), 0);

        let ex1 = U512::from_be_hex("1454b04302ca6d736d14ef4731e4fa4523575b01b70fa655bb3689add5be0940adfae6b0876c77cde9bcbca7395f8bc5bf4978515358e0d9a4b0c912585be75b");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex1))), -1);

        let ex2 = U512::from_be_hex("04df799937699b8fd86229aa453e293fec5b358d4069d0f55a56ed18f83ae4b1e54614414731b8ac35b3c21e1ae3ed0eff89941ad4f16117868cfc7bca2a2dcd");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex2))), 1);

        let ex3 = U512::from_be_hex("480a8f23512bad8926e1b36a5e03407667bde3ecd8351cf849bdfa1be30b3dc936c52b668da76e86471ef958e5bb50252ce4bf3603a399dd369bcb9064f059c2");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex3))), -1);

        let ex4 = U512::from_be_hex("5af6b23c068158872f03883e10b3026fef2fd5708755c493cee605c7e6c3243718135bac57505704e59d241b71b10e0ad58bff17bb42adfdb40d428bcf2febe0");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex4))), 1);

        let ex5 = U512::from_be_hex("0ccd53b639e77335d0bbee7bc9d9a7eb1ecdcb0b7132c5dcb9655277207ee4b4370a5c5aa459c52811243ecf8e5e6fab2a7d5fbcfe239e4a7902f0355ff4c893");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex5))), 1);

        let ex6 = U512::from_be_hex("3ee6fa6705cfa82647462b61857472f4ee1ab628eac244e35c3e5be9e82702db1030d0f4cf6f9dfbd05ac0eeffa43e12c830593a2ea5297247c1c370d6e1a1ee");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex6))), -1);

        let ex7 = U512::from_be_hex("5d669e71a6c70a600ac4ae886d52d6ee736d7658e7f88a2400649a47e129953de7b626786c4bd019bf2d4d5a267c7e75418ccd0f681ff25bbb4c44e76ff2a3c5");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex7))), 1);

        let ex8 = U512::from_be_hex("20e773f68bb987735d4f22b593f787c60b0ccd7b4e9af3f96690d2beef34f913021af413846cee6c90548dad4f82938d410e7e861c393fe54d889d530add3245");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex8))), -1);

        let ex9 = U512::from_be_hex("4804801fffc3028778a6c739ef17c78788f873d256ff372404942621fdf7440711dff50259a62cd75bd4a7e1afe89a95de178b5bd61bcb10c32c44a6920e4470");
        assert_eq!(jacobi_vartime(&(ModP::new(&ex9))), 1);
    }
}
