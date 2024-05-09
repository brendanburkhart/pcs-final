use crypto_bigint::{Encoding, Random, U320, U512};
use rand::{thread_rng, Rng};
use rug::float::Round;
use rug::integer::Order;
use rug::Integer;

use crate::constants::{
    ModClassGroup, ModP, BASIS, NUM_PRIMES, ORTHO_NORMS, POOL, PRIMES, PRIMES16,
};
use crate::lattice::{add_slice, dlw_reduce, dot, l1};
use crate::montgomery::{MontgomeryCurve, Point};

#[derive(Debug, Clone, PartialEq)]
pub struct ClassGroupElement {
    pub value: ModClassGroup,
}

impl ClassGroupElement {
    pub fn new(v: ModClassGroup) -> ClassGroupElement {
        ClassGroupElement { value: v }
    }

    pub fn sample() -> ClassGroupElement {
        ClassGroupElement {
            value: ModClassGroup::random(&mut thread_rng()),
        }
    }

    pub fn reduce(&self) -> ReducedClassGroupElement {
        let pool_size = 7500usize;
        let mut b: [Integer; NUM_PRIMES] = [Integer::ZERO; 74];
        let p_mp: Integer = Integer::from_digits(
            &self.value.retrieve().as_limbs().map(|x| x.0.to_le()),
            Order::LsfLe,
        );
        b[0].clone_from(&p_mp);
        for basis_idx in (0..NUM_PRIMES).rev() {
            let (c, _) = (dot(&b, basis_idx) / &ORTHO_NORMS[basis_idx])
                .to_integer_round(Round::Nearest)
                .unwrap();
            let slice = &BASIS[basis_idx * NUM_PRIMES..(basis_idx + 1) * NUM_PRIMES];
            for dim_idx in 0..b.len() {
                b[dim_idx] -= &c * slice[dim_idx];
            }
        }
        let mut e_prime = dlw_reduce(b.map(|x| x.to_i8().unwrap()), pool_size);
        let mut best_len = l1(&e_prime);

        for _ in 0..2 {
            let shifted = {
                let ridx = thread_rng().gen_range(0..10000);
                let ridx2 = thread_rng().gen_range(0..10000);
                add_slice(
                    &add_slice(
                        &e_prime,
                        POOL[ridx * NUM_PRIMES..ridx * NUM_PRIMES + NUM_PRIMES]
                            .try_into()
                            .unwrap(),
                    ),
                    POOL[ridx2 * NUM_PRIMES..ridx2 * NUM_PRIMES + NUM_PRIMES]
                        .try_into()
                        .unwrap(),
                )
            };
            let t = dlw_reduce(shifted, pool_size);
            let norm_t = l1(&t);
            if norm_t < best_len {
                best_len = norm_t;
                e_prime = t;
            }
        }
        ReducedClassGroupElement::new(e_prime)
    }

    pub fn reduce_get_exp(&self) -> [i8; 74] {
        let pool_size = 10000usize;
        let mut b: [Integer; NUM_PRIMES] = [Integer::ZERO; 74];
        let p_mp: Integer = Integer::from_digits(
            &self.value.retrieve().as_limbs().map(|x| x.0.to_le()),
            Order::LsfLe,
        );
        b[0].clone_from(&p_mp);
        for basis_idx in (0..NUM_PRIMES).rev() {
            let (c, _) = (dot(&b, basis_idx) / &ORTHO_NORMS[basis_idx])
                .to_integer_round(Round::Nearest)
                .unwrap();
            let slice = &BASIS[basis_idx * NUM_PRIMES..(basis_idx + 1) * NUM_PRIMES];
            for dim_idx in 0..b.len() {
                b[dim_idx] -= &c * slice[dim_idx];
            }
        }
        let mut e_prime = dlw_reduce(b.map(|x| x.to_i8().unwrap()), pool_size);
        let mut best_len = l1(&e_prime);

        for _ in 0..2 {
            let shifted = {
                let ridx = thread_rng().gen_range(0..10000);
                let ridx2 = thread_rng().gen_range(0..10000);
                add_slice(
                    &add_slice(
                        &e_prime,
                        POOL[ridx * NUM_PRIMES..ridx * NUM_PRIMES + NUM_PRIMES]
                            .try_into()
                            .unwrap(),
                    ),
                    POOL[ridx2 * NUM_PRIMES..ridx2 * NUM_PRIMES + NUM_PRIMES]
                        .try_into()
                        .unwrap(),
                )
            };
            let t = dlw_reduce(shifted, pool_size);
            let norm_t = l1(&t);
            if norm_t < best_len {
                best_len = norm_t;
                e_prime = t;
            }
        }
        e_prime
    }

    pub fn to_bytes(&self) -> [u8; 40] {
        self.value.retrieve().to_be_bytes()
    }

    pub fn from_be_slice(a: &[u8]) -> ClassGroupElement {
        ClassGroupElement {
            value: ModClassGroup::new(&U320::from_be_slice(a)),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReducedClassGroupElement {
    // exponent for each ideal, separated into positive/negative exponents
    // exponents[0][i] is positive exponent for ideal i
    exponents: [[u8; NUM_PRIMES]; 2],
    // (p+1)/k, where k is product of primes which have non-zero exponent
    // in the corresponding positive/negative_exponents
    kernel_order_factor: [ModP; 2],
}

impl ReducedClassGroupElement {
    pub fn new(exponents: [i8; NUM_PRIMES]) -> ReducedClassGroupElement {
        let mut positive_exponents = [0u8; NUM_PRIMES];
        let mut negative_exponents = [0u8; NUM_PRIMES];
        const FOUR: U512 = U512::from_u16(4u16);
        let mut positive_order = ModP::new(&FOUR);
        let mut negative_order = ModP::new(&FOUR);

        for idx in 0..exponents.len() {
            assert!(idx < PRIMES.len());
            assert!(idx < positive_exponents.len());
            assert!(idx < negative_exponents.len());

            let prime = PRIMES[idx];
            match exponents[idx] {
                d if d > 0 => {
                    positive_exponents[idx] = exponents[idx] as u8;
                    negative_exponents[idx] = 0u8;
                    negative_order *= prime;
                }
                d if d < 0 => {
                    positive_exponents[idx] = 0u8;
                    negative_exponents[idx] = (-exponents[idx]) as u8;
                    positive_order *= prime;
                }
                _ => {
                    positive_exponents[idx] = 0u8;
                    negative_exponents[idx] = 0u8;
                    positive_order *= prime;
                    negative_order *= prime;
                }
            }
        }

        ReducedClassGroupElement {
            exponents: [positive_exponents, negative_exponents],
            kernel_order_factor: [positive_order, negative_order],
        }
    }

    //TODO: CONSTANT TIME? SEE https://eprint.iacr.org/2018/1198.pdf.
    pub fn act_on(&self, e: &MontgomeryCurve) -> MontgomeryCurve {
        let mut exponents = self.exponents;
        let mut kernel_order_factor = self.kernel_order_factor;
        let mut done: [bool; 2] = [false, false];

        assert_eq!(e.a.z, ModP::ONE);
        let mut e = MontgomeryCurve::new(e.a.x);

        while !done[0] || !done[1] {
            let p = Point::random();
            //any point is on a curve or its twist
            let sign: usize = 1 - e.on_curve(&p) as usize;
            if done[sign] {
                continue;
            }
            // q = [(p+1)/k]p
            let mut q = e.mult(&p, &kernel_order_factor[sign]);
            done[sign] = true;
            for idx in (0..exponents[sign].len()).rev() {
                assert!(idx < PRIMES.len());
                assert!(idx < exponents[sign].len());
                if exponents[sign][idx] == 0 {
                    continue;
                }
                // compute k/ell
                let mut kl = ModP::ONE;
                for j in 0..idx {
                    if exponents[sign][j] != 0 {
                        kl *= PRIMES[j];
                    }
                }
                // k = [k/ell]q = [(p+1)/ell]p
                let k = e.mult(&q, &kl);
                // need kernel to have order ell, skip if not
                if bool::from(k.is_zero()) {
                    done[sign] &= exponents[sign][idx] == 0;
                    continue;
                }
                // compute isogeny corresponding to action of one ideal
                (q, e) = e.isogeny(&k, PRIMES16[idx] as usize, &q);
                // reduce exponent for that ideal, if it is zero we add it to common kernel order
                exponents[sign][idx] -= 1;
                if exponents[sign][idx] == 0 {
                    kernel_order_factor[sign] *= PRIMES[idx];
                } else {
                    done[sign] = false;
                }
            }

            e = e.normalize();
        }

        return e;
    }
}

// Test ReducedCLassGroupElement
#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::BASE_CURVE;
    
    #[test]
    fn act_on() {
        let e = MontgomeryCurve::new(ModP::ZERO);
        let exp: [i8; NUM_PRIMES] = [
            -5, 2, 0, -3, 4, -4, -5, 3, 5, -1, -2, -4, 0, -2, -3, 3, 1, -2, 5, 3, 4, 3, -4, 2, 2,
            3, -1, 0, 1, -3, 0, 1, -5, -2, 0, 2, 0, 0, -5, 5, 4, 5, 0, -5, 0, -1, 0, 1, 5, 1, 1,
            -3, 0, 5, 1, 2, -1, 1, -5, 0, 1, 5, 3, 2, -1, -5, 4, 2, 1, 2, -2, 0, 1, 5,
        ];

        let g1 = ReducedClassGroupElement::new(exp);
        let b = g1.act_on(&e);

        let correct1 = ModP::new(&U512::from_be_hex("2D3F42F31F984ACE1F45E62D35F7C9936BA51863A204A7AF9562DF7822E01323EAECAB2D86BBA42CB9B1DAA7DAA565800BD5BF35A0297218E8CBDB0399618180"));

        assert_eq!(b.normalize().a.x, correct1);

        let exp: [i8; NUM_PRIMES] = [
            1, -2, 5, 1, 2, 4, -1, 0, -2, -1, 2, 5, -3, 3, 3, -1, -2, -1, 0, -5, -1, -1, -5, 4, 2,
            -1, -1, -5, -4, -3, 4, 1, 4, -2, 4, -5, 3, -1, 1, 2, 0, 4, 1, -5, 4, 1, 4, -1, 0, -5,
            3, -2, -3, 0, -1, 4, 3, -2, -5, -5, 4, 3, 2, 1, -2, 3, 3, -2, -3, -5, 5, 3, -5, 2,
        ];

        let g2 = ReducedClassGroupElement::new(exp);
        let b = g2.act_on(&e);

        let correct2 = ModP::new(&U512::from_be_hex("09EB001955B4E84ECFFE86806E0C8313800D0475CFF3519FAF30DC5F3A060E97AE258051DABED0245406DF3BD41B4A03F3C7756C2DE8DE4AD28AC8CD8D506695"));

        assert_eq!(b.normalize().a.x, correct2);

        let e1 = MontgomeryCurve::new(correct1);
        let c = g2.act_on(&e1);

        let correct3 = ModP::new(&U512::from_be_hex("2BA3EBCD76B29349F525D3B73BA841065926870C3A1F23902EF53652D880BCF6E8D2705B2F94E23551BBFE9F4FD9A4DA1EADF24EA62DC2A7F425A8EB901E31A6"));

        assert_eq!(c.normalize().a.x, correct3);

        let e2 = MontgomeryCurve::new(correct2);
        let c = g1.act_on(&e2);
        assert_eq!(c.normalize().a.x, correct3);

        let a = ModP::new(&U512::from_be_hex("5EB2AEEF49060ED93CC067CC83EDDA45D2494F1CF0EB19F41DA034D00A61CBFDFE7C05C0E2730E14EE51B1C0DD5F10CD4958FB9567E9125410860FADDE6D5306"));
        let e3 = MontgomeryCurve::new(a);
        let exp: [i8; NUM_PRIMES] = [
            -5, 2, -5, -1, -4, -3, 5, 4, -2, 5, 3, -4, 4, 4, 5, 5, -5, -1, -2, -1, 2, -3, 1, -5,
            -2, 5, 5, 5, -2, 2, 3, 4, 2, -5, 4, 2, 1, 4, -3, 1, -3, 0, 5, 4, -4, 0, 0, -3, 1, 3, 0,
            -1, -4, -4, -5, -4, -5, 3, -3, 0, -4, 2, -1, 4, 5, 0, -3, 3, -4, -1, -2, -2, 2, -5,
        ];

        let g = ReducedClassGroupElement::new(exp);
        let b = g.act_on(&e3);

        let correct = ModP::new(&U512::from_be_hex("244BFECEC58AB059E806D5E001BFA1230F5FD3735C2D78EA8F901E4C3FDE881D2FBE39781C948436C538EFA0E2C54B650E390D2B519BD6A3BA6026AFCC819DCB"));
        assert_eq!(correct, b.normalize().a.x);

        let b = g.act_on(&e);
        let correct = ModP::new(&U512::from_be_hex("0313B6847C6679D3E73A9DD53E2C48E7E1279BE4749D519B2CC13FF5F7D8B235944A1994761C0DFD8306A899567D1DE98ECE0F2431C907EAC61CD5E1F34E0E9E"));
        assert_eq!(correct, b.normalize().a.x);
    }

    #[test]
    fn reduce_basic() {
        let e = MontgomeryCurve::new(ModP::ZERO);
        let action1 =
            ClassGroupElement::new(ModClassGroup::new(&U320::from(8792348923892u64))).reduce();
        let e = action1.act_on(&e);
        let action2 =
            ClassGroupElement::new(ModClassGroup::new(&U320::from(8438348348348u64))).reduce();
        let res1 = action2.act_on(&e);

        let e = MontgomeryCurve::new(ModP::ZERO);
        let action3 =
            ClassGroupElement::new(ModClassGroup::new(&U320::from(17230697272240u64))).reduce();
        let res2 = action3.act_on(&e);
        assert_eq!(res1, res2);
    }

    #[test]
    fn reduce() {
        let test = ModClassGroup::new(&U320::from_be_hex(
            "0000000000000000C4FD5F914373FF3AAA5B084A20DC5230D0A53604E44E0C01CF85E801AEE7BD71",
        ));
        let reduced = ClassGroupElement::new(test).reduce_get_exp();
        let correct = [
            -4, 1, 0, 2, 5, 0, 1, -5, 4, -1, 3, 2, 8, -1, 3, -4, 0, -13, -7, 5, 2, -1, 0, 1, 3, -4,
            11, 2, -7, -5, 0, 10, -2, 3, 3, -2, 5, -8, 0, 7, 5, 0, -1, -1, 3, -4, 0, -4, 0, -1, 4,
            0, 4, 1, 2, -2, 4, 1, 0, -4, 0, 1, 2, 8, 0, -3, 0, 0, -3, 0, 0, -2, -1, 3,
        ];
        println!("{} {}", l1(&reduced), l1(&correct));
        assert_eq!(ReducedClassGroupElement::new(reduced).act_on(&BASE_CURVE), ReducedClassGroupElement::new(correct).act_on(&BASE_CURVE));
        let test = ModClassGroup::new(&U320::from_be_hex(
            "0000000000000000D1600A146193CEC7FDC14F739B7D48579B84E1FC56E4884438DBD16D723E105A",
        ));
        let reduced = ClassGroupElement::new(test).reduce_get_exp();
        let correct = [
            5, 1, -1, -1, 0, 4, -2, 1, 1, 3, 0, 0, -11, 1, -2, 2, 0, -1, -3, 0, -10, 0, -5, -2, -6,
            1, 1, 6, 1, 0, -7, 3, 2, 3, 7, -3, -12, -3, 1, 0, -2, 2, -4, 2, 3, 8, -6, -2, 2, 4, 0,
            1, -4, 4, 6, 6, -1, -1, 8, 0, -1, 3, -1, 3, 8, 2, 7, 2, 0, -1, -5, 3, -2, -3,
        ];
        println!("{} {}", l1(&reduced), l1(&correct));
        assert_eq!(ReducedClassGroupElement::new(reduced).act_on(&BASE_CURVE), ReducedClassGroupElement::new(correct).act_on(&BASE_CURVE));
        
    }
}
