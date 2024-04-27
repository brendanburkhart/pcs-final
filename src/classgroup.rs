use lazy_static::lazy_static;
use num_bigint::{BigUint, RandBigInt};
use num_modular::ModularCoreOps;
use num_traits::{One, Zero};

use crate::montgomery::{MontgomeryCurve, Point};

// class group structure for the CSIDH-512 class group,
// data from CSI-FiSH: https://eprint.iacr.org/2019/498

pub const PRIMES: [u16; 74] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
];

lazy_static! {
    // p = 4 * prod(PRIMES) - 1
    pub static ref P: BigUint = (BigUint::from(0x1b81b90533c6c87bu64) * BigUint::from(2u32).pow(0 * 64)) +
    (BigUint::from(0xc2721bf457aca835u64) * BigUint::from(2u32).pow(1 * 64)) +
    (BigUint::from(0x516730cc1f0b4f25u64) * BigUint::from(2u32).pow(2 * 64)) +
    (BigUint::from(0xa7aac6c567f35507u64) * BigUint::from(2u32).pow(3 * 64)) +
    (BigUint::from(0x5afbfcc69322c9cdu64) * BigUint::from(2u32).pow(4 * 64)) +
    (BigUint::from(0xb42d083aedc88c42u64) * BigUint::from(2u32).pow(5 * 64)) +
    (BigUint::from(0xfc8ab0d15e3e4c4au64) * BigUint::from(2u32).pow(6 * 64)) +
    (BigUint::from(0x65b48e8f740f89bfu64) * BigUint::from(2u32).pow(7 * 64));

    pub static ref BASE_CURVE: MontgomeryCurve = MontgomeryCurve::new(BigUint::zero());

    // length of Hasse-interval = 2 * 2*sqrt(p)
    pub static ref HASSE_INTERVAL: BigUint = BigUint::from(0x17895e71e1a20b3fu64) * BigUint::from(2u32).pow(0 * 64) +
                                             BigUint::from(0x38d0cd95f8636a56u64) * BigUint::from(2u32).pow(1 * 64) +
                                             BigUint::from(0x142b9541e59682cdu64) * BigUint::from(2u32).pow(2 * 64) +
                                             BigUint::from(0x856f1399d91d6592u64) * BigUint::from(2u32).pow(3 * 64) +
                                             BigUint::from(0x02u64) * BigUint::from(2u32).pow(4 * 64);
}

#[derive(Debug)]
pub struct ClassGroupElement {
    pub value: BigUint,
}

impl ClassGroupElement {
    pub fn sample() -> ClassGroupElement {
        return ClassGroupElement {
            value: rand::thread_rng().gen_biguint_below(&P),
        };
    }

    // TODO
    pub fn reduce(&self) -> ReducedClassGroupElement {
        return ReducedClassGroupElement::new([0i8; PRIMES.len()]);
    }
}

#[derive(Debug)]
pub struct ReducedClassGroupElement {
    // exponent for each ideal, separated into positive/negative exponents
    // exponents[0][i] is positive exponent for ideal i
    exponents: [[u8; PRIMES.len()]; 2],

    // (p+1)/k, where k is product of primes which have non-zero exponent
    // in the corresponding positive/negative_exponents
    kernel_order_factor: [BigUint; 2],
}

impl ReducedClassGroupElement {
    pub fn new(exponents: [i8; PRIMES.len()]) -> ReducedClassGroupElement {
        let mut positive_exponents = [0u8; PRIMES.len()];
        let mut negative_exponents = [0u8; PRIMES.len()];

        let mut positive_order = BigUint::from(4u32);
        let mut negative_order = BigUint::from(4u32);

        for idx in 0..PRIMES.len() {
            let prime = BigUint::from(PRIMES[idx]);
            if exponents[idx] > 0 {
                positive_exponents[idx] = exponents[idx] as u8;
                negative_exponents[idx] = 0u8;
                negative_order = negative_order.mulm(&prime, &P);
            } else if exponents[idx] < 0 {
                positive_exponents[idx] = 0u8;
                negative_exponents[idx] = (-exponents[idx]) as u8;
                positive_order = positive_order.mulm(&prime, &P);
            } else {
                positive_exponents[idx] = 0u8;
                negative_exponents[idx] = 0u8;
                positive_order = positive_order.mulm(&prime, &P);
                negative_order = negative_order.mulm(&prime, &P);
            }
        }

        return ReducedClassGroupElement {
            exponents: [positive_exponents, negative_exponents],
            kernel_order_factor: [positive_order, negative_order],
        };
    }

    pub fn act_on(&self, e: &MontgomeryCurve) -> MontgomeryCurve {
        let mut exponents = self.exponents.clone();
        let mut kernel_order_factor = self.kernel_order_factor.clone();
        let mut done: [bool; 2] = [false, false];

        assert!(e.a.z.is_one());
        let mut e = MontgomeryCurve::new(e.a.x.clone());

        while !done[0] || !done[1] {
            let p = Point::random();
            let sign: usize = if e.on_curve(&p) { 0 } else { 1 };

            if done[sign] {
                continue;
            }

            // q = [(p+1)/k]p
            let mut q = e.mult(&p, &kernel_order_factor[sign]);

            done[sign] = true;
            for idx in (0..PRIMES.len()).rev() {
                if exponents[sign][idx] == 0 {
                    continue;
                }

                let ell = PRIMES[idx] as usize;

                // compute k/ell
                let mut kl = BigUint::one();
                for j in 0..idx {
                    if exponents[sign][j] != 0 {
                        let prime = BigUint::from(PRIMES[j]);
                        kl = kl.mulm(prime, &P);
                    }
                }

                // k = [k/ell]q = [(p+1)/ell]p
                let k = e.mult(&q, &kl);
                // need kernel to have order ell, skip if not
                if k.is_zero() {
                    done[sign] &= exponents[sign][idx] == 0;
                    continue;
                }

                // compute isogeny corresponding to action of one ideal
                (q, e) = e.isogeny(&k, ell, &q);

                // reduce exponent for that ideal, if it is zero we add it to common kernel order
                exponents[sign][idx] -= 1;
                if exponents[sign][idx] == 0 {
                    let ell = BigUint::from(ell);
                    kernel_order_factor[sign] = (&kernel_order_factor[sign]).mulm(ell, &P);
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
    use std::str::FromStr;

    use sha3::digest::core_api::CoreProxy;

    use super::*;

    #[test]
    fn act_on() {
        let e = MontgomeryCurve::new(BigUint::zero());

        let exp: [i8; PRIMES.len()] = [-5, 2, 0, -3, 4, -4, -5, 3, 5, -1, -2, -4, 0, -2, -3, 3, 1, -2, 5, 3, 4, 3, -4, 2, 2, 3, -1, 0, 1, -3, 0, 1, -5, -2, 0, 2, 0, 0, -5, 5, 4, 5, 0, -5, 0, -1, 0, 1, 5, 1, 1, -3, 0, 5, 1, 2, -1, 1, -5, 0, 1, 5, 3, 2, -1, -5, 4, 2, 1, 2, -2, 0, 1, 5];

        let g1 = ReducedClassGroupElement::new(exp);
        let b = g1.act_on(&e);

        let correct1 = BigUint::from_str("2369783717237495611979301873784587263452377304093138963788722599679539491451202199649554417336783264170944788809733323781398987848062812298930717909418368").unwrap();

        assert_eq!(b.normalize().a.x, correct1);

        let exp: [i8; PRIMES.len()] = [1, -2, 5, 1, 2, 4, -1, 0, -2, -1, 2, 5, -3, 3, 3, -1, -2, -1, 0, -5, -1, -1, -5, 4, 2, -1, -1, -5, -4, -3, 4, 1, 4, -2, 4, -5, 3, -1, 1, 2, 0, 4, 1, -5, 4, 1, 4, -1, 0, -5, 3, -2, -3, 0, -1, 4, 3, -2, -5, -5, 4, 3, 2, 1, -2, 3, 3, -2, -3, -5, 5, 3, -5, 2];

        let g2 = ReducedClassGroupElement::new(exp);
        let b = g2.act_on(&e);

        let correct2 = BigUint::from_str("519446251179368208408483229789219744648346503395090107025372993121474090347141973291186283710582000417720250680500039355089180430360119430964937843828373").unwrap();

        assert_eq!(b.normalize().a.x, correct2);

        let e1 = MontgomeryCurve::new(correct1);
        let c = g2.act_on(&e1);

        let correct3 = BigUint::from_str("2285628850849164625571393133696472708332361107883389592548326297859453099536834098169003372179749949769887257494210901352194177027170150449257705589584294").unwrap();

        assert_eq!(c.normalize().a.x, correct3);

        let e2 = MontgomeryCurve::new(correct2);
        let c = g1.act_on(&e2);
        assert_eq!(c.normalize().a.x, correct3);

        let a = BigUint::from_str("4959735746944445429437167063372787866221260859298686717034957233168920087384110868789056610865328902674826851391106090304366031309421559147593361324331782").unwrap();
        let e3 = MontgomeryCurve::new(a);
        let exp: [i8; PRIMES.len()] = [
            -5, 2, -5, -1, -4, -3, 5, 4, -2, 5, 3, -4, 4, 4, 5, 5, -5, -1, -2, -1, 2, -3, 1, -5, -2, 5, 5, 5, -2, 2, 3, 4, 2, -5, 4, 2, 1, 4, -3, 1, -3, 0, 5, 4, -4, 0, 0, -3, 1, 3, 0, -1, -4, -4, -5, -4, -5, 3, -3, 0, -4, 2, -1, 4, 5, 0, -3, 3, -4, -1, -2, -2, 2, -5
       ];
       let g = ReducedClassGroupElement::new(exp);
       let b = g.act_on(&e3);

       let correct = BigUint::from_str("1901020642689517378383562530566654825675635738058812160599579539988231749004280847525732806286957784867696001219607874569066119129809664916689905063665099").unwrap();
       assert_eq!(correct, b.normalize().a.x);

       let b = g.act_on(&e);
       let correct = BigUint::from_str("161155762622134743756219226521144378586403701443189697279645241145897980207434329778060118292207019642201946026095259577408887611036821226090547555208862").unwrap();
       assert_eq!(correct, b.normalize().a.x);
    }
}