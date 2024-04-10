use lazy_static::lazy_static;
use num_bigint::BigUint;

// class group structure for the CSIDH-512 class group,
// data from CSI-FiSH: https://eprint.iacr.org/2019/498

pub const PRIMES: [u16; 74] = [
     3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
    61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
   139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
   229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
   317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
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
}
