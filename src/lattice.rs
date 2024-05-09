use rug::{Float, Integer};

use crate::constants::{FLOAT_PRECISION, NUM_PRIMES, ORTHO_BASIS, POOL};

pub fn dot(b: &[Integer; 74], basis_idx: usize) -> Float {
    let slice = &ORTHO_BASIS[NUM_PRIMES * basis_idx..NUM_PRIMES * (basis_idx + 1)];
    slice
        .iter()
        .zip(b)
        .fold(Float::with_val(FLOAT_PRECISION, 0.0), |accum, (s, b)| {
            accum + Float::with_val(FLOAT_PRECISION, s * b)
        })
}

pub fn sub_slice(a: &[i8], b: &[i8]) -> [i8; 74] {
    let mut result = [0i8; 74];
    for idx in 0..result.len() {
        assert!(idx < a.len());
        assert!(idx < b.len());
        result[idx] = a[idx] - b[idx];
    }
    result
}

pub fn add_slice(a: &[i8], b: &[i8]) -> [i8; 74] {
    let mut result = [0i8; 74];
    for idx in 0..result.len() {
        assert!(idx < a.len());
        assert!(idx < b.len());
        result[idx] = a[idx] + b[idx];
    }
    result
}

pub fn l1(a: &[i8; 74]) -> u16 {
    a.iter().fold(0u16, |s, x| s + (x.unsigned_abs() as u16))
}

pub fn dlw_reduce(e: [i8; 74], pool_size: usize) -> [i8; 74] {
    let mut e_prime = e;
    let mut stalled = false;
    let mut best_norm = l1(&e_prime);
    while !stalled {
        stalled = true;
        for idx in 0..pool_size {
            let s: &[i8] = &POOL[idx * NUM_PRIMES..idx * NUM_PRIMES + 74];
            let diff = &sub_slice(&e_prime, s);
            let sum = &add_slice(&e_prime, s);
            let l1diff = l1(diff);
            let l1sum = l1(sum);
            if l1sum < best_norm {
                best_norm = l1sum;
                e_prime = *sum;
                stalled = false;
            }
            if l1diff < best_norm {
                best_norm = l1diff;
                e_prime = *diff;
                stalled = false;
            }
        }
    }
    e_prime
}
