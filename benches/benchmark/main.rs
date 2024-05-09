use std::time::Duration;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, RngCore};

use pcs_final::signatures::KeyPair;

fn bench_keygen_256_13(c: &mut Criterion) {
    c.bench_function("bench_keygen_256_13", |b| {
        b.iter(|| {
            std::hint::black_box({ KeyPair::<256, 13>::generate() });
        });
    });
}

fn bench_sign_256_13(c: &mut Criterion) {
    let keypair = KeyPair::<256, 13>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_256_13", |b| {
        b.iter(|| {
            std::hint::black_box({ keypair.sign(&msg) });
        });
    });
}

fn bench_verify_256_13(c: &mut Criterion) {
    let keypair = KeyPair::<256, 13>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    let signature = keypair.sign(&msg);
    let pk = keypair.public_key();
    c.bench_function("bench_verify_256_13", |b| {
        b.iter(|| {
            std::hint::black_box({ pk.verify(&signature, &msg) });
        });
    });
}

fn bench_keygen_1024_11(c: &mut Criterion) {
    c.bench_function("bench_keygen_1024_11", |b| {
        b.iter(|| {
            std::hint::black_box({ KeyPair::<1024, 11>::generate() });
        });
    });
}

fn bench_sign_1024_11(c: &mut Criterion) {
    let keypair = KeyPair::<1024, 11>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_1024_11", |b| {
        b.iter(|| {
            std::hint::black_box({ keypair.sign(&msg) });
        });
    });
}

fn bench_verify_1024_11(c: &mut Criterion) {
    let keypair = KeyPair::<1024, 11>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    let signature = keypair.sign(&msg);
    let pk = keypair.public_key();
    c.bench_function("bench_verify_1024_11", |b| {
        b.iter(|| {
            std::hint::black_box({ pk.verify(&signature, &msg) });
        });
    });
}

fn bench_keygen_4096_9(c: &mut Criterion) {
    c.bench_function("bench_keygen_4096_9", |b| {
        b.iter(|| {
            std::hint::black_box({ KeyPair::<4096, 9>::generate() });
        });
    });
}

fn bench_sign_4096_9(c: &mut Criterion) {
    let keypair = KeyPair::<4096, 9>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_4096_9", |b| {
        b.iter(|| {
            std::hint::black_box({ keypair.sign(&msg) });
        });
    });
}

fn bench_verify_4096_9(c: &mut Criterion) {
    let keypair = KeyPair::<4096, 9>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    let signature = keypair.sign(&msg);
    let pk = keypair.public_key();
    c.bench_function("bench_verify_4096_9", |b| {
        b.iter(|| {
            std::hint::black_box({ pk.verify(&signature, &msg) });
        });
    });
}

fn bench_keygen_32768_7(c: &mut Criterion) {
    c.bench_function("bench_keygen_32768_7", |b| {
        b.iter(|| {
            std::hint::black_box({ KeyPair::<32768, 7>::generate() });
        });
    });
}

fn bench_sign_32768_7(c: &mut Criterion) {
    let keypair = KeyPair::<32768, 7>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_32768_7", |b| {
        b.iter(|| {
            std::hint::black_box({ keypair.sign(&msg) });
        });
    });
}

fn bench_verify_32768_7(c: &mut Criterion) {
    let keypair = KeyPair::<32768, 7>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    let signature = keypair.sign(&msg);
    let pk = keypair.public_key();
    c.bench_function("bench_verify_32768_7", |b| {
        b.iter(|| {
            std::hint::black_box({ pk.verify(&signature, &msg) });
        });
    });
}

fn bench_keygen_262144_6(c: &mut Criterion) {
    c.bench_function("bench_keygen_262144_6", |b| {
        b.iter(|| {
            std::hint::black_box({ KeyPair::<262144, 6>::generate() });
        });
    });
}

fn bench_sign_262144_6(c: &mut Criterion) {
    let keypair = KeyPair::<262144, 6>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_262144_6", |b| {
        b.iter(|| {
            std::hint::black_box({ keypair.sign(&msg) });
        });
    });
}

fn bench_verify_262144_6(c: &mut Criterion) {
    let keypair = KeyPair::<262144, 6>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    let signature = keypair.sign(&msg);
    let pk = keypair.public_key();
    c.bench_function("bench_verify_262144_6", |b| {
        b.iter(|| {
            std::hint::black_box({ pk.verify(&signature, &msg) });
        });
    });
}


criterion_group!(name = benches;
    config = Criterion::default();//.sample_size(50);
    targets =
        bench_sign_256_13,
        bench_verify_256_13,
        bench_sign_1024_11,
        bench_verify_1024_11,
        bench_sign_4096_9,
        bench_verify_4096_9,
        bench_sign_32768_7,
        bench_verify_32768_7,
        bench_sign_262144_6,
        bench_verify_262144_6
);
criterion_main!(benches);
