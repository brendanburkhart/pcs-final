use criterion::{Criterion, criterion_group, criterion_main};
use rand::{RngCore, thread_rng};

use pcs_final::signatures::KeyPair;

fn bench_keygen_256_13(c: &mut Criterion) {
    c.bench_function("bench_keygen_256_13", |b| {
        b.iter(|| {
            std::hint::black_box({
                KeyPair::<256, 13>::generate()
            });
        });
    });
}

fn bench_sign_256_13(c: &mut Criterion) {
    let keypair = KeyPair::<256, 13>::generate();
    let mut msg = [0u8; 1024];
    thread_rng().fill_bytes(&mut msg);
    c.bench_function("bench_sign_256_13", |b| {
        b.iter(|| {
            std::hint::black_box({
                keypair.sign(&msg)
            });
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
            std::hint::black_box({
                pk.verify(&signature, &msg)
            });
        });
    });
}

criterion_group!(name = benches;
    config = Criterion::default();

    targets =
        bench_keygen_256_13,
        bench_sign_256_13,
        bench_verify_256_13
);
criterion_main!(benches);
