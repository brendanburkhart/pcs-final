# Practical Cryptographic Systems Final Project

Rust implementation of the CSIDH and CSI-FiSH post-quantum cryptographic schemes, in partial fulfillment of the final project for EN.601.445 Practical Cryptographic Systems at Johns Hopkins University.

### Building & installing

You will need to build GMP, MPFR, and MPC from source &mdash; this should be relatively fast. In particular, we use
- GMP v6.3.0, this can be downloaded from https://gmplib.org/#DOWNLOAD
- MPFR v4.2.1, which can be downloaded from https://www.mpfr.org/mpfr-current/#download,
- MPC v1.3.1, which can be downloaded from https://www.multiprecision.org/mpc/download.html

You will need also Rust and Cargo (the Rust build system) installed. Clone this repository, and run `cargo build --release` to compile.

Unit/integration tests can be done via `cargo test --release`.

### Use

This project provides both binary and library targets, and so can be used either as a library via the public API, or as a command-line tool to create/verify CSI-FiSH file signatures.

Assuming the standard release target has been built, the command-line interface is:
```
./target/release/pcs_final --keygen <public key file> <secret key file>
./target/release/pcs_final --sign <file to sign> <secret key>
./target/release/pcs_final --verify <file to verify> <signature file> <public key file>
```

The included statistical performance benchmarks can be run via `cargo bench`.
