mod classgroup;
mod constants;
mod lattice;
mod modular;
mod montgomery;

mod signatures;

mod merkle;

use pcs_final::signatures::{KeyPair, Signature, SigningKey, VerifyingKey};
use std::env;
use std::fs;
use std::fs::File;
use std::io::Write;

const CURVES: usize = 256;
const ROUNDS: usize = 13;

fn main() {
    // Get command-line arguments
    let args: Vec<String> = env::args().collect();

    // Check if the correct number of arguments are provided
    if args.len() < 3 {
        println!(
            "Usage: {} --keygen <public key file> <secret key file> \nUsage: {} --sign <file to sign> <secret key> \nUsage: {} --verify <file to verify> <signature file> <public key file>",
            args[0], args[0], args[0]
        );
        return;
    }

    // Check the flag and validate arguments accordingly
    let flag = &args[1];
    match flag.as_str() {
        "--keygen" => {
            if args.len() != 4 {
                println!(
                    "Usage: {} --keygen <public key file> <secret key file>",
                    args[0]
                );
                return;
            }
            let pk_path = &args[2];
            let sk_path = &args[3];
            let keys = KeyPair::<CURVES, ROUNDS>::generate();
            let mut pk_file = File::create(pk_path).expect("Couldn't create public key file");
            let mut sk_file = File::create(sk_path).expect("Couldn't create secret key file");
            pk_file
                .write_all(&keys.public_key().serialize())
                .expect("Couldn't write public key file");
            sk_file
                .write_all(&keys.serialized_secret_key())
                .expect("Couldn't write secret key file");
            println!("Wrote public key to {pk_path}");
            println!("Wrote secret key to {sk_path}");
        }
        "--sign" => {
            if args.len() != 4 {
                println!("Usage: {} --sign <file to sign> <secret key> ", args[0]);
                return;
            }
            let mut file_path = &args[2];
            if !validate_path(file_path) {
                println!("Invalid <file to sign>");
            }
            let sk_path = &args[3];
            if !validate_path(sk_path) {
                println!("Invalid <secret key>");
            }
            let file_bytes = fs::read(file_path).expect("Failed to read <file path>");
            let sk =
                SigningKey::deserialize(fs::read(sk_path).expect("Failed to read <secret key>"))
                    .expect("Corrupted <secret key>");
            let sig_path = format!("{}.sig", file_path);
            let mut sig_file = File::create(&sig_path).expect("Couldn't create signature file");
            sig_file
                .write_all(
                    &KeyPair::<CURVES, ROUNDS>::from_signing_key(sk)
                        .sign(&file_bytes)
                        .serialize(),
                )
                .expect("Couldn't write signature file");
            println!("Wrote signature to {sig_path}");
        }
        "--verify" => {
            if args.len() != 5 {
                panic!(
                    "Usage: {} --verify <file to verify> <signature file> <public key file>",
                    args[0]
                );
            }
            let verify_path = &args[2];
            let sig_path = &args[3];
            let pk_path = &args[4];
            if !validate_path(verify_path) {
                println!("Invalid <file to verify>");
            }
            if !validate_path(sig_path) {
                println!("Invalid <signature file>");
            }
            if !validate_path(pk_path) {
                println!("Invalid <public key file>");
            }
            let verify_bytes = fs::read(verify_path).expect("Failed to read <file path>");
            let sig_bytes = fs::read(sig_path).expect("Failed to read <file path>");
            let sig = Signature::deserialize(&sig_bytes).expect("Corrupted <signature>");
            let pk_bytes = fs::read(pk_path).expect("Failed to read <file path>");
            let pk = VerifyingKey::deserialize(&pk_bytes).expect("Corrupted <public key>");
            if pk.verify(&sig, &verify_bytes).is_ok() {
                println!("VALID");
                return;
            };
            println!("INVALID");
        }
        _ => {
            println!(
                "Usage: {} --keygen <public key file> <secret key file> \nUsage: {} --sign <file to sign> <secret key> \nUsage: {} --verify <file to verify> <signature file> <public key file>",
                args[0], args[0], args[0]
            );
        }
    }
}

fn validate_path(path: &String) -> bool {
    // Use metadata to check if the path exists and is a file
    match fs::metadata(path) {
        Ok(metadata) => metadata.is_file(),
        Err(_) => false,
    }
}
