use std::array::TryFromSliceError;
use std::fmt::Error;
use std::io::Read;

use crypto_bigint::Encoding;
use rayon::prelude::*;
use sha3::digest::{ExtendableOutput, Update};

use crate::classgroup::ClassGroupElement;
use crate::constants::{DeserializationError, HashType, Hasher, VerificationFailed};
use crate::constants::{BASE_CURVE, HASH_SIZE};
use crate::merkle::{ClassGroupMerkleProof, ClassGroupMerkleTree};
use crate::montgomery::MontgomeryCurve;

#[derive(Debug, PartialEq)]
pub struct SigningKey {
    proof_tree: ClassGroupMerkleTree,
    public_curves: Vec<MontgomeryCurve>,
    _secret_actions: Vec<ClassGroupElement>,
}

impl SigningKey {
    pub fn serialize(&self) -> Vec<u8> {
        let num_curves = (self.public_curves.len() as u32).to_be_bytes().to_vec();
        let serialized_curves = self
            .public_curves
            .iter()
            .flat_map(|curve| curve.to_bytes().expect("bad curve"))
            .collect::<Vec<u8>>();

        let num_secret_cge = (self._secret_actions.len() as u32).to_be_bytes().to_vec();
        let serialized_secret_cge = self
            ._secret_actions
            .iter()
            .flat_map(|cge| cge.value.retrieve().to_be_bytes())
            .collect::<Vec<u8>>();

        let serialized_proof_tree = self.proof_tree.serialize();
        let serialized_secret_cge_len = (serialized_secret_cge.len() as u32 + 4)
            .to_be_bytes()
            .to_vec();
        let serialized_curves_len = (serialized_curves.len() as u32 + 4).to_be_bytes().to_vec();
        let serialized_proof_len = (serialized_proof_tree.len() as u32).to_be_bytes().to_vec();
        [
            serialized_secret_cge_len,
            serialized_curves_len,
            serialized_proof_len,
            num_secret_cge,
            serialized_secret_cge,
            num_curves,
            serialized_curves,
            serialized_proof_tree,
        ]
        .concat()
    }

    pub fn deserialize(key: Vec<u8>) -> Result<SigningKey, DeserializationError> {
        let serialized_secret_cge_len = <u32>::from_be_bytes(key[..4].try_into()?) as usize;
        let serialized_curves_len = <u32>::from_be_bytes(key[4..8].try_into()?) as usize;
        let serialized_proof_len = <u32>::from_be_bytes(key[8..12].try_into()?) as usize;
        let sig = &key[12..];
        let serialized_secret_cge_bytes = &sig[..serialized_secret_cge_len];
        let sig = &sig[serialized_secret_cge_len..];
        let serialized_curves_bytes = &sig[..serialized_curves_len];
        let sig = &sig[serialized_curves_len..];
        let serialized_proof_bytes = &sig[..serialized_proof_len];
        if serialized_secret_cge_bytes.len() <= 4 {
            return Err(DeserializationError);
        }
        let num_cge = <u32>::from_be_bytes(serialized_secret_cge_bytes[0..4].try_into()?);
        let num_cge_bytes = serialized_secret_cge_bytes.len() - 4;
        if !(num_cge_bytes % 40 == 0 && ((num_cge_bytes as u32 / 40) == num_cge)) {
            return Err(DeserializationError);
        }
        let deserialized_cge: Vec<ClassGroupElement> = serialized_secret_cge_bytes[4..]
            .chunks_exact(40)
            .map(ClassGroupElement::from_be_slice)
            .collect();
        if serialized_curves_bytes.len() <= 4 {
            return Err(DeserializationError);
        }
        let num_public_curves = <u32>::from_be_bytes(serialized_curves_bytes[0..4].try_into()?);
        let num_curves_bytes = serialized_curves_bytes.len() - 4;
        if !(num_curves_bytes % 64 == 0 && ((num_curves_bytes as u32 / 64) == num_public_curves)) {
            return Err(DeserializationError);
        }
        let deserialized_curves: Vec<MontgomeryCurve> = serialized_curves_bytes[4..]
            .chunks_exact(64)
            .map(MontgomeryCurve::from_be_slice)
            .collect();
        let deserialized_proof_tree = ClassGroupMerkleTree::deserialize(serialized_proof_bytes)?;

        Ok(SigningKey {
            proof_tree: deserialized_proof_tree,
            public_curves: deserialized_curves,
            _secret_actions: deserialized_cge,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VerifyingKey {
    root: HashType,
    merkle: HashType,
}

impl VerifyingKey {
    pub fn serialize(&self) -> Vec<u8> {
        [
            HASH_SIZE.to_be_bytes().to_vec(),
            self.root.to_vec(),
            self.merkle.to_vec(),
        ]
        .concat()
    }
    pub fn deserialize(vk: &[u8]) -> Result<VerifyingKey, DeserializationError> {
        let hash_size = <u32>::from_be_bytes(vk[..4].try_into()?) as usize;
        let root: HashType = vk[4..4 + hash_size].try_into()?;
        let merkle: HashType = vk[4 + hash_size..4 + hash_size + hash_size]
            .try_into()
            .unwrap();
        Ok(VerifyingKey { root, merkle })
    }
}

pub struct KeyPair<const CURVES: usize, const ROUNDS: usize> {
    sk: SigningKey,
    pk: VerifyingKey,
}

impl<const CURVES: usize, const ROUNDS: usize> KeyPair<CURVES, ROUNDS> {
    pub fn generate() -> KeyPair<CURVES, ROUNDS> {
        let (cge, curves): (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) =
            Self::random_curves(CURVES);
        let _tree = ClassGroupMerkleTree::from_leaves(&curves);
        KeyPair {
            pk: VerifyingKey {
                root: _tree.root(),
                merkle: _tree.merkle_key(),
            },
            sk: SigningKey {
                proof_tree: _tree,
                public_curves: curves,
                _secret_actions: cge,
            },
        }
    }

    pub fn from_signing_key(sk: SigningKey) -> KeyPair<CURVES, ROUNDS> {
        KeyPair {
            pk: VerifyingKey {
                root: sk.proof_tree.root(),
                merkle: sk.proof_tree.merkle_key(),
            },
            sk,
        }
    }

    fn random_curves(num_curves: usize) -> (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) {
        (0..num_curves)
            .into_par_iter()
            .map(|_| {
                let r = ClassGroupElement::sample();
                let curve = r.reduce().act_on(&BASE_CURVE);
                (r, curve)
            })
            .unzip()
    }

    pub fn sign(&self, message: &[u8]) -> Signature {
        let (b, ephemeral_curves): (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) =
            Self::random_curves(ROUNDS);
        let mut hasher = Hasher::default();
        for curve in &ephemeral_curves {
            hasher.update(&curve.to_bytes().unwrap());
        }

        hasher.update(message);
        let mut reader = hasher.finalize_xof();
        let mut challenges: Vec<u8> = Vec::new();
        let mut ephemeral_cge: Vec<ClassGroupElement> = Vec::new();
        let mut opened_curve_indices: Vec<usize> = Vec::new();
        let mut opened_curves: Vec<MontgomeryCurve> = Vec::new();

        for ephem_cge in b {
            let mut challenge = [0u8; 4];
            sha3::digest::XofReader::read(&mut reader, &mut challenge);
            challenges.extend_from_slice(&challenge);
            let n = <i32>::from_be_bytes(challenge);
            let curve_num = (n.unsigned_abs() as usize) % CURVES;
            // ri = bi − sign(ci)a|ci| mod ClassNum
            let s = if n > 0 {
                opened_curves.push(self.sk.public_curves[curve_num].clone());
                ephem_cge.value - self.sk._secret_actions[curve_num].value
            } else {
                opened_curves.push(self.sk.public_curves[curve_num].twist());
                ephem_cge.value + self.sk._secret_actions[curve_num].value
            };
            ephemeral_cge.push(ClassGroupElement::new(s));
            opened_curve_indices.push(curve_num);
        }

        Signature {
            num_curves: CURVES as u32,
            challenges,
            ephemeral_cge,
            opened_curves,
            proof: self
                .sk
                .proof_tree
                .proof_from_leaf_indices(&opened_curve_indices),
        }
    }

    pub fn public_key(&self) -> VerifyingKey {
        self.pk
    }

    pub fn serialized_secret_key(&self) -> Vec<u8> {
        self.sk.serialize()
    }
}

impl VerifyingKey {
    pub fn verify(&self, sig: &Signature, message: &[u8]) -> Result<bool, VerificationFailed> {
        // this can panic,should be changed to gracefully fail
        let challenges = sig
            .challenges
            .chunks_exact(4)
            .map(|x| <i32>::from_be_bytes(x.try_into().unwrap()))
            .collect::<Vec<i32>>();
        if challenges.len() != sig.opened_curves.len() {
            return Err(VerificationFailed);
        }
        let leaf_hashes: &Vec<(u32, HashType)> = &challenges
            .iter()
            .zip(&sig.opened_curves)
            .map(|(c, curve)| {
                let mut hasher = Hasher::default();
                let curve_num = c.unsigned_abs() % sig.num_curves;
                if c.is_negative() {
                    hasher.update(&curve.twist().a.x.retrieve().to_be_bytes());
                } else {
                    hasher.update(&curve.a.x.retrieve().to_be_bytes());
                }
                let label = curve_num + sig.num_curves;
                hasher.update(&label.to_be_bytes());
                hasher.update(&self.merkle);
                let mut result: HashType = HashType::default();
                sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
                (label, result)
            })
            .collect();
        sig.proof.verify(&self.merkle, &self.root, &leaf_hashes)?;
        let mut hasher = Hasher::default();
        let ephemeral_curves = &sig
            .ephemeral_cge
            .par_iter()
            .zip(&sig.opened_curves)
            .map(|(ephem_cge, curve)| ephem_cge.reduce().act_on(curve))
            .collect::<Vec<MontgomeryCurve>>();
        for pk_curve in ephemeral_curves {
            hasher.update(&pk_curve.to_bytes().unwrap());
        }
        hasher.update(message);
        let mut derived = sig.challenges.clone();
        let mut reader = hasher.finalize_xof();
        reader
            .read_exact(derived.as_mut_slice())
            .or(Err(VerificationFailed))?;
        if derived == sig.challenges {
            return Ok(true);
        }
        Err(VerificationFailed)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Signature {
    //size of first 4: 28 + num_rounds * 4 + num_rounds * 40 bytes + num_rounds * 64 bytes
    num_curves: u32,
    challenges: Vec<u8>,                   // num_rounds * 4 bytes
    ephemeral_cge: Vec<ClassGroupElement>, // num_rounds * 40 bytes
    opened_curves: Vec<MontgomeryCurve>,   // num_rounds * 64 bytes
    proof: ClassGroupMerkleProof,
}

impl Signature {
    pub fn serialize(&self) -> Vec<u8> {
        let serialized_num_curves = self.num_curves.to_be_bytes().to_vec();
        let serialized_challenges = self.challenges.clone();
        let num_ephemeral_cge = (self.ephemeral_cge.len() as u32).to_be_bytes().to_vec();
        let serialized_ephemeral_cge = self
            .ephemeral_cge
            .iter()
            .flat_map(|cge| cge.value.retrieve().to_be_bytes())
            .collect::<Vec<u8>>();

        let num_curves = (self.opened_curves.len() as u32).to_be_bytes().to_vec();
        let serialized_curves = self
            .opened_curves
            .iter()
            .flat_map(|curve| curve.to_bytes().expect("bad curve"))
            .collect::<Vec<u8>>();

        let serialized_proof = self.proof.serialize();
        let serialized_challenges_len = (serialized_challenges.len() as u32).to_be_bytes().to_vec();
        let serialized_ephemeral_cge_len = (serialized_ephemeral_cge.len() as u32 + 4)
            .to_be_bytes()
            .to_vec();
        let serialized_curves_len = (serialized_curves.len() as u32 + 4).to_be_bytes().to_vec();
        let serialized_proof_len = (serialized_proof.len() as u32).to_be_bytes().to_vec();
        [
            serialized_num_curves,
            serialized_challenges_len,
            serialized_ephemeral_cge_len,
            serialized_curves_len,
            serialized_proof_len,
            serialized_challenges,
            num_ephemeral_cge,
            serialized_ephemeral_cge,
            num_curves,
            serialized_curves,
            serialized_proof,
        ]
        .concat()
    }

    pub fn deserialize(sig: &Vec<u8>) -> Result<Signature, DeserializationError> {
        let num_curves = <u32>::from_be_bytes(sig[..4].try_into()?);
        let sig = &sig[4..];
        let serialized_challenges_len = <u32>::from_be_bytes(sig[..4].try_into()?) as usize;
        let serialized_ephemeral_cge_len = <u32>::from_be_bytes(sig[4..8].try_into()?) as usize;
        let serialized_curves_len = <u32>::from_be_bytes(sig[8..12].try_into()?) as usize;
        let serialized_proof_len = <u32>::from_be_bytes(sig[12..16].try_into()?) as usize;
        let sig = &sig[16..];
        let serialized_challenges_bytes = &sig[..serialized_challenges_len];
        let sig = &sig[serialized_challenges_len..];
        let serialized_ephemeral_cge_bytes = &sig[..serialized_ephemeral_cge_len];
        let sig = &sig[serialized_ephemeral_cge_len..];
        let serialized_curves_bytes = &sig[..serialized_curves_len];
        let sig = &sig[serialized_curves_len..];
        let serialized_proof_bytes = &sig[..serialized_proof_len];

        let deserialized_challenge = serialized_challenges_bytes.to_vec();
        let deserialized_proof = ClassGroupMerkleProof::deserialize(serialized_proof_bytes)?;

        if serialized_ephemeral_cge_bytes.len() <= 4 {
            return Err(DeserializationError);
        }

        let num_cge = <u32>::from_be_bytes(serialized_ephemeral_cge_bytes[0..4].try_into()?);
        let num_cge_bytes = serialized_ephemeral_cge_bytes.len() - 4;
        if !(num_cge_bytes % 40 == 0 && ((num_cge_bytes as u32 / 40) == num_cge)) {
            return Err(DeserializationError);
        }
        let deserialized_cge: Vec<ClassGroupElement> = serialized_ephemeral_cge_bytes[4..]
            .chunks_exact(40)
            .map(ClassGroupElement::from_be_slice)
            .collect();
        if serialized_curves_bytes.len() <= 4 {
            return Err(DeserializationError);
        }
        let num_opened_curves = <u32>::from_be_bytes(serialized_curves_bytes[0..4].try_into()?);
        let num_curves_bytes = serialized_curves_bytes.len() - 4;
        if !(num_curves_bytes % 64 == 0 && ((num_curves_bytes as u32 / 64) == num_opened_curves)) {
            return Err(DeserializationError);
        }
        let deserialized_curves: Vec<MontgomeryCurve> = serialized_curves_bytes[4..]
            .chunks_exact(64)
            .map(MontgomeryCurve::from_be_slice)
            .collect();

        // Validate public key curves
        for curve in &deserialized_curves {
            let valid = curve.is_nonsingular() && curve.is_supersingular();
            if !valid {
                return Err(DeserializationError);
            }
        }

        Ok(Signature {
            num_curves,
            challenges: deserialized_challenge,
            ephemeral_cge: deserialized_cge,
            proof: deserialized_proof,
            opened_curves: deserialized_curves,
        })
    }
}

mod tests {
    use rand::{thread_rng, RngCore};

    use super::*;

    //I know these are pretty dumb tests, but I don't really have time for more, so, ....
    #[test]
    fn serialization() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = KeyPair::<256, 7>::generate();
        let signature = j.sign(&msg);
        let serialized = signature.serialize();
        println!(
            "Signature Size (256 curves, 7 rounds): {}",
            serialized.len()
        );
        let desig = Signature::deserialize(&serialized).unwrap();
        assert_eq!(signature, desig);
    }

    //#[test]
    fn sig_len() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = KeyPair::<256, 7>::generate();
        let mean_length = (0..1000)
            .map(|_| j.sign(&msg).serialize().len())
            .reduce(|a, b| a + b)
            .unwrap()
            / 1000;
        println!("Signature Size (256 curves, 7 rounds): {}", mean_length);
    }

    #[test]
    fn verify() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = KeyPair::<256, 7>::generate();
        for _ in 0..100 {
            let signature = j.sign(&msg);
            j.pk.verify(&signature, &msg).unwrap();
        }
    }

    #[test]
    #[should_panic]
    fn maul_message() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = KeyPair::<256, 7>::generate();
        let signature = j.sign(&msg);
        msg[512] += 1;
        j.pk.verify(&signature, &msg).unwrap();
    }

    #[test]
    #[should_panic]
    fn maul_signature() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = KeyPair::<256, 7>::generate();
        let signature = j.sign(&msg);
        let mut sig = signature.serialize();
        sig[30] += 1;
        let signature = Signature::deserialize(&sig);
        println!("{}", signature.is_ok());
        match signature {
            Ok(sig) => j.pk.verify(&sig, &msg).unwrap(),
            Err(DeserializationError) => true,
        };
    }

    #[test]
    fn serialize_signing_key() {
        let keypair = KeyPair::<32, 7>::generate();
        let bytes = keypair.sk.serialize();
        let keypair2 = SigningKey::deserialize(bytes).unwrap();
        assert_eq!(keypair.sk, keypair2);
    }

    #[test]
    fn serialize_verifying_key() {
        let keypair = KeyPair::<256, 7>::generate();
        let bytes = keypair.pk.serialize();
        let keypair2 = VerifyingKey::deserialize(&bytes).unwrap();
        assert_eq!(keypair.pk, keypair2);
    }
}
