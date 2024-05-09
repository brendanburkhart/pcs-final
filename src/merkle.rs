use std::collections::{HashMap, HashSet, VecDeque};

use crypto_bigint::{rand_core::RngCore, Encoding};
use rand::thread_rng;
use sha3::digest::{ExtendableOutput, Update};

use crate::constants::{DeserializationError, HashType, Hasher, VerificationFailed, HASH_SIZE};
use crate::montgomery::MontgomeryCurve;

#[derive(Debug, Clone, PartialEq)]
pub struct ClassGroupMerkleTree {
    root: HashType,
    merkle_key: HashType,
    layers: Vec<Vec<(u32, HashType)>>,
}

impl ClassGroupMerkleTree {
    pub fn from_leaves(leaves: &[MontgomeryCurve]) -> Self {
        let n_leaves = leaves.len();
        assert!(n_leaves.is_power_of_two());
        let depth = n_leaves.ilog2();
        let mut merkle_key = [0u8; 32];
        thread_rng().fill_bytes(&mut merkle_key);

        let mut layers: Vec<Vec<(u32, HashType)>> = Vec::new();
        layers.push(
            leaves
                .iter()
                .enumerate()
                .map(|(i, curve)| {
                    let label = (1 << (depth - layers.len() as u32)) + i as u32;
                    let mut result: HashType = HashType::default();
                    let mut hasher = Hasher::default();
                    hasher.update(&curve.a.x.retrieve().to_be_bytes());
                    hasher.update(&label.to_be_bytes());
                    hasher.update(&merkle_key);
                    sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
                    (label, result)
                })
                .collect::<Vec<(u32, HashType)>>(),
        );
        while layers[layers.len() - 1].len() != 1 {
            layers.push(
                layers[layers.len() - 1]
                    .chunks(2)
                    .enumerate()
                    .map(|(i, chunk)| {
                        let mut result: HashType = HashType::default();
                        let mut hasher = Hasher::default();
                        let label = (1 << (depth - layers.len() as u32)) + i as u32;
                        hasher.update(&chunk[0].1);
                        hasher.update(&chunk[1].1);
                        hasher.update(&label.to_be_bytes());
                        hasher.update(&merkle_key);
                        sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
                        (label, result)
                    })
                    .collect::<Vec<(u32, HashType)>>(),
            );
        }
        ClassGroupMerkleTree {
            root: layers.last().unwrap()[0].1,
            merkle_key,
            layers,
        }
    }

    //TODO: remove redundant internal nodes.
    pub fn proof_from_leaf_indices(&self, leaf_indices: &[usize]) -> ClassGroupMerkleProof {
        let mut proof_tree: Vec<Vec<(u32, HashType)>> = Vec::new();
        let mut indices = leaf_indices
            .iter()
            .map(|x| (*x as u32, (x + (1 << self.depth())) as u32))
            .collect::<Vec<(u32, u32)>>()
            .to_vec();

        loop {
            let to_remove = HashSet::<(u32, u32)>::from_iter(indices.iter().cloned());
            let mut siblings = indices
                .iter()
                .map(|(x, label)| {
                    if (x % 2) == 0 {
                        (x + 1, label + 1)
                    } else {
                        (x - 1, label - 1)
                    }
                })
                .collect::<Vec<(u32, u32)>>();
            siblings.dedup();
            siblings.retain(|x| !to_remove.contains(x));
            proof_tree.push(
                siblings
                    .iter()
                    .map(|(x, label)| {
                        (
                            *label,
                            self.layers[self.depth() - label.ilog2() as usize][*x as usize].1,
                        )
                    })
                    .collect::<Vec<(u32, HashType)>>(),
            );

            indices = indices
                .iter()
                .map(|(x, label)| {
                    if (x % 2) == 0 {
                        (x / 2, label / 2)
                    } else {
                        ((x - 1) / 2, (label - 1) / 2)
                    }
                })
                .collect();
            indices.dedup();
            if indices[0].1 == 1 {
                break;
            }
        }
        let proof: Vec<(u32, HashType)> = proof_tree.iter().flatten().copied().collect();
        ClassGroupMerkleProof { proof }
    }

    pub fn depth(&self) -> usize {
        self.layers.len() - 1
    }

    pub fn root(&self) -> HashType {
        self.root
    }

    pub fn merkle_key(&self) -> HashType {
        self.merkle_key
    }

    pub fn leaves(&self) -> Vec<(u32, HashType)> {
        self.layers[0].clone()
    }

    pub fn serialize(&self) -> Vec<u8> {
        let hash_size = HASH_SIZE.to_be_bytes().to_vec();
        let serialized_proof_tree = self
            .layers
            .iter()
            .flat_map(|vec| {
                let flattened = vec
                    .iter()
                    .flat_map(|(idx, hash)| [idx.to_be_bytes().to_vec(), hash.to_vec()].concat())
                    .collect::<Vec<u8>>();
                [(flattened.len() as u32).to_be_bytes().to_vec(), flattened].concat()
            })
            .collect::<Vec<u8>>();
        [
            hash_size,
            self.root.to_vec(),
            self.merkle_key.to_vec(),
            serialized_proof_tree,
        ]
        .concat()
    }

    pub fn deserialize(key: &[u8]) -> Result<ClassGroupMerkleTree, DeserializationError> {
        let hash_size = <u32>::from_be_bytes(key[..4].try_into()?) as usize;
        let root: HashType = key[4..4 + hash_size].try_into()?;
        let merkle_key: HashType = key[4 + hash_size..4 + hash_size + hash_size]
            .try_into()
            .unwrap();
        let mut key = &key[4 + hash_size + hash_size..];
        let mut layers: Vec<Vec<(u32, HashType)>> = Vec::new();
        loop {
            let level_size = <u32>::from_be_bytes(key[..4].try_into()?) as usize;
            layers.push(
                key[4..4 + level_size]
                    .chunks_exact(4 + HASH_SIZE as usize)
                    .map(|x| {
                        (
                            <u32>::from_be_bytes(x[..4].try_into().unwrap()),
                            <HashType>::try_from(&x[4..]).unwrap(),
                        )
                    })
                    .collect::<Vec<(u32, HashType)>>(),
            );
            key = &key[4 + level_size..];
            if (level_size as u32) == (HASH_SIZE + 4) {
                break;
            }
        }
        Ok(ClassGroupMerkleTree {
            root,
            merkle_key,
            layers,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClassGroupMerkleProof {
    proof: Vec<(u32, HashType)>,
}

impl ClassGroupMerkleProof {
    pub fn verify(
        &self,
        merkle_key: &HashType,
        root: &HashType,
        leaf_hashes: &[(u32, HashType)],
    ) -> Result<(), VerificationFailed> {
        let mut level: VecDeque<(u32, HashType)> = VecDeque::from(leaf_hashes.to_vec());
        let mut tree: HashMap<u32, HashType> = HashMap::new();
        for (lab, hash) in &self.proof {
            tree.insert(*lab, *hash);
        }
        loop {
            //some 3am code right here...
            assert!(!level.is_empty());
            let (label, hash) = level.pop_front().unwrap();
            let mut result: HashType = HashType::default();
            let mut hasher = Hasher::default();
            if (label % 2) == 0 {
                hasher.update(&hash);
                let sibling = tree.get(&(label + 1));
                match sibling {
                    Some(h) => hasher.update(h),
                    None => {
                        match level.iter().find(|(l, _)| *l == (label + 1)) {
                            Some((_, hash)) => hasher.update(hash),
                            None => return Err(VerificationFailed),
                        };
                    }
                }
            } else {
                let sibling = tree.get(&(label - 1));
                match sibling {
                    Some(h) => hasher.update(h),
                    None => {
                        match level.iter().find(|(l, _)| *l == (label - 1)) {
                            Some((_, hash)) => hasher.update(hash),
                            None => return Err(VerificationFailed),
                        };
                    }
                }
                hasher.update(&hash);
            }
            hasher.update(&(label / 2).to_be_bytes());
            hasher.update(merkle_key);
            sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
            let next = (label / 2, result);
            level.push_back(next);
            tree.insert(next.0, next.1);
            tree.insert(label, hash);
            if next.0 == 1 {
                break;
            }
        }
        let putative_root = tree.get(&1).unwrap();
        if putative_root == root {
            Ok(())
        } else {
            Err(VerificationFailed)
        }
    }

    pub fn serialize(&self) -> Vec<u8> {
        let mut result = (self.proof.len() as u32).to_be_bytes().to_vec();
        let mut proof = self
            .proof
            .iter()
            .flat_map(|(node, hash)| [node.to_be_bytes().as_slice(), hash].concat())
            .collect::<Vec<u8>>();
        result.append(&mut proof);
        result
    }

    pub fn deserialize(
        putative_proof: &[u8],
    ) -> Result<ClassGroupMerkleProof, DeserializationError> {
        let proof_len = putative_proof.len() as u32;
        let data_len = proof_len - 4;
        if proof_len > 4 && ((data_len % (HASH_SIZE + 4)) == 0) {
            let num_nodes = <u32>::from_be_bytes(putative_proof[0..4].try_into()?);
            if (data_len / 36) == num_nodes {
                let proof = ClassGroupMerkleProof {
                    proof: putative_proof[4..]
                        .chunks_exact((HASH_SIZE + 4) as usize)
                        .map(|c| {
                            let node_idx = <u32>::from_be_bytes(c[0..4].try_into().unwrap());
                            let node_hash: HashType = c[4..].try_into().unwrap();
                            (node_idx, node_hash)
                        })
                        .collect::<Vec<(u32, HashType)>>(),
                };
                return Ok(proof);
            }
        }
        Err(DeserializationError)
    }
}

mod tests {
    use crypto_bigint::Random;

    use crate::constants::ModP;

    use super::*;

    //I know these are pretty dumb tests, but I don't really have time for more, so, ....

    #[test]
    fn test() {
        let j: Vec<MontgomeryCurve> = (0..4)
            .map(|_| MontgomeryCurve::new(ModP::random(&mut thread_rng())))
            .collect();
        let mt = ClassGroupMerkleTree::from_leaves(&j);
        let proof = mt.proof_from_leaf_indices(&[0, 3]);
        assert_eq!(proof.proof[0].1, {
            let mut hasher = Hasher::default();
            let mut result = HashType::default();
            hasher.update(&j[1].a.x.retrieve().to_be_bytes());
            hasher.update(&5u32.to_be_bytes());
            hasher.update(&mt.merkle_key);
            sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
            result
        });
        assert_eq!(proof.proof[1].1, {
            let mut hasher = Hasher::default();
            let mut result = HashType::default();
            hasher.update(&j[2].a.x.retrieve().to_be_bytes());
            hasher.update(&6u32.to_be_bytes());
            hasher.update(&mt.merkle_key);
            sha3::digest::XofReader::read(&mut hasher.finalize_xof(), &mut result);
            result
        });
        let serialized = proof.serialize();
        let proof = ClassGroupMerkleProof::deserialize(serialized.as_slice()).unwrap();
        let result = proof.verify(
            &mt.merkle_key,
            &mt.root(),
            &Vec::from([mt.leaves()[0], mt.leaves()[3]]),
        );
        result.unwrap();
    }

    #[test]
    fn serialize_merkle_tree() {
        let j: Vec<MontgomeryCurve> = (0..4)
            .map(|_| MontgomeryCurve::new(ModP::random(&mut thread_rng())))
            .collect();
        let mt = ClassGroupMerkleTree::from_leaves(&j);
        let serialized = mt.serialize();
        let deserialized = ClassGroupMerkleTree::deserialize(&serialized).unwrap();
        assert_eq!(mt, deserialized);
    }
}
