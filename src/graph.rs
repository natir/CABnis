/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use itertools::Itertools;

fn build_kmermasks(deep: u8, k: u8) -> Vec<u64> {
    let mut kmermasks = Vec::new();

    let mut mask = (1 << (k as u64 * 2)) - 1;

    for _ in 0..(deep) {
        mask >>= 2;
        kmermasks.push(mask)
    }

    kmermasks
}

fn build_subkmer(deep: u8) -> Vec<Vec<u64>> {
    let mut kseq = Vec::new();

    let nucs = "ACTG";
    for i in 0..deep {
        let length = i + 1;
        kseq.push(Vec::new());
        for multi_nucs in (0..length).map(|_| nucs.bytes()).multi_cartesian_product() {
            kseq[i as usize].push(cocktail::kmer::seq2bit(multi_nucs.into_iter().as_slice()));
        }
    }

    kseq
}

pub struct Graph {
    solidity: bv::BitVec<u8>,
    kmermasks: Vec<u64>,
    subkmer: Vec<Vec<u64>>,
    max_deep: u8,
    k: u8,
}

impl Graph {
    pub fn new(solidity: bv::BitVec<u8>, k: u8, max_deep: u8) -> Self {
        Graph {
            solidity,
            kmermasks: build_kmermasks(max_deep, k),
            subkmer: build_subkmer(max_deep),
            max_deep,
            k,
        }
    }

    pub fn is_solid(&self, kmer: u64) -> bool {
        self.solidity.get(cocktail::kmer::remove_first_bit(
            cocktail::kmer::cannonical(kmer, self.k),
        ))
    }

    pub fn successors(&self, kmer: u64) -> Option<(Vec<u64>, u8)> {
        for deep in 0..self.max_deep {
            let prefix = (kmer & self.kmermasks[deep as usize]) << (2 * (deep + 1));

            let mut exist_kmer = Vec::new();

            for suffix in self.subkmer[deep as usize].iter() {
                let next_kmer = prefix ^ suffix;
                if self.is_solid(next_kmer) {
                    exist_kmer.push(next_kmer);
                }
            }

            if !exist_kmer.is_empty() {
                return Some((exist_kmer, deep + 1));
            }
        }

        None
    }

    pub fn predecessors(&self, kmer: u64) -> Option<(Vec<u64>, u8)> {
        for deep in 0..self.max_deep {
            let suffix = kmer >> (2 * (deep + 1));

            let mut exist_kmer = Vec::new();

            for prefix in self.subkmer[deep as usize]
                .iter()
                .map(|x| x << (2 * (self.k - (deep + 1))))
            {
                let next_kmer = prefix ^ suffix;
                if next_kmer == kmer {
                    continue;
                }

                if self.is_solid(next_kmer) {
                    exist_kmer.push(next_kmer);
                }
            }

            if !exist_kmer.is_empty() {
                return Some((exist_kmer, deep + 1));
            }
        }

        None
    }
}

pub struct Viewed {
    bitvec: bv::BitVec<u8>,
    k: u8,
}

impl Viewed {
    pub fn new(len: u64, k: u8) -> Self {
        Viewed {
            bitvec: bv::BitVec::new_fill(false, len),
            k,
        }
    }

    pub fn contains(&self, kmer: u64) -> bool {
        self.bitvec.get(cocktail::kmer::remove_first_bit(
            cocktail::kmer::cannonical(kmer, self.k),
        ))
    }

    pub fn insert(&mut self, kmer: u64) {
        self.bitvec.set(
            cocktail::kmer::remove_first_bit(cocktail::kmer::cannonical(kmer, self.k)),
            true,
        );
    }
}
