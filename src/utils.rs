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

pub struct Kmer {
    solidity: bv::BitVec<u8>,
    kmermasks: Vec<u64>,
    subkmer: Vec<Vec<u64>>,
    max_deep: u8,
    k: u8,
}

impl Kmer {
    pub fn new(solidity: bv::BitVec<u8>, k: u8, max_deep: u8) -> Self {
        Kmer {
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
	    k
        }
    }

    pub fn contains(&self, kmer: u64) -> bool {
        self.bitvec.get(cocktail::kmer::remove_first_bit(cocktail::kmer::cannonical(kmer, self.k)))
    }

    pub fn insert(&mut self, kmer: u64) {
        self.bitvec
            .set(cocktail::kmer::remove_first_bit(cocktail::kmer::cannonical(kmer, self.k)), true);
    }
}

pub fn build_tig(kmer: u64, k: u8, solid: &Kmer, visited: &mut Viewed) -> (String, u64, u64) {
    let mut tig = std::collections::VecDeque::new();
    
    let mut current = kmer;
    for n in cocktail::kmer::kmer2seq(current, k).chars() {
	tig.push_back(n);
    }

    while let Some((pred, ovl_len)) = solid.predecessors(current) {
	if pred.len() != 1 {
	    warn!("prev break on {:?} because multi pred", cocktail::kmer::kmer2seq(current, k));
	    for n in pred {
		warn!("{:?}", cocktail::kmer::kmer2seq(n, k));
	    }
	    break;
	}

	if let Some((succ, _)) = solid.successors(current) {
	    if succ.len() != 1 {
		warn!("prev break on {:?} because multi succ", cocktail::kmer::kmer2seq(current, k));
		for n in succ {
		    warn!("{:?}", cocktail::kmer::kmer2seq(n, k));
		}
		break;
	    }
	}

	add_kmer_in_tig(pred[0], k, ovl_len, &mut tig, true);
	current = pred[0];
	visited.insert(current);
    }
    let begin = current;

    current = kmer;
    
    while let Some((succ, ovl_len)) = solid.successors(current) {
	if succ.len() != 1 {
	    warn!("succ break on {:?} because multi succ ", cocktail::kmer::kmer2seq(current, k));
	    for n in succ {
		warn!("{:?}", cocktail::kmer::kmer2seq(n, k));
	    }
	    break;
	}

	if let Some((pred, _)) = solid.predecessors(current) {
	    if pred.len() != 1 {
		warn!("succ break on {:?} because multi pred", cocktail::kmer::kmer2seq(current, k));
		for n in pred {
		    warn!("{:?}", cocktail::kmer::kmer2seq(n, k));
		}
		break;
	    }
	}

	add_kmer_in_tig(succ[0], k, ovl_len, &mut tig, false);
	current = succ[0];
	visited.insert(current);
    }
    
    let ret_tig = tig.iter().map(|x| *x).collect::<String>();

    (ret_tig, cocktail::kmer::cannonical(begin, k), cocktail::kmer::cannonical(current, k))
}


fn add_kmer_in_tig(kmer: u64, k: u8, not_ovl_len: u8, tig: &mut std::collections::VecDeque<char>, in_front: bool) {
    if in_front {
	let seq = cocktail::kmer::kmer2seq(kmer, k);
	for n in seq[..not_ovl_len as usize].chars().rev() {
	    tig.push_front(n);
	}
    } else {
	let seq = cocktail::kmer::kmer2seq(kmer, k);
	for n in seq[(k - not_ovl_len) as usize..].chars() {
	    tig.push_back(n);
	}
    }
} 

pub fn create_link(tig: &usize, end: &u64, first_ends: &std::collections::HashMap<u64, Vec<usize>>, second_ends: &std::collections::HashMap<u64, Vec<usize>>, other_before: bool, paralelle_tig: &std::collections::HashSet<(usize, usize)>) -> Vec<(usize, char, usize, char)> {
    let mut ret = Vec::new();
    
    if let Some(other_tigs) = first_ends.get(end) {
	for other_tig in other_tigs {
	    if paralelle_tig.contains(&normalize_usize_2tuple((*other_tig, *tig))) {
		continue;
	    }

	    if other_tig == tig {
		continue;
	    }
	    
	    if other_before {
		ret.push((*other_tig, '+', *tig, '+'));
	    } else {
		ret.push((*tig, '-', *other_tig, '+'));
	    }
	}
    }
	
    if let Some(other_tigs) = second_ends.get(end) {
	for other_tig in other_tigs {
	    if paralelle_tig.contains(&normalize_usize_2tuple((*other_tig, *tig))) {
		continue;
	    }
	    
	    if other_tig == tig {
		continue;
	    }

	    if other_before {
		ret.push((*other_tig, '+', *tig, '-'));
	    } else {
		ret.push((*tig, '+', *other_tig, '+'));
	    }
	}
    }

    ret
}

pub fn normalize_u64_2tuple(mut a: (u64, u64)) -> (u64, u64) {
    if a.0 > a.1 { std::mem::swap(&mut a.0, &mut a.1); }    

    a
}

pub fn normalize_usize_2tuple(mut a: (usize, usize)) -> (usize, usize) {
    if a.0 > a.1 { std::mem::swap(&mut a.0, &mut a.1); }

    a
}
