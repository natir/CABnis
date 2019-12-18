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
}

impl Viewed {
    pub fn new(len: u64) -> Self {
        Viewed {
            bitvec: bv::BitVec::new_fill(false, len),
        }
    }

    pub fn contains(&self, kmer: u64) -> bool {
        self.bitvec.get(cocktail::kmer::remove_first_bit(kmer))
    }

    pub fn insert(&mut self, kmer: u64) {
        self.bitvec
            .set(cocktail::kmer::remove_first_bit(kmer), true);
    }
}

pub fn invert_ori(orientation: char) -> char {
    if orientation == '+' {
        '-'
    } else {
        '+'
    }
}

pub fn invert_dir(direction: petgraph::Direction) -> petgraph::Direction {
    match direction {
        petgraph::Direction::Incoming => petgraph::Direction::Outgoing,
        petgraph::Direction::Outgoing => petgraph::Direction::Incoming,
    }
}

pub fn successors(
    current: u64,
    orientation: char,
    graph: &petgraph::graphmap::DiGraphMap<u64, (u8, char, char)>,
    direction: petgraph::Direction,
) -> Vec<u64> {
    graph
        .neighbors_directed(current, direction)
        .filter(|x| graph.edge_weight(current, *x).unwrap().1 == orientation)
        .collect()
}

pub fn simple_path(
    origin: u64,
    orientation: char,
    graph: &petgraph::graphmap::DiGraphMap<u64, (u8, char, char)>,
    direction: petgraph::Direction,
    visited: &mut Viewed,
) -> Vec<(u64, char, u8)> {
    let mut nodes = Vec::new();

    let mut current = origin;
    let mut ori = orientation;

    let mut nexts: Vec<u64> = successors(current, ori, graph, direction);
    let mut preds: Vec<u64> = successors(current, invert_ori(ori), graph, invert_dir(direction));

    while nexts.len() == 1 && preds.len() == 1 {
        let next = nexts[0];

        if visited.contains(next) {
            break;
        }
        visited.insert(next);

        let edge = graph.edge_weight(current, next).unwrap();
        ori = edge.2;
        current = next;

        nodes.push((current, ori, edge.0));

        nexts = successors(current, ori, graph, direction);
        preds = successors(current, invert_ori(ori), graph, invert_dir(direction));
    }

    nodes
}

pub fn get_links(
    tig_id: usize,
    node: u64,
    ori: char,
    graph: &petgraph::graphmap::DiGraphMap<u64, (u8, char, char)>,
    direction: petgraph::Direction,
    ext2tig: &std::collections::HashMap<u64, usize>,
    unitigs: &[(String, u64, u64)],
) -> Vec<(usize, char, usize, char, u8)> {
    let mut links: Vec<(usize, char, usize, char, u8)> = Vec::new();

    for next in successors(node, ori, graph, direction) {
        if let Some(next_id) = ext2tig.get(&next) {
            if let Some((weight, _, _)) = graph.edge_weight(node, next) {
                if let Some((_, next_b, next_e)) = unitigs.get(*next_id) {
                    if *next_b == next {
                        match direction {
                            petgraph::Direction::Incoming => {
                                links.push((tig_id, '-', *next_id, '+', *weight))
                            }
                            petgraph::Direction::Outgoing => {
                                links.push((tig_id, '+', *next_id, '+', *weight))
                            }
                        }
                    } else if *next_e == next {
                        match direction {
                            petgraph::Direction::Incoming => {
                                links.push((*next_id, '+', tig_id, '+', *weight))
                            }
                            petgraph::Direction::Outgoing => {
                                links.push((tig_id, '+', *next_id, '-', *weight))
                            }
                        }
                    }
                }
            }
        }
    }

    links
}
