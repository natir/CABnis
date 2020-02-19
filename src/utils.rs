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
use anyhow::{Context, Result};
use niffler;

/* local mod */
use crate::cli;
use crate::error::Error;
use crate::graph;

pub fn build_tig(
    kmer: u64,
    k: u8,
    solid: &graph::kmer::Graph,
    visited: &mut graph::kmer::Viewed,
) -> Option<(std::collections::VecDeque<u8>, u64, u64)> {
    let mut tig = std::collections::VecDeque::new();

    let mut current = kmer;
    for n in cocktail::kmer::kmer2seq(current, k).bytes() {
        tig.push_back(n);
    }

    /* if a unitig with size equal to k and nb_pred < 2 || nb_succ < 2 it's not a valid unitig */
    let mut nb_pred = 0;
    let mut nb_succ = 0;

    while let Some((pred, ovl_len)) = solid.predecessors(current) {
        nb_pred = pred.len();
        if pred.len() != 1 {
            break;
        }

        if let Some((succ, _)) = solid.successors(current) {
            if succ.len() != 1 {
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
        nb_succ = succ.len();
        if succ.len() != 1 {
            break;
        }

        if let Some((pred, _)) = solid.predecessors(current) {
            if pred.len() != 1 {
                break;
            }
        }

        add_kmer_in_tig(succ[0], k, ovl_len, &mut tig, false);
        current = succ[0];
        visited.insert(current);
    }

    if current == begin && (nb_pred < 2 || nb_succ < 2) {
        return None;
    }

    Some((
        tig,
        cocktail::kmer::cannonical(begin, k),
        cocktail::kmer::cannonical(current, k),
    ))
}

fn add_kmer_in_tig(
    kmer: u64,
    k: u8,
    not_ovl_len: u8,
    tig: &mut std::collections::VecDeque<u8>,
    in_front: bool,
) {
    if in_front {
        let seq = cocktail::kmer::kmer2seq(kmer, k);
        for n in seq[..not_ovl_len as usize].bytes().rev() {
            tig.push_front(n);
        }
    } else {
        let seq = cocktail::kmer::kmer2seq(kmer, k);
        for n in seq[(k - not_ovl_len) as usize..].bytes() {
            tig.push_back(n);
        }
    }
}

pub fn normalize_u64_2tuple(mut a: (u64, u64)) -> (u64, u64) {
    if a.0 > a.1 {
        std::mem::swap(&mut a.0, &mut a.1);
    }

    a
}

pub fn normalize_usize_2tuple(mut a: (usize, usize)) -> (usize, usize) {
    if a.0 > a.1 {
        std::mem::swap(&mut a.0, &mut a.1);
    }

    a
}

pub fn get_count(params: &cli::Command) -> Result<(u8, bv::BitVec<u8>)> {
    match &params.subcmd {
        cli::SubCommand::Count(subcmd_params) => {
            info!("Begin of read solidity information");

            let (k, data) = cocktail::io::read_solidity_bitfield(
                std::io::BufReader::new(std::fs::File::open(&subcmd_params.input).with_context(
                    || Error::CantReadFile {
                        filename: subcmd_params.input.clone(),
                    },
                )?),
                std::fs::metadata(&subcmd_params.input).unwrap().len(),
            );

            info!("End of read solidity information");

            Ok((k, data))
        }
        cli::SubCommand::Reads(subcmd_params) => {
            info!("Begin of kmer counting");

            let mut count = pcon::count::Count::new(subcmd_params.kmer_size, 8);

            let (reader, _) = niffler::get_reader(Box::new(std::io::BufReader::new(
                std::fs::File::open(&subcmd_params.input).with_context(|| Error::CantReadFile {
                    filename: subcmd_params.input.clone(),
                })?,
            )))?;

            let fasta_reader = bio::io::fasta::Reader::new(reader);

            for record in fasta_reader.records() {
                let result = record.with_context(|| Error::ReadingError {
                    filename: subcmd_params.input.clone(),
                })?;

                count.add_sequence(result.seq());
            }

            count.clean_buckets();

            info!("End of kmer counting");

            Ok((
                subcmd_params.kmer_size,
                count.generate_bitfield(subcmd_params.abundance_min),
            ))
        }
    }
}
