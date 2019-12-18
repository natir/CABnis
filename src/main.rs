/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

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

/* crate declaration */
/* cli management*/
#[macro_use]
extern crate structopt;

/* error and logging */
#[macro_use]
extern crate anyhow;
#[macro_use]
extern crate log;
extern crate thiserror;

extern crate bio;
extern crate bv;
extern crate cocktail;
extern crate itertools;
extern crate petgraph;

/* local mod */
mod cli;
mod error;
mod utils;

/* std use */
use std::io::Write;

/* crate use */
use anyhow::{Context, Result};
use structopt::StructOpt;

/* local use */
use error::Error;

fn main() -> Result<()> {
    env_logger::init();

    let params = cli::Command::from_args();

    info!("Begin of read solidity information");
    let (k, data) = cocktail::io::read_solidity_bitfield(
        std::io::BufReader::new(std::fs::File::open(&params.solidity).with_context(|| {
            Error::CantReadFile {
                filename: params.solidity.clone(),
            }
        })?),
        std::fs::metadata(&params.solidity).unwrap().len(),
    );
    info!("End of read solidity information");

    let solid = utils::Kmer::new(data, k, params.edge_threshold);
    let mut graph: petgraph::graphmap::DiGraphMap<u64, (u8, char, char)> =
        petgraph::graphmap::GraphMap::new();

    info!("Begin of kmer graph building");
    let mut kmer_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.kmer).with_context(|| {
            Error::CantWriteFile {
                filename: params.kmer.clone(),
            }
        })?);

    writeln!(kmer_writer, "H\tVN:Z:1.0")?;
    let mut kmer_set = std::collections::HashSet::new();
    for kmer in 0..(cocktail::kmer::get_kmer_space_size(k) << 1) {
        if !solid.is_solid(kmer) {
            continue;
        }

        let cano = cocktail::kmer::cannonical(kmer, k);

        if let Some((preds, ovl_len)) = solid.predecessors(cano) {
            for predecessor in preds {
                let pred_sign = if cocktail::kmer::parity_even(predecessor) {
                    ('+', '-')
                } else {
                    ('-', '+')
                };

                writeln!(
                    kmer_writer,
                    "L\t{}\t{}\t{}\t+\t{}M",
                    cocktail::kmer::cannonical(predecessor, k),
                    pred_sign.0,
                    cano,
                    k - ovl_len
                )?;

                graph.add_edge(
                    cocktail::kmer::cannonical(predecessor, k),
                    cano,
                    (k - ovl_len, pred_sign.0, '+'),
                );
                graph.add_edge(
                    cano,
                    cocktail::kmer::cannonical(predecessor, k),
                    (k - ovl_len, '-', pred_sign.1),
                );

                kmer_set.insert(cocktail::kmer::cannonical(predecessor, k));
                kmer_set.insert(cano);
            }
        }

        if let Some((succs, ovl_len)) = solid.successors(cano) {
            for successor in succs {
                let succ_sign = if cocktail::kmer::parity_even(successor) {
                    ('+', '-')
                } else {
                    ('-', '+')
                };

                writeln!(
                    kmer_writer,
                    "L\t{}\t+\t{}\t{}\t{}M",
                    cano,
                    cocktail::kmer::cannonical(successor, k),
                    succ_sign.0,
                    k - ovl_len
                )?;

                graph.add_edge(
                    cano,
                    cocktail::kmer::cannonical(successor, k),
                    (k - ovl_len, '+', succ_sign.0),
                );
                graph.add_edge(
                    cocktail::kmer::cannonical(successor, k),
                    cano,
                    (k - ovl_len, succ_sign.1, '-'),
                );

                kmer_set.insert(cocktail::kmer::cannonical(successor, k));
                kmer_set.insert(cano);
            }
        }
    }

    for kmer in kmer_set.iter() {
        writeln!(
            kmer_writer,
            "S\t{}\t{}",
            kmer,
            cocktail::kmer::kmer2seq(*kmer, k)
        )?;
    }
    info!("End of kmer graph building");

    //println!("{:?}", petgraph::dot::Dot::new(&graph));

    info!("Begin of unitig building");
    let mut unitigs = Vec::new();
    let mut visited = utils::Viewed::new(cocktail::kmer::get_kmer_space_size(k));
    let mut ext2tig = std::collections::HashMap::new();

    for origin in graph.nodes() {
        //println!("node {}", node);
        if visited.contains(origin) {
            continue;
        }

        let mut unitig = cocktail::kmer::kmer2seq(origin, k);

        let mut begin = origin;
        let preds = utils::simple_path(
            origin,
            '+',
            &graph,
            petgraph::Direction::Incoming,
            &mut visited,
        );
        for pred in preds {
            let next = if pred.1 == '-' {
                cocktail::kmer::kmer2seq(cocktail::kmer::revcomp(pred.0, k), k)
            } else {
                cocktail::kmer::kmer2seq(pred.0, k)
            };

            unitig += &next[pred.2 as usize..];
            begin = pred.0;
        }

        let mut end = origin;
        let succs = utils::simple_path(
            origin,
            '-',
            &graph,
            petgraph::Direction::Incoming,
            &mut visited,
        );
        for succ in succs {
            let next = if succ.1 == '+' {
                cocktail::kmer::kmer2seq(cocktail::kmer::revcomp(succ.0, k), k)
            } else {
                cocktail::kmer::kmer2seq(succ.0, k)
            };

            unitig = vec![&next[..(k as usize - succ.2 as usize)], &unitig].join("");
            end = succ.0;
        }

        ext2tig.insert(begin, unitigs.len());
        ext2tig.insert(end, unitigs.len());

        unitigs.push((unitig, end, begin));
    }
    info!("End of unitig building");

    info!("Begin of write unitig");
    let mut fasta_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.fasta).with_context(|| {
            Error::CantWriteFile {
                filename: params.fasta.clone(),
            }
        })?);

    for (i, (unitig, begin, end)) in unitigs.iter().enumerate() {
        writeln!(
            fasta_writer,
            ">{} begin:{} end:{} len:{}\n{}",
            i,
            begin,
            end,
            unitig.len(),
            unitig
        )?;
    }
    info!("End of write unitig");

    info!("Begin of unitig graph writing");
    let mut graph_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.gfa).with_context(|| {
            Error::CantWriteFile {
                filename: params.gfa.clone(),
            }
        })?);

    for (i, (unitig, begin, end)) in unitigs.iter().enumerate() {
        writeln!(
            graph_writer,
            "S\t{}\t{}\tbe:Z:{}\ten:Z:{}",
            i, unitig, begin, end
        )?;
    }

    for (i, (_, begin, end)) in unitigs.iter().enumerate() {
        for (source, s_ori, target, t_ori, weight) in utils::get_links(
            i,
            *begin,
            '+',
            &graph,
            petgraph::Direction::Incoming,
            &ext2tig,
            &unitigs,
        ) {
            writeln!(
                graph_writer,
                "L\t{}\t{}\t{}\t{}\t{}M",
                source, s_ori, target, t_ori, weight
            )?;
        }

        for (source, s_ori, target, t_ori, weight) in utils::get_links(
            i,
            *begin,
            '-',
            &graph,
            petgraph::Direction::Incoming,
            &ext2tig,
            &unitigs,
        ) {
            writeln!(
                graph_writer,
                "L\t{}\t{}\t{}\t{}\t{}M",
                source, s_ori, target, t_ori, weight
            )?;
        }

        for (source, s_ori, target, t_ori, weight) in utils::get_links(
            i,
            *end,
            '+',
            &graph,
            petgraph::Direction::Outgoing,
            &ext2tig,
            &unitigs,
        ) {
            writeln!(
                graph_writer,
                "L\t{}\t{}\t{}\t{}\t{}M",
                source, s_ori, target, t_ori, weight
            )?;
        }

        for (source, s_ori, target, t_ori, weight) in utils::get_links(
            i,
            *end,
            '-',
            &graph,
            petgraph::Direction::Outgoing,
            &ext2tig,
            &unitigs,
        ) {
            writeln!(
                graph_writer,
                "L\t{}\t{}\t{}\t{}\t{}M",
                source, s_ori, target, t_ori, weight
            )?;
        }
    }

    info!("End of unitig graph writing");

    Ok(())
}
