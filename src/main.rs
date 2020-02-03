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
extern crate niffler;
extern crate pcon;
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
use itertools::Itertools;

/* local use */
use error::Error;

fn main() -> Result<()> {
    env_logger::init();

    let params = cli::Command::from_args();

    let (k, data) = match &params.subcmd {
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

            (k, data)
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

            (
                subcmd_params.kmer_size,
                count.generate_bitfield(subcmd_params.abundance_min),
            )
        }
    };

    let solid = utils::Kmer::new(data, k, params.edge_threshold);
    let mut visited = utils::Viewed::new(cocktail::kmer::get_kmer_space_size(k), k);
    
    let mut unitigs = Vec::new();
    let mut ext2tig_beg = std::collections::HashMap::new();
    let mut ext2tig_end = std::collections::HashMap::new();
    let mut ends2tig = std::collections::HashMap::new();

    info!("Begin of unitig building");
    for kmer in 0..(cocktail::kmer::get_kmer_space_size(k) << 1) {
        if !solid.is_solid(kmer) {
            continue;
        }

	if visited.contains(kmer) {
	    continue;
	}
	
	visited.insert(kmer);
	let (tig, begin, end) = utils::build_tig(kmer, k, &solid, &mut visited);
	
	ext2tig_beg.entry(begin).or_insert(Vec::new()).push(unitigs.len());
	ext2tig_end.entry(end).or_insert(Vec::new()).push(unitigs.len());
	ends2tig.entry(utils::normalize_u64_2tuple((begin, end))).or_insert(Vec::new()).push(unitigs.len());
	unitigs.push(tig);
    }
    info!("End of unitig building");

    info!("Begin of unitig writting");
    let mut unitigs_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.unitigs).with_context(|| {
            Error::CantWriteFile {
                filename: params.unitigs.clone(),
            }
        })?);

    for (i, unitig) in unitigs.iter().enumerate() {
	writeln!(unitigs_writer, ">{}\tln:i:{}\n{}", i, unitig.len(), unitig)?;
    }
    info!("End of unitig writting");

    for (ext, tig) in ext2tig_beg.iter() {
	warn!("ext {:?} tig {:?}", cocktail::kmer::kmer2seq(*ext, k), tig);
    }
    for (ext, tig) in ext2tig_end.iter() {
	warn!("ext {:?} tig {:?}", cocktail::kmer::kmer2seq(*ext, k), tig);
    }
    
    info!("Begin of unitig graph writting");
    let mut graph_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.graph).with_context(|| {
            Error::CantWriteFile {
                filename: params.graph.clone(),
            }
        })?);


    let mut paralelle_tig = std::collections::HashSet::new();
    for tigs in ends2tig.values() {
	if tigs.len() > 1 {
	    for tigs2 in tigs.iter().combinations(2) {
		paralelle_tig.insert(utils::normalize_usize_2tuple((*tigs2[0], *tigs2[1])));
	    }	
	} 
    }
	
    writeln!(graph_writer, "H\tVN:Z:1.0")?;
    for (i, tig) in unitigs.iter().enumerate() {
	writeln!(graph_writer,"S\t{}\t{}", i, tig)?;
    }
    
    for (begin, tigs) in ext2tig_beg.iter() {
	for tig in tigs {
	    for link in utils::create_link(tig, begin, &ext2tig_end, &ext2tig_beg, true, &paralelle_tig) {
		writeln!(graph_writer, "L\t{}\t{}\t{}\t{}\t14M", link.0, link.1, link.2, link.3)?;
	    }
	}
    }

    for (end, tigs) in ext2tig_end.iter() {
	for tig in tigs { 
	    for link in utils::create_link(tig, end, &ext2tig_end, &ext2tig_beg, false, &paralelle_tig) {
		writeln!(graph_writer, "L\t{}\t{}\t{}\t{}\t14M", link.0, link.1, link.2, link.3)?;
	    }
	}
    }
    info!("End of unitig graph writting");
    
    Ok(())
}
