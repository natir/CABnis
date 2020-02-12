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
mod graph;
mod utils;

/* std use */
use std::io::Write;

/* crate use */
use anyhow::{Context, Result};
use itertools::Itertools;
use structopt::StructOpt;

/* local use */
use error::Error;

fn main() -> Result<()> {
    env_logger::init();

    let params = cli::Command::from_args();

    let (k, data) = utils::get_count(&params)?;

    let solid = graph::kmer::Graph::new(data, k, params.edge_threshold);

    if let Some(out_path) = &params.kmer {
        info!("Begin of kmer graph building");
        let mut kmer_writer =
            std::io::BufWriter::new(std::fs::File::create(&out_path).with_context(|| {
                Error::CantWriteFile {
                    filename: out_path.to_string(),
                }
            })?);

        graph::kmer::write_kmer_graph(&mut kmer_writer, k, &solid)?;
        info!("End of kmer graph building");
    }

    info!("Begin of unitig building");

    let mut unitigs_writer =
        std::io::BufWriter::new(std::fs::File::create(&params.unitigs).with_context(|| {
            Error::CantWriteFile {
                filename: params.unitigs.clone(),
            }
        })?);

    let (ends2tig, mut unitig_graph) = graph::unitig::write_unitig(&mut unitigs_writer, k, &solid)?;
    info!("End of unitig building");

    info!("Begin of unitg graph building");
    unitig_graph = graph::unitig::add_missing_edge(solid, k, unitig_graph);
    info!("End of unitig graph building");

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

    info!("\tBegin of S record writing");
    for node in unitig_graph.nodes() {
        if let graph::unitig::Node::Tig(n) = node {
            writeln!(
                graph_writer,
                "S\t{}\t*\tLN:i:{}\tcircular:Z:{}",
                n.id, n.len, n.circular
            )?;
        }
    }

    info!("\tEnd of S record writing");

    info!("\tBegin of L record writing");
    for node in unitig_graph.nodes() {
        if let graph::unitig::Node::Tig(n) = node {
            if n.circular {
                writeln!(graph_writer, "L\t{}\t-\t{}\t+\t14M", n.id, n.id)?;
            }
        }
    }

    for link in graph::unitig::tig_kmer_tig(&unitig_graph) {
        if paralelle_tig.contains(&utils::normalize_usize_2tuple((link.0, link.2))) {
            continue;
        }

        writeln!(
            graph_writer,
            "L\t{}\t{}\t{}\t{}\t14M",
            link.0, link.1, link.2, link.3
        )?;
    }

    for link in graph::unitig::tig_kmer_kmer_tig(&unitig_graph) {
        if paralelle_tig.contains(&utils::normalize_usize_2tuple((link.0, link.2))) {
            continue;
        }

        writeln!(
            graph_writer,
            "L\t{}\t{}\t{}\t{}\t14M",
            link.0, link.1, link.2, link.3
        )?;
    }

    info!("\tEnd of L record writing");
    info!("End of unitig graph writting");

    Ok(())
}
