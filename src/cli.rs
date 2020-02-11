/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(
    version = "0.1",
    author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>",
    name = "cabanis",
    about = "Use solid kmer and sequence to build a compacted ABruijn graph in gfa format."
)]
pub struct Command {
    #[structopt(
        short = "g",
        long = "graph",
        required = true,
        help = "path of gfa output file"
    )]
    pub graph: String,

    #[structopt(
        short = "u",
        long = "unitigs",
        required = true,
        help = "path of fasta output file"
    )]
    pub unitigs: String,

    #[structopt(short = "k", long = "kmer", help = "path of kmer graph output file")]
    pub kmer: Option<String>,

    #[structopt(
        short = "t",
        long = "edge-weight-threshold",
        default_value = "5",
        help = "All edge in graph are lower than this value"
    )]
    pub edge_threshold: u8,

    #[structopt(subcommand)]
    pub subcmd: SubCommand,
}

#[derive(StructOpt, Debug)]
pub enum SubCommand {
    #[structopt(about = "Generate unitig graph from pcon count")]
    Count(Count),
    #[structopt(about = "Generate unitig graph from reads")]
    Reads(Reads),
}

#[derive(StructOpt, Debug)]
pub struct Count {
    #[structopt(
        short = "i",
        long = "input",
        required = true,
        help = "path to pcon solidity file"
    )]
    pub input: String,
}

#[derive(StructOpt, Debug)]
pub struct Reads {
    #[structopt(
        short = "i",
        long = "input",
        required = true,
        help = "path to reads file"
    )]
    pub input: String,

    #[structopt(
        short = "k",
        long = "kmer-size",
        required = true,
        help = "kmer size, if kmer size is even real value is equal to k-1, max value 31"
    )]
    pub kmer_size: u8,

    #[structopt(
        short = "a",
        long = "abudance-min",
        default_value = "1",
        required = true,
        help = "write only kmer with abudance is higher than this parametre"
    )]
    pub abundance_min: u8,
}
