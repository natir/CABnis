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

/* project use */
use crate::graph;
use crate::utils;

/* crate use */
use anyhow::Result;

#[derive(Debug, PartialEq, PartialOrd, Clone, Hash)]
pub enum Edge {
    Kmer,
    Begin,
    End,
    Both,
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Hash, Copy, Eq, Ord)]
pub enum Node {
    Tig(Tig),
    Kmer(Kmer),
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Hash, Copy, Eq, Ord)]
pub struct Tig {
    pub id: usize,
    pub len: usize,
    pub circular: bool,
}

#[derive(Debug, PartialEq, PartialOrd, Clone, Hash, Copy, Eq, Ord)]
pub struct Kmer {
    pub id: u64,
}

fn build_link(s1: Node, t1: Node, s2: Node, t2:Node, graph: &petgraph::graphmap::UnGraphMap<Node, Edge>) -> Option<(usize, char, usize, char)> {
    if let Node::Tig(first) = s1 {
	if let Node::Tig(second) = t2 {
            if let Some(e1) = graph.edge_weight(s1, t1) {
		if let Some(e2) = graph.edge_weight(s2, t2) {
                    if e1 == &Edge::Begin && e2 == &Edge::Begin {
			Some((first.id, '-', second.id, '+'))
                    } else if e1 == &Edge::Begin && e2 == &Edge::End {
			Some((first.id, '-', second.id, '-'))
                    } else if e1 == &Edge::End && e2 == &Edge::Begin {
			Some((first.id, '+', second.id, '+'))
                    } else if e1 == &Edge::End && e2 == &Edge::End {
			Some((first.id, '+', second.id, '-'))
                    } else {
			None
		    }
		} else {
		    None
		}
            } else {
		None
	    }
	} else {
	    None
	}
    } else {
	None
    }
}

pub fn tig_kmer_tig(
    graph: &petgraph::graphmap::UnGraphMap<Node, Edge>,
) -> std::collections::HashSet<(usize, char, usize, char)> {
    let mut ret = std::collections::HashSet::new();

    for node in graph.nodes() {
        if let Node::Tig(_) = node {
            for nnode in graph.neighbors(node) {
                if let Node::Kmer(_) = nnode {
                    for nnnode in graph.neighbors(nnode) {
                        if nnnode == node {
                            continue;
                        }

			if let Some(link) = build_link(node, nnode, nnode, nnnode, graph) {
			    ret.insert(link);
			}
                    }
                }
            }
        }
    }

    ret
}

pub fn tig_kmer_kmer_tig(
    graph: &petgraph::graphmap::UnGraphMap<Node, Edge>,
) -> std::collections::HashSet<(usize, char, usize, char)> {
    let mut ret = std::collections::HashSet::new();

    for node in graph.nodes() {
        if let Node::Tig(_) = node {
            for nnode in graph.neighbors(node) {
                if let Node::Kmer(_) = nnode {
                    for nnnode in graph.neighbors(nnode) {
                        if let Node::Kmer(_) = nnnode {
                            for nnnnode in graph.neighbors(nnnode) {
                                if nnnnode == node {
                                    continue;
                                }
	
				if let Some(link) = build_link(node, nnode, nnnode, nnnnode, graph) {
				    ret.insert(link);
				}
                            }
                        }
                    }
                }
            }
        }
    }

    ret
}

pub fn write_unitig<W>(
    writer: &mut W,
    k: u8,
    solid: &graph::kmer::Graph,
) -> Result<(
    std::collections::HashMap<(u64, u64), Vec<usize>>,
    Vec<Node>,
    Vec<Node>,
    petgraph::graphmap::UnGraphMap<Node, Edge>,
)>
where
    W: std::io::Write,
{
    let mut tig_counter = 0;
    let mut ext_nodes = Vec::new();
    let mut tig_nodes = Vec::new();
    let mut visited = graph::kmer::Viewed::new(cocktail::kmer::get_kmer_space_size(k), k);
    let mut ends2tig: std::collections::HashMap<(u64, u64), Vec<usize>> =
        std::collections::HashMap::new();
    let mut unitig_graph = petgraph::graphmap::UnGraphMap::new();

    for kmer in 0..cocktail::kmer::get_kmer_space_size(k) {
        if !solid.is_solid(kmer) {
            continue;
        }

        if visited.contains(kmer) {
            continue;
        }

        visited.insert(kmer);
        if let Some((tig, begin, end)) = utils::build_tig(kmer, k, &solid, &mut visited) {
            ends2tig
                .entry(crate::utils::normalize_u64_2tuple((begin, end)))
                .or_insert_with(Vec::new)
                .push(tig_counter);

            let node_tig = graph::unitig::Node::Tig(graph::unitig::Tig {
                id: tig_counter,
                len: tig.len(),
                circular: begin == end,
            });
            let node_begin = graph::unitig::Node::Kmer(graph::unitig::Kmer { id: begin });
            let node_end = graph::unitig::Node::Kmer(graph::unitig::Kmer { id: end });

            unitig_graph.add_node(node_tig);
            tig_nodes.push(node_tig);

            unitig_graph.add_node(node_begin);
            ext_nodes.push(node_begin);

            unitig_graph.add_node(node_end);
            ext_nodes.push(node_end);

            if let Some(edge) = unitig_graph.edge_weight(node_tig, node_begin) {
                if edge == &graph::unitig::Edge::Begin {
                    unitig_graph.add_edge(node_tig, node_begin, graph::unitig::Edge::Both);
                }
            } else {
                unitig_graph.add_edge(node_tig, node_begin, graph::unitig::Edge::Begin);
            }

            if let Some(edge) = unitig_graph.edge_weight(node_tig, node_end) {
                if edge == &graph::unitig::Edge::End {
                    unitig_graph.add_edge(node_tig, node_end, graph::unitig::Edge::End);
                }
            } else {
                unitig_graph.add_edge(node_tig, node_end, graph::unitig::Edge::End);
            }

            writeln!(
                writer,
                ">{} LN:i:{} circular:Z:{} begin:i:{} end:i:{}\n{}",
                tig_counter,
                tig.len(),
                begin == end,
                begin,
                end,
                tig
            )?;

            tig_counter += 1;
        } else {
            continue;
        }
    }

    Ok((ends2tig, ext_nodes, tig_nodes, unitig_graph))
}

pub fn add_missing_edge(
    ext_nodes: Vec<Node>,
    solid: graph::kmer::Graph,
    k: u8,
    mut unitig_graph: petgraph::graphmap::UnGraphMap<Node, Edge>,
) -> petgraph::graphmap::UnGraphMap<Node, Edge> {
    for node in ext_nodes {
        if let graph::unitig::Node::Kmer(n) = node {
            if let Some((succs, _)) = solid.successors(n.id) {
                for succ in succs {
                    let cano = cocktail::kmer::cannonical(succ, k);
                    let node_succ = graph::unitig::Node::Kmer(graph::unitig::Kmer { id: cano });
                    if unitig_graph.contains_node(node_succ) {
                        unitig_graph.add_edge(
                            graph::unitig::Node::Kmer(n),
                            node_succ,
                            graph::unitig::Edge::Kmer,
                        );
                    }
                }
            }

            if let Some((preds, _)) = solid.predecessors(n.id) {
                for pred in preds {
                    let cano = cocktail::kmer::cannonical(pred, k);
                    let node_pred = graph::unitig::Node::Kmer(graph::unitig::Kmer { id: cano });
                    if unitig_graph.contains_node(node_pred) {
                        unitig_graph.add_edge(
                            graph::unitig::Node::Kmer(n),
                            node_pred,
                            graph::unitig::Edge::Kmer,
                        );
                    }
                }
            }
        }
    }

    unitig_graph
}
