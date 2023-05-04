extern crate itertools;
use itertools::Itertools;
use std::collections::HashSet;

// Checks that a set of edges are paiwise disjoint
fn _is_a_matching(edges: &Vec<&(usize, usize)>) -> bool {
    let mut witness: HashSet<usize> = HashSet::new();
    for (a, b) in edges.iter() {
        witness.insert(*a);
        witness.insert(*b);
    }
    return witness.len() == edges.len() * 2;
}

fn _powerset<T>(s: &[T]) -> Vec<Vec<&T>> {
    (0..2usize.pow(s.len() as u32))
        .map(|i| {
            s.iter()
                .enumerate()
                .filter(|&(t, _)| (i >> t) % 2 == 1)
                .map(|(_, element)| element)
                .collect()
        })
        .collect()
}
pub struct Topology {
    pub nvertices: usize,
    pub block_size: usize,
    pub edges: Vec<(usize, usize)>,
}

impl Topology {
    /// Returns all the possible matchings in the topology
    pub fn get_all_matchings(&self) -> Vec<Vec<(usize, usize)>> {
        let mut all_matchings = Vec::new();
        for matching_size in 1..(self.nvertices / 2 + 1) {
            let matchings: Vec<Vec<(usize, usize)>> = self
                .edges
                .iter()
                .combinations(matching_size)
                .filter(_is_a_matching)
                .map(|matching| {
                    matching
                        .iter()
                        .map(|&edge| *edge)
                        .collect::<Vec<(usize, usize)>>()
                })
                .collect();
            all_matchings.extend_from_slice(&matchings);
        }
        return all_matchings;
    }

    /// Returns all depth-1 CNOT circuits available in the topology
    /// We first generate all possible matchings, and for each matching,
    /// all possible CNOT orientations
    pub fn get_all_operators(&self) -> Vec<Vec<(usize, usize)>> {
        let mut all_operators = Vec::new();
        for matching in self.get_all_matchings() {
            let k = matching.len();
            let operator: Vec<Vec<(usize, usize)>> = (0..(1 << k))
                .map(|i| {
                    matching
                        .iter()
                        .enumerate()
                        .map(|(index, (a, b))| {
                            if ((i >> index) & 1) == 0 {
                                (*b, *a)
                            } else {
                                (*a, *b)
                            }
                        })
                        .collect()
                })
                .collect();
            all_operators.extend_from_slice(&operator);
        }
        return all_operators;
    }
    /// These are all the topologies benchmarked in the paper
    /// It is quite simple to add a new one :)
    pub fn from_string(name: &str) -> Self {
        match name {
            "line_2" => Self {
                nvertices: 4,
                block_size: 2,
                edges: vec![(0, 1), (1, 2), (2, 3)],
            },
            "ladder_2" => Self {
                nvertices: 4,
                block_size: 2,
                edges: vec![(0, 1), (0, 2), (1, 3), (2, 3)],
            },
            "ladder_2_diagonals" => Self {
                nvertices: 4,
                block_size: 2,
                edges: vec![(0, 1), (0, 2), (1, 3), (2, 3), (0, 3), (1, 2)],
            },
            "ladder_3" => Self {
                nvertices: 6,
                block_size: 3,
                edges: vec![(0, 1), (0, 3), (1, 4), (1, 2), (2, 5), (3, 4), (4, 5)],
            },
            "ladder_3_diagonals" => Self {
                nvertices: 6,
                block_size: 3,
                edges: vec![
                    (0, 1),
                    (0, 3),
                    (1, 4),
                    (1, 2),
                    (2, 5),
                    (3, 4),
                    (4, 5),
                    (0, 4),
                    (1, 3),
                    (1, 5),
                    (2, 4),
                ],
            },
            "all_to_all_3" => Self {
                nvertices: 6,
                block_size: 3,
                edges: vec![
                    (0, 1),
                    (0, 2),
                    (0, 3),
                    (0, 4),
                    (0, 5),
                    (1, 2),
                    (1, 3),
                    (1, 4),
                    (1, 5),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                    (3, 4),
                    (3, 5),
                    (4, 5),
                ],
            },
            "ladder_4" => Self {
                nvertices: 8,
                block_size: 4,
                edges: vec![
                    (0, 1),
                    (0, 4),
                    (1, 2),
                    (1, 5),
                    (2, 3),
                    (2, 6),
                    (3, 7),
                    (4, 5),
                    (5, 6),
                    (6, 7),
                ],
            },
            "ladder_4_diagonals" => Self {
                nvertices: 8,
                block_size: 4,
                edges: vec![
                    (0, 1),
                    (0, 4),
                    (1, 2),
                    (1, 5),
                    (2, 3),
                    (2, 6),
                    (3, 7),
                    (4, 5),
                    (5, 6),
                    (6, 7),
                    (0, 5),
                    (1, 4),
                    (1, 6),
                    (2, 5),
                    (2, 7),
                    (3, 6),
                ],
            },
            "all_to_all_4" => Self {
                nvertices: 8,
                block_size: 4,
                edges: vec![
                    (0, 1),
                    (0, 2),
                    (0, 3),
                    (0, 4),
                    (0, 5),
                    (0, 6),
                    (0, 7),
                    (1, 2),
                    (1, 3),
                    (1, 4),
                    (1, 5),
                    (1, 6),
                    (1, 7),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                    (2, 6),
                    (2, 7),
                    (3, 4),
                    (3, 5),
                    (3, 6),
                    (3, 7),
                    (4, 5),
                    (4, 6),
                    (4, 7),
                    (5, 6),
                    (5, 7),
                    (6, 7),
                ],
            },
            "grid" => Self {
                nvertices: 8,
                block_size: 4,
                edges: vec![
                    (0, 1),
                    (0, 2),
                    (1, 3),
                    (2, 3),
                    (2, 4),
                    (3, 5),
                    (4, 5),
                    (4, 6),
                    (5, 7),
                    (6, 7),
                ],
            },
            "grid_diagonals" => Self {
                nvertices: 8,
                block_size: 4,
                edges: vec![
                    (0, 1),
                    (0, 2),
                    (1, 3),
                    (2, 3),
                    (2, 4),
                    (3, 5),
                    (4, 5),
                    (4, 6),
                    (5, 7),
                    (6, 7),
                    (0, 3),
                    (1, 2),
                    (2, 5),
                    (3, 4),
                    (4, 7),
                    (5, 6),
                ],
            },
            _ => panic!("Unknown topology!"),
        }
    }
    /// Returns the subgraph induced by a single block in the topology
    pub fn get_block_topology(&self) -> Self {
        Self {
            nvertices: self.nvertices / 2,
            block_size: self.block_size,
            edges: self
                .edges
                .iter()
                .filter(|(a, b)| (*a < self.nvertices / 2) & (*b < self.nvertices / 2))
                .map(|a| *a)
                .collect(),
        }
    }
}
