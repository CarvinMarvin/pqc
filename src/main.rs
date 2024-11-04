mod linalg;
mod basic_types;
mod io;
mod parity_compiler;
mod graph;

extern crate nalgebra as na;

use io::{read_edge_list, print_constraints};
use graph::SimpleUndirHyperGraph;
use parity_compiler::{compute_k_body_constraints_basis, find_cycles};
use std::env;


fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2
    {
        panic! ("Not enough arguments: provide a file name for the graph file as the first argument.\n\
                 Example graph files can be found in the 'graphs' directory of this project.");
    }
    let edge_list = read_edge_list(&args[1]).unwrap_or_else (|e| {panic!("{e}")});

    let graph = SimpleUndirHyperGraph::from_edge_list(&edge_list);
    let (mat_G,
         mat_P_transp,
         k_body_constraints_basis) = compute_k_body_constraints_basis(&graph);
    let all_possible_parity_constraints: Vec<Vec<usize>> = find_cycles(&graph, 3, 4);

    print!("All possible 3- or 4-body constraints: \n");
    print_constraints (&all_possible_parity_constraints, &edge_list);
    print!("\nA basis for all k-body constraints: \n");
    print_constraints (&k_body_constraints_basis, &edge_list);
    print!("\nLower bound on the number of independent constraints to compile the problem: {:?}\n", k_body_constraints_basis.len());
    print!("\nGenerator matrix: {mat_G}");
    let mat_P = mat_P_transp.transpose();
    print!("Parity matrix: {mat_P}" );
}
