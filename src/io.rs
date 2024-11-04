//! Provides helper functions to read input from files and write output to the stdout.
//! 
//! Nothing fancy and not intended to be super-user friendly or fool-proof, but just
//! enough for the assignment.

use std::fs::File;
use std::io;
use std::io::BufRead;

/// Reads a hyper-graph edge list from a text file.
/// 
/// Note that this function is not intended to be a function that handles all possible
/// errors in the input file and provide super detailed error messages. Instead this is
/// should be considered a simple implementation that allows us to easily provide input
/// to the core functions for the assignment.
/// 
/// # Arguments
/// * `file_name` - The path of the file to read from. This must be a text file of which
///                 each line represents a (unique) hyper edge and contains the indices
///                 of the vertices incident to that edge, separated by white space.
/// 
/// # Returns
/// On success an [Ok] is returned containing a vector that has for each hyper edge a
/// vector of the indices of the vertices incident to that edge. Otherwise an [Err] is
/// returned containing a [String] describing the error.
/// 
pub fn read_edge_list (file_name: &String) -> Result<Vec<Vec<u32>>, String>
{
    let file = match File::open(file_name)  {
        Ok(file) => file,
        Err(_) => return Err(format!("Could not open {file_name}"))
    };

    let reader = io::BufReader::new(file);

    let mut result = Vec::new();
    for line in reader.lines() {
        let line_str = match line {
            Ok(line_str) => line_str,
            Err(_) => return Err("Error while reading line from input file".to_string())
        };
        let mut hyper_edge = Vec::new();

        for vertex_str in line_str.split(' ') {
            let vertex_id = match vertex_str.parse::<u32>() {
                Ok(n) => n,
                Err(_) => return Err("Input file contained a line that could not be parsed as a list of integers separated by spaces".to_string()),
            };
            hyper_edge.push (vertex_id);
        }
        result.push(hyper_edge);
    }
    Ok(result)
}

/// Prints the provided parity constraints to the stdout.
/// 
/// The contraints are printed per line. Each line contains a list of edges
/// that the corresponding constraint comprises, where each edge is printed as
/// a list surrounded by round brackets of the indices of the vertices incident
/// to that edge.
/// 
/// # Arguments
/// * `parity_constraints` - A vector that contains for each constraint a vector of
///                          indices of edges that it comprises. These indices are
///                          indices of the `edge` vector.
/// * `edge` - A vector that contains for edge a vector of the indices of the vertices
///            that are incident to that vector.
pub fn print_constraints(parity_constraints: &Vec<Vec<usize>>, edges: &Vec<Vec<u32>>) {
    for constr in parity_constraints {
        for e in constr {
            print!("(");
            for vertex in &edges[*e] {
                print! ("{vertex} ");
            }
            print!(") ");
        }
        print! ("\n");
    }
}