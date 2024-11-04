//! Contains the core algorithms for the assignment, which are specific for
//! the application of a parity compiler. 
//! Please note that [crate::linalg] also contains core algorithms for the assignment,
//! but these could also be used outside the application of a parity compiler.

use crate::basic_types::ZMod2;
use crate::basic_types::*;
use crate::linalg::null_space;

use na::base::DMatrix;
use std::iter::zip;

use crate::graph::SimpleUndirHyperGraph;

/// Computes the generator matrix corresponding to the parity transformation
/// for the undirectional hyper graph `graph`.
/// 
/// For a definition of the generator matrix see the
/// [Ender et al., "Parity Quantum Optimization: Compiler"](https://doi.org/10.22331/q-2023-03-17-950)
/// paper.
/// 
/// # Returns
/// The generator matrix, where an element at row *i* and column *j*, corresponds to
/// the participation of the vertex with index *i* in the edge with index *j* in
/// `graph`. Observe that this matrix has `graph.num_edges()` rows and
/// `graph.num_vertices()` columns.
pub fn generator_matrix (graph: &SimpleUndirHyperGraph) -> DMatrix<ZMod2>
{
    let mut mat_G = DMatrix::<ZMod2>::zeros(graph.num_vertices(), graph.num_edges());

    for (hyper_edge, mut col) in zip(graph.edges_iter(), mat_G.column_iter_mut()) {
        for vertex in hyper_edge.incident_vertices_iter() {
            col[vertex.get_index()] = ZMod2::one()
        }
    }

    mat_G
}


/// Converts a parity check matrix to a vector of parity constraints.
/// 
/// # Arguments
/// * `mat_P_transp` -  A matrix where an element at row *j* and column *i*,
///                     corresponds to the participation of the edge with
///                     index *j* in *i*'th constraint/cycle. Here we are
///                     referring to edges and cycles of the original graph
///                     this transposed parity check matrix was constructed
///                     for.
/// 
/// # Returns
/// Returns a vector for which the element with index *i* is a vector of
/// indices of the edges that the *i*'th constrains/cycle consists of.
pub fn extract_parity_constraints (mat_P_transp: &DMatrix<ZMod2>) -> Vec<Vec<usize>>
{ 
    let mut result = Vec::new();

    for col in mat_P_transp.column_iter() {
        let mut parity_constr = Vec::new();
        for r in 0 .. mat_P_transp.nrows() {
            if col[r].is_one() {
                parity_constr.push(r);
            }
        }
        result.push(parity_constr);
    }
    result
}

/// Finds all simple cycles in an undirectional hypergraph `graph` of which the
/// number of edges is in the interval (`min_cycle_len` ... `max_cycle_len`).
/// 
/// Here a *cycle* of an undirectional graph *g* is defined as any subset *C* of
/// the set of all edges of *g* such that each vertex of *g* is incident to an
/// even number of edges of *C*.
/// A *simple cycle* of an undirectional graph *g* is defined as any cycle *C*
/// that cannot be paritioned into subsets *C<sub>1</sub>*,...,*C<sub>m</sub>*
/// (for some *m* > 1) such that *C<sub>1</sub>*,...,*C<sub>m</sub>* are also
/// cycles.
/// Note that for `max_cycle_len` <= 4, the set of cycles of a graph is equal to
/// the set of simple cycles.
/// 
/// # Arguments
/// * `graph` - The input graph.
/// * `min_cycle_len` - The minimum length of simple cycle.
/// * `max_cycle_len` - The maximum length of simple cycle.
/// 
/// # Returns
/// Returns a vector that contains for each found simple cycle a vector of indices
/// of the edges belonging to that simple cycle.
pub fn find_cycles (graph: &SimpleUndirHyperGraph,
                    min_cycle_len: usize,
                    max_cycle_len: usize) -> Vec<Vec<usize>> {
    // See the documentation of find_cycles_from_path for a explanation of the
    // following variables:
    let mut cycles = Vec::new();
    let mut helper = vec![false; graph.num_vertices()];
    let mut is_not_satisfied = vec![false; graph.num_vertices()];
    let mut in_path = vec![false; graph.num_edges()];
    let mut exclude = vec![false; graph.num_edges()];
    let mut path = Vec::new();
    
    // For each edge e find all cycles containing that edge and that do not
    // contain an edge with a smaller index than e:
    for e in graph.edges_iter() {
        // Add e to path:
        path.push(e.get_index());
        in_path[e.get_index()] = true;

        // Mark all vertices incident to e as not-satisfied, because they
        // occur an uneven number of times in path:
        let mut not_satisfied = Vec::from_iter(e.incident_vertices_iter().map(|v| v.get_index()));
        for v in not_satisfied.iter() {
            is_not_satisfied[*v] = true;
        }
        
        // Find all cycles containing that edge and that do not
        // contain an edge with a smaller index than e:
        find_cycles_from_path(graph, min_cycle_len, max_cycle_len,
                              &mut not_satisfied, &mut is_not_satisfied,
                              &mut path, &mut in_path, &mut exclude, &mut cycles,
                              &mut helper);
        
        // Clear the path and reset the other 'bookkeeping' variables,
        // so we can reuse them at the next iteration:
        path.clear();
        in_path[e.get_index()] = false;
        for v in not_satisfied.iter() {
            is_not_satisfied[*v] = false;
        }
    }

    cycles
}


/// Finds all simple cycles in an undirectional hypergraph `graph` of which the
/// number of edges is in the interval (`min_cycle_len` ... `max_cycle_len`),
/// that contain the edges in `path`, do not contain the edges that have been
/// marked as visited in `visited` and of which the first edge index in `path`
/// is the smallest edge index in the cycle.
/// 
/// For a definition of simple cycles, see the documentation of [find_cycles].
/// Please note the preconditions for this function in the descriptions of the
/// arguments.
/// 
/// # Arguments
/// * `graph` - The input graph.
/// * `min_cycle_len` - The minimum length of simple cycle.
/// * `max_cycle_len` - The maximum length of simple cycle.
/// * `not_satisfied` - A vector that contains the incides of all vertices of
///                     `graph` which are incident to an odd number of edges in
///                     `path`.
/// * `is_not_satisfied` - A vector of length `graph.num_vertices()` where for each
///                        vertex index `v`, `is_not_satistfied[v]` indicates 
///                        whether the vertex with index `v` is incident to an odd
///                        number of edges in `path`.
/// * `path` - A vector that contains indices of edges of `graph`. The function will
///            search for simple cycles that contain these edges.
/// * `in_path` - A vector of length `graph.num_edges()` where for each edge index
///               `e`, `is_path[e]` is true iff `e` is in `path`.
/// * `exclude` - A vector of length `graph.num_edges()` where for each edge index
///               `e`, `cycle[e] == true` indicates that the corresponding edge must
///               not be part of a found simple cycle.
/// * `found_cycles` - A vector to which this function will append for each found
///                    simple cycle, a vector of indices of the edges belonging to
///                    that simple cycle.
/// * `helper` - A helper vector of length `graph.num_vertices()` of which each
///              element is `false` before and after calling this function.
fn find_cycles_from_path (graph: &SimpleUndirHyperGraph,
                          min_cycle_len: usize,
                          max_cycle_len: usize,
                          not_satisfied: &mut Vec<usize>,
                          is_not_satisfied: &mut Vec<bool>,
                          path: &mut Vec<usize>,
                          in_path: &mut Vec<bool>,
                          exclude: &mut Vec<bool>,
                          found_cycles: &mut Vec<Vec<usize>>,
                          helper: &mut Vec<bool>) { //  f(d) = d*D_e*D_v*f(d+1) = L! * (D_e*D_v)^L 
	let in_newly_satisfied = helper;
    
    if not_satisfied.len() == 0 {
        // There are no vertices that are incident to an odd number
        // of edges in path, so we have found a cycle.
        if min_cycle_len <= path.len() {
            found_cycles.push (path.clone());
        }
    } else if path.len() < max_cycle_len {
        let mut newly_satisfied = Vec::new();
        let mut newly_dissatisfied = Vec::new();
        let mut new_not_satisfied = Vec::new();
        let mut excluded_from_path = Vec::new();

        // With the following two nested loops and if-statement we iterate
        // over all edges incident to not-satisfied vertices in `path`, excluding
        // the edges in `path`, edges marked in vector `exclude` and edges that
        // have a lower index than the first edge in `path`.
        // Adding such an edge to `path` will make the corresponding not-satisfied
        // vertex satisfied because it is then incident to an even number of edges
        // in path.
		for v in not_satisfied.iter() { // { O(d*D_e) }
			for e in graph.get_vertex(*v).incident_edges_iter() { // { O(D_v) }
				if !in_path[e.get_index()] &&
                   !exclude[e.get_index()] &&
                   e.get_index() > path[0]
                {
                    // Note that edge e may be incident to multiple not-satisfied
                    // vertices of path. To avoid trying to find simple cycles containing
                    // path and e multiple times, we only do this if v is the vertex
                    // with the smallest index amongst the not-satisfied vertices
                    // incident to e. Hence the following check:
					let mut v_is_smallest_unsatisfied_of_e = true;
					for v2 in e.incident_vertices_iter() { // { O(D_e) }
						if is_not_satisfied[v2.get_index()] &&  v2.get_index() < *v {
							v_is_smallest_unsatisfied_of_e = false;
							break;
                        }
                    }
					if v_is_smallest_unsatisfied_of_e {
						path.push(e.get_index());
						in_path[e.get_index()] = true;
						
                        // Update the is_not_satisfied value of each vertex incident to e
                        // after adding e to path. Also we keep track of all the vertices
                        // that turn from satisfied into not-satisfied and visa versa,
                        // so that we can restore the original is_not_satisfied values
                        // after removing e from path again.
                        newly_dissatisfied.clear();
                        newly_satisfied.clear();
						for v2 in e.incident_vertices_iter().map(|x| x.get_index()) { // { O(D_e) }
							if is_not_satisfied[v2] {
                                is_not_satisfied[v2] = false;
                                newly_satisfied.push(v2);
                                in_newly_satisfied[v2] = true;
                            } else {
								is_not_satisfied[v2] = true;
								newly_dissatisfied.push(v2);
							}
                        }
                        // We fill new vector, new_not_satisfied, that contains all
                        // indices of vertices that are not-satisfied after adding e
                        // to path. We use this vector as the new is_satisfied vector
                        // in the recursive to this function.
                        new_not_satisfied.clear();
						for v2 in not_satisfied.iter() { // { O(d*D_e) }
							if in_newly_satisfied[*v2] {
                                // to make sure all elements of in_newly_satisfied are false before 
                                // find_cycles_from_path() is called.
								in_newly_satisfied[*v2] = false; 
                            } else {
								new_not_satisfied.push(*v2);
                            }
                        }
						for v2 in newly_dissatisfied.iter() { // { O(D_e) } 
							new_not_satisfied.push(*v2);
                        }
						
						find_cycles_from_path(graph,
                                              min_cycle_len,
                                              max_cycle_len,
                                              &mut new_not_satisfied,
                                              is_not_satisfied,
                                              path,
                                              in_path,
                                              exclude,
                                              found_cycles,
                                              in_newly_satisfied); // {O(f(d+1))}
                        
                        // Remove e from path:
						path.pop();
						in_path[e.get_index()] = false;

                        // Restore the values of is_not_satisfied after removing e from
                        // path again:
						for v2 in newly_dissatisfied.iter() { // { O(D_e) }
							is_not_satisfied[*v2] = false;
                        }
						for v2 in newly_satisfied.iter() { // { O(D_e) }
							is_not_satisfied[*v2] = true;
                        }

                        // Now after we found all simple cycles that contain path ++ [e]
                        // we want to exclude e from being added to path in future searches
                        // within the current excution of this function, because we already
                        // found these cycles as the order of the edges in path is irrelevant:
                        exclude[e.get_index()] = true;
                        // We push it to the following vector so that we can restore values
                        // in `exclude` at the end of this function execution.
                        excluded_from_path.push(e.get_index());
                    }
                }
            }
        }
        for e in excluded_from_path.iter() {
            exclude[*e] = false;
        }
    }
}

/// Given an undirectional hyper graph `graph`, this function computes a set of
/// linearly independent cycles (of any length) of `graph` that span the set of
/// all cycles of `graph`. That is, the returned set is a basis for the cycles
/// set.
/// 
/// # Returns
/// Returns a tuple `(mat_G, mat_P_transp, k_body_constraints_basis)` where:
/// * `mat_G` - Is the generator matrix *G* corresponding to the parity
///             transformation of graph `graph`.
/// * `mat_P_transp` - Is the transpose of a parity check matrix *P* for
///                    generator matrix *G*.
/// * `k_body_constraints_basis` - A vector containing for each cycle of the
///                                computed basis, a vector containing the
///                                indices of the edges, that the cycle consists
///                                of. The *i*'th element in this vector corresponds
///                                to the *i*'th column of `mat_P_transp`.
/// 
/// For a definition of the generator matrix *G* and parity check matrix *P*,
/// see the [Ender et al., "Parity Quantum Optimization: Compiler"](https://doi.org/10.22331/q-2023-03-17-950)
/// paper.
pub fn compute_k_body_constraints_basis (graph: &SimpleUndirHyperGraph) -> (DMatrix<ZMod2>, DMatrix<ZMod2>, Vec<Vec<usize>>)
{
    
    let mat_G = generator_matrix(graph);
    let mat_P_transp = null_space(&mut mat_G.clone());
    let k_body_constraints_basis = extract_parity_constraints (&mat_P_transp);

    (mat_G, mat_P_transp, k_body_constraints_basis)
}





#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn compute_k_body_constraints_basis_on_empty_graph ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_2 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![4, 5, 6]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_3 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,3],
                                                                                        vec![3,4]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_4 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,3],
                                                                                        vec![3,0]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_5 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,3],
                                                                                        vec![3,0]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }


    #[test]
    fn compute_k_body_constraints_basis_6 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,0],
                                                                                        vec![3,4],
                                                                                        vec![4,5],
                                                                                        vec![5,3]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_7 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,0],
                                                                                        vec![3,4],
                                                                                        vec![4,5],
                                                                                        vec![5,3],
                                                                                        vec![2,3]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }

    #[test]
    fn compute_k_body_constraints_basis_8 ()
    {
        let graph = SimpleUndirHyperGraph::from_edge_list(& vec![vec![0,1,2,3,4,5],
                                                                                        vec![0,1],
                                                                                        vec![1,2],
                                                                                        vec![2,3],
                                                                                        vec![3,4],
                                                                                        vec![4,5],
                                                                                        vec![5,0]]);
        let (g,p,c) = compute_k_body_constraints_basis(&graph);
        print!("{g}, {p}, {:?}", c);
    }
}