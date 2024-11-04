//! Provides a simple struct and corresponding methods to represent undirectional
//! hyper graphs.

#[doc(hidden)]
type UndirHyperEdgeData = Vec<usize>;
#[doc(hidden)]
type UndirHyperVertexData = Vec<usize>;

/// A simple struct to represent undirectional hyper graphs. It is simple in the
/// sense that once a [SimpleUndirHyperGraph] instance has been constructed, it is
/// not manipulated afterwards. Furthermore, only functionality necessary within 
/// the scope of the assignment is provided.
/// In this graph each vertex is identified by an index in the interval [0...*n*-1]
/// where *n* is the number of vertices. Each hyper edge is identified by an index
/// in the interval [0...*m*-1], where *m* is the number of hyper edges.
pub struct SimpleUndirHyperGraph {
    #[doc(hidden)]
    edges: Vec<UndirHyperEdgeData>,
    #[doc(hidden)]
    vertices: Vec<UndirHyperVertexData>,
}

/// A reference to a vertex in a [SimpleUndirHyperGraph].
#[derive(Copy, Clone)]
pub struct UndirHyperVertexRef<'a> {
    #[doc(hidden)]
    parent_graph: &'a SimpleUndirHyperGraph,
    #[doc(hidden)]
    index: usize,
}

/// A reference to a vertex in a [SimpleUndirHyperGraph].
#[derive(Copy, Clone)]
pub struct UndirHyperEdgeRef<'a> {
    #[doc(hidden)]
    parent_graph: &'a SimpleUndirHyperGraph,
    #[doc(hidden)]
    index: usize,
}

impl SimpleUndirHyperGraph {
    /// Returns the number of vertices in the hyper graph.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Returns the number of hyper edges in the hyper graph.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Returns a reference to the vertex with index `index`.
    pub fn get_vertex(&self, index: usize) -> UndirHyperVertexRef {
        UndirHyperVertexRef {
            parent_graph: self,
            index,
        }
    }

    /// Returns a reference to the hyper edge with index `index`.
    pub fn get_edge(&self, index: usize) -> UndirHyperEdgeRef {
        UndirHyperEdgeRef {
            parent_graph: self,
            index,
        }
    }

    /// Returns an iterator over references to all vertices of the hyper graph.
    pub fn vertices_iter<'a>(&'a self) -> impl Iterator<Item = UndirHyperVertexRef<'a>> {
        (0 .. self.vertices.len())
            .map(|vertex_idx| UndirHyperVertexRef {
                parent_graph: self,
                index: vertex_idx,
            })
    }

    /// Returns an iterator over references to all hyper edges of the hyper graph.
    pub fn edges_iter<'a>(&'a self) -> impl Iterator<Item = UndirHyperEdgeRef<'a>> {
        (0 .. self.edges.len())
            .map(|edge_idx| UndirHyperEdgeRef {
                parent_graph: self,
                index: edge_idx,
            })
    }

    /// Constructs a [SimpleUndirHyperGraph] from an edge list.
    /// 
    /// # Arguments
    /// * `edge_list` - A vector that contains for each hyper edge, a vector of indices
    ///                 of the vertices incident to that hyper edge. This vector must not
    ///                 contain any duplicate hyper edges.
    /// 
    /// # Returns
    /// A [SimpleUndirHyperGraph] such that
    /// * it has *n* vertices, where *n*-1 is the highest vertex index occuring in
    ///   `edge_list`,
    /// * the indices of the vertices correspond to the vertex indices referred  to in
    ///   `edge_list`,
    /// * it has `edge_list.len()` hyper edges,
    /// * the indices of the hyper edges correspond to the indices of the hyper edges in
    ///   `edge_list`.
    pub fn from_edge_list (edge_list: &Vec<Vec<u32>>) -> Self
    {
        let mut edges = Vec::with_capacity(edge_list.len());
        let num_vertices = 1 + *edge_list.iter().map(|v| v.iter().max().unwrap_or(&0)).max().unwrap_or(&0) as usize;

        let mut vertices = vec![vec![]; num_vertices];
        for e in edge_list {
            let mut edge = Vec::with_capacity(e.len());
            for v in e {
                edge.push (*v as usize);
                vertices[*v as usize].push(edges.len());
            }
            edges.push(edge);
        }

        Self{edges, vertices}
    }
}

impl<'a> UndirHyperVertexRef<'a> {
    /// Returns an iterator over references to all hyper edges incident to the vertex.
    pub fn incident_edges_iter(&'a self) -> impl Iterator<Item = UndirHyperEdgeRef<'a>>  {
        self.parent_graph.vertices[self.index]
            .iter()
            .map(move |&edge_idx| UndirHyperEdgeRef { // TODO: not sure about move here
                parent_graph: self.parent_graph,
                index: edge_idx,
            })
    }

    pub fn get_index(&self) -> usize {
        self.index
    }
}

impl<'a> UndirHyperEdgeRef<'a> {
    /// Returns an iterator over references to all vertices incident to the hyper edge.
    pub fn incident_vertices_iter(&'a self) -> impl Iterator<Item = UndirHyperVertexRef<'a>> {
        self.parent_graph.edges[self.index]
            .iter()
            .map(move |&vertex_idx| UndirHyperVertexRef { // TODO: not sure about move here
                parent_graph: self.parent_graph,
                index: vertex_idx,
            })
    }

    pub fn get_index(&self) -> usize {
        self.index
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_vertex_iterator ()
    {
        let graph = SimpleUndirHyperGraph{edges: vec![vec![0, 1], vec![1,2], vec![2,0]],
                                          vertices: vec![vec![0,2], vec![0,1], vec![1,2], vec![]]};
        
        for v in graph.vertices_iter() {
            print!("{:?}:\n", v.get_index());
            for e in v.incident_edges_iter() {
                print!("  {:?}:", e.get_index());
                for v in e.incident_vertices_iter() {
                    print!(" {:?}", v.get_index())
                }
                print!("\n");
            }
        }
    }

    #[test]
    fn test_edge_iterator ()
    {
        let graph = SimpleUndirHyperGraph{edges: vec![vec![0, 1], vec![1,2], vec![2,0]],
                                          vertices: vec![vec![0,2], vec![0,1], vec![1,2], vec![]]};
        
        for e in graph.edges_iter() {
            print!("{:?}:\n", e.get_index());
            for v in e.incident_vertices_iter() {
                print!("  {:?}:", v.get_index());
                for e in v.incident_edges_iter() {
                    print!(" {:?}", e.get_index())
                }
                print!("\n");
            }
        }
    }
}