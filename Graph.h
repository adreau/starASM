#ifndef DEF_GRAPH
#define DEF_GRAPH

#include <vector>

// An edge points to a node id (a contig part), and describes how the contig parts are connected
struct Edge {
  size_t     nodeId;
  Link_types link_type;

  Edge (size_t n, Link_types l): nodeId(n), link_type(l) {}

  bool operator== (const Edge &edge) {
    return ((nodeId == edge.nodeId) && link_type == edge.nodeId);
  }

  bool operator!= (const Edge &edge) {
    return !(*this == edge);
  }
};

// A node id is a reference to a contig part
struct NodeId {
  size_t contigId;
  size_t contigPartId;

  NodeId (size_t c1, size_t c2): contigId(c1), contigPartId(c2) {}

  void unset () {
    contigId = unset_value;
  }

  bool is_set () {
    return (contigId == unset_value);
  }
};

// A node is a reference to a contig part, and a set of edges
struct Node: public NodeId {
  std::vector < Edge > edges;

  Node (size_t c1, size_t c2): NodeId(c1, c2) {}

  int get_n_edges () {
    int n_edges = 0;
    for (Edge &edge: edges) {
      if (edge.nodeId != unset_value) {
        ++n_edges;
      }
    }
    return n_edges;
  }

  Edge &get_first_edge () {
    for (Edge &edge: edges) {
      if (edge.nodeId != unset_value) {
        return edge;
      }
    }
    std::cerr << "Error! Trying to find the first edge whereas it seems that there is none.\n";
    exit(EXIT_FAILURE);
  }

  Edge &get_other_edge (Edge &edge1) {
    for (Edge &edge2: edges) {
      if ((edge2.nodeId != unset_value) && (edge1 != edge2)) {
        return edge2;
      }
    }
    std::cerr << "Error! Trying to find the other edge whereas node seems to contain more than 2 edges.\n";
    exit(EXIT_FAILURE);
  }
};


struct Graph {
  std::vector < Node > nodes;

  void add_node (size_t contigId, size_t contigPartId) {
    nodes.emplace_back(contigId, contigPartId);
  }

  void add_edge (size_t nodeId1, size_t nodeId2, Link_types link_type) {
    nodes[nodeId1].edges.emplace_back(nodeId2, link_type);
    nodes[nodeId2].edges.emplace_back(nodeId1, reverse_link_type[link_type]);
  }

  // Each edge is repeated: if there is an edge from node1 to node2, there should be an edge from node2 to node1
  // This fonction finds the reciprocal edge
  Edge &get_reciprocal_edge (size_t nodeId1, Edge &edge1) {
    size_t nodeId2 = edge1.nodeId;
    Link_types link_type = reverse_link_type[edge1.link_type];
    for (Edge &edge2: nodes[nodeId2].edges) {
      if ((edge2.nodeId == nodeId1) && (edge2.link_type == link_type)) {
        return edge2;
      }
    }
    std::cerr << "Error!  Cannot find reciprocal edge.\n";
    exit(EXIT_FAILURE);
  }
};

#endif
