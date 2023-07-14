#ifndef NODE_H
#define NODE_H

#include <iostream>


// An edge points to a node id (a contig part), and describes how the contig parts are connected
struct Edge {
  size_t       nodeId;
  Link_types   link_type;
  double       jaccard;
  unsigned int n_reads;
  unsigned int n_min;
  unsigned int n_max;

  Edge (size_t n, Link_types l, double j, unsigned int r, unsigned int i, unsigned a): nodeId(n), link_type(l), jaccard(j), n_reads(r), n_min(i), n_max(a) {}

  void unset () {
    nodeId = unset_value;
  }

  bool is_set () {
    return (nodeId != unset_value);
  }

  bool operator== (const Edge &edge) {
    return ((nodeId == edge.nodeId) && (link_type == edge.link_type));
  }

  bool operator!= (const Edge &edge) {
    return !(*this == edge);
  }

  friend std::ostream &operator<< (std::ostream &out, Edge &e);
};

inline std::ostream &operator<< (std::ostream &out, Edge &e) {
  if (e.is_set()) {
    out << e.nodeId << "-" << e.link_type;
  }
  else {
    out << "--";
  }
  return out;
}



// A node id is a reference to a contig part
struct NodeId {
  size_t contigId;
  size_t contigPartId;

  NodeId (size_t c1, size_t c2): contigId(c1), contigPartId(c2) {}

  void unset () {
    contigId = unset_value;
  }

  bool is_set () {
    return (contigId != unset_value);
  }

  friend std::ostream &operator<< (std::ostream &out, NodeId &n);
};

inline std::ostream &operator<< (std::ostream &out, NodeId &n) {
  if (n.is_set()) {
    out << n.contigId << "-" << n.contigPartId;
  }
  else {
    out << "--";
  }
  return out;
}


// A node is a reference to a contig part, and a set of edges
struct Node: public NodeId {
  std::vector < Edge > edges;

  Node (size_t c1, size_t c2): NodeId(c1, c2) {}

  int get_n_edges () {
    int n_edges = 0;
    for (Edge &edge: edges) {
      if (edge.is_set()) {
        ++n_edges;
      }
    }
    return n_edges;
  }

  Edge &get_first_edge () {
    for (Edge &edge: edges) {
      if (edge.is_set()) {
        return edge;
      }
    }
    std::cerr << "Error! Trying to find the first edge whereas it seems that there is none.\n";
    exit(EXIT_FAILURE);
  }

  Edge &get_other_edge (size_t nodeId, Link_types link_type) {
    for (Edge &edge: edges) {
      if ((edge.nodeId != nodeId) || (edge.link_type != link_type)) {
        return edge;
      }
    }
    std::cerr << "Error! Trying to find the other edge whereas node seems to contain less than 2 edges.\n";
    exit(EXIT_FAILURE);
  }

  friend std::ostream &operator<< (std::ostream &out, Node &n);
};

inline std::ostream &operator<< (std::ostream &out, Node &n) {
  if (n.is_set()) {
    out << static_cast < NodeId& > (n);
    for (Edge &e: n.edges) {
      if (e.is_set()) {
        out << " -> " << e;
      }
    }
  }
  else {
    out << "--";
  }
  return out;
}

#endif
