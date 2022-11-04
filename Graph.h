#ifndef DEF_GRAPH
#define DEF_GRAPH

#include <vector>


struct Edge {
  size_t     nodeId;
  Link_types link_type;

  Edge (size_t n, Link_types l): nodeId(n), link_type(l) {}
};

struct Node {
  size_t               contigId;
  size_t               contigPartId;
  std::vector < Edge > edges;

  Node (size_t c1, size_t c2): contigId(c1), contigPartId(c2) {}
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
};

#endif
