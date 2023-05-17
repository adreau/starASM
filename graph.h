#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

#include "contig.h"

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
    std::cerr << "Error!  Cannot find reciprocal edge for " << nodeId1 << " -> " << edge1.nodeId << "-" << edge1.link_type << ".\n";
    exit(EXIT_FAILURE);
  }
};


void create_nodes (Contigs &contigs, Graph &graph);
void create_arcs (Contigs &contigs, Graph &graph);
void write_graph (Contigs &contigs, Graph &graph);
void remove_bifurcations (Graph &graph);

#endif
