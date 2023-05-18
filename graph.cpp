#include <cassert>
#include <fstream>

#include "globals.h"
#include "contig.h"
#include "graph.h"

void create_nodes (Contigs &contigs, Graph &graph) {
  for (size_t contigId = 0; contigId < contigs.size(); ++contigId) {
    Contig &contig = contigs[contigId];
    for (size_t contigPartId = 0; contigPartId < contig.contigParts.size(); ++contigPartId) {
      graph.add_node(contigId, contigPartId);
    }
  }
}

void create_cis_arcs (Contigs &contigs, Graph &graph) {
  std::cerr << "Reconnecting split contigs.\n";
  unsigned int n_edges = 0;
  unsigned int n_possible_edges = 0;
  for (size_t nodeId = 1; nodeId < graph.nodes.size(); ++nodeId) {
    Node &prevNode = graph.nodes[nodeId - 1];
    Node &nextNode = graph.nodes[nodeId];
    if (prevNode.contigId == nextNode.contigId) {
      ++n_possible_edges;
      ContigPart prevContigPart = get_contig_part(contigs, prevNode);
      ContigPart nextContigPart = get_contig_part(contigs, nextNode);
      if (intersectMoleculesSize(prevContigPart.all_barcodes, nextContigPart.all_barcodes) >= Globals::min_n_reads) {
        graph.add_edge(nodeId - 1, nodeId, Link_types::EB);
        ++n_edges;
      }
      if (nodeId % 100 == 0) std::cerr << "\tNode #" << nodeId << " / " << graph.nodes.size() << ".\r" << std::flush;
    }
  }
  std::cerr << "\tNode #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
  std::cerr << "\t" << n_edges << " edges added (out of " << n_possible_edges << " possible edges).\n";
}


void create_trans_arcs (Contigs &contigs, Graph &graph) {
  std::cerr << "Connecting distant contigs.\n";
  unsigned int edge_id       = 0;
  unsigned int n_edges       = 0;
  unsigned int total_n_edges = graph.nodes.size() * (graph.nodes.size() - 1) / 2;
  for (size_t nodeId1 = 0; nodeId1 < graph.nodes.size(); ++nodeId1) {
    Node       &node1       = graph.nodes[nodeId1];
    ContigPart &contigPart1 = get_contig_part(contigs, node1);
    for (size_t nodeId2 = nodeId1 + 1; nodeId2 < graph.nodes.size(); ++nodeId2) {
      Node       &node2       = graph.nodes[nodeId2];
      ContigPart &contigPart2 = get_contig_part(contigs, node2);
      if (intersectMoleculesSize(contigPart1.barcodes_beg, contigPart2.barcodes_beg) >= Globals::min_n_reads) {
        graph.add_edge(nodeId1, nodeId2, Link_types::BB);
        ++n_edges;
      }
      if (intersectMoleculesSize(contigPart1.barcodes_beg, contigPart2.barcodes_end) >= Globals::min_n_reads) {
        graph.add_edge(nodeId1, nodeId2, Link_types::BE);
        ++n_edges;
      }
      if (intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_beg) >= Globals::min_n_reads) {
        graph.add_edge(nodeId1, nodeId2, Link_types::EB);
        ++n_edges;
      }
      if (intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_end) >= Globals::min_n_reads) {
        graph.add_edge(nodeId1, nodeId2, Link_types::EE);
        ++n_edges;
      }
      if (edge_id % 1000 == 0) std::cerr << "\tEdge #" << edge_id << " / " << total_n_edges << ".\r" << std::flush;
      ++edge_id;
    }
  }
  std::cerr << "\tEdge #" << total_n_edges << " / " << total_n_edges << ".\n";
  std::cerr << "\t" << n_edges << " edges added (out of " << total_n_edges << " possible edges).\n";
}

void create_arcs (Contigs &contigs, Graph &graph) {
  create_cis_arcs(contigs, graph);
  create_trans_arcs(contigs, graph);
}


void get_node_name (Contigs &contigs, Graph &graph, size_t nodeId, std::string &name) {
  Node       &node       = graph.nodes[nodeId];
  Contig     &contig     = contigs[node.contigId];
  ContigPart &contigPart = contig.contigParts[node.contigPartId];
  name = Globals::chrs[node.contigId] + "_" + std::to_string(contigPart.start) + "_" + std::to_string(contigPart.end);
}


void write_graph (Contigs &contigs, Graph &graph) {
  std::ofstream graph_file (Globals::graph_file_name, std::ofstream::out);
  std::string ctg1, ctg2;
  if (! graph_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::graph_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }
  graph_file << "H\tVN:Z:1.0\n";
  for (unsigned int i = 0; i < contigs.size(); ++i) {
    Contig &contig = contigs[i];
    for (ContigPart &contigPart: contig.contigParts) {
      graph_file << "S" << "\t" << Globals::chrs[i] << "_" << contigPart.start << "_" << contigPart.end << "\t*\tLN:i:" << contigPart.getSize() << "\n";
    }
  }
  for (size_t nodeId1 = 0; nodeId1 < graph.nodes.size(); ++nodeId1) {
    Node &node = graph.nodes[nodeId1];
    get_node_name(contigs, graph, nodeId1, ctg1);
    for (Edge &edge: node.edges) {
      size_t nodeId2 = edge.nodeId;
      if (nodeId1 < nodeId2) {
        get_node_name(contigs, graph, nodeId2, ctg2);
        switch (edge.link_type) {
          case Link_types::BB:
            graph_file << "L\t" << ctg2 << "\t-\t" << ctg1 << "\t+\n";
            break;
          case Link_types::BE:
            graph_file << "L\t" << ctg2 << "\t+\t" << ctg1 << "\t+\n";
            break;
          case Link_types::EB:
            graph_file << "L\t" << ctg1 << "\t+\t" << ctg2 << "\t+\n";
            break;
          case Link_types::EE:
            graph_file << "L\t" << ctg1 << "\t+\t" << ctg2 << "\t-\n";
            break;
        }
      }
    }
  }
  graph_file.close();
}


// Remove arcs when bifurcations are observed
// Arcs are actually not removed, but replaced by -1
void remove_bifurcations (Graph &graph) {

  unsigned int n_removed = 0;
  for (size_t nodeId = 0; nodeId < graph.nodes.size(); ++nodeId) {
    Node &node  = graph.nodes[nodeId];
    int n_begin = 0;
    int n_end   = 0;
    for (Edge &edge: node.edges) {
      if ((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::BE)) {
        ++n_begin;
      }
      else {
        ++n_end;
      }
      // Found a bifurcation at the beginning of the contig
      if (n_begin > 1) {
        for (Edge &edge: node.edges) {
          if (edge.is_set()) {
            if ((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::BE)) {
              // Find the corresponding edge in the other node (each edge is present twice)
              graph.get_reciprocal_edge(nodeId, edge).unset();
              // Discard this edge (do it after having found the reciprocal edge!)
              edge.unset();
            }
          }
        }
        ++n_removed;
      }
      // Found a bifurcation at the end of the contig
      if (n_end > 1) {
        for (Edge &edge: node.edges) {
          if (edge.is_set()) {
            if ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE)) {
              graph.get_reciprocal_edge(nodeId, edge).unset();
              edge.unset();
            }
          }
        }
        ++n_removed;
      }
    }
  }
  std::cerr << "Removed " << n_removed << " / " << (2 * graph.nodes.size()) << " bifurcations.\n";
}
