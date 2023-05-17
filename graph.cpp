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
  size_t nodeId = 0;
  unsigned int n_edges = 0;
  unsigned int n_possible_edges = 0;
  for (size_t contigId = 0; contigId < contigs.size(); ++contigId) {
    Contig &contig = contigs[contigId];
    for (size_t contigPartId = 0; contigPartId < contig.contigParts.size() - 1; ++contigPartId) {
      ContigPart &contigPart1 = contig.contigParts[contigPartId];
      ContigPart &contigPart2 = contig.contigParts[contigPartId + 1];
      assert(graph.nodes[nodeId].contigId       == contigId);
      assert(graph.nodes[nodeId].contigPartId   == contigPartId);
      assert(graph.nodes[nodeId+1].contigId     == contigId);
      assert(graph.nodes[nodeId+1].contigPartId == contigPartId + 1);
      if (intersectMoleculesSize(contigPart1.all_barcodes, contigPart2.all_barcodes) >= Globals::min_n_reads) {
        graph.add_edge(nodeId, nodeId + 1, Link_types::EB);
        ++n_edges;
      }
      if (nodeId % 100 == 0) std::cerr << "\tNode #" << nodeId << " / " << graph.nodes.size() << ".\r" << std::flush;
      ++nodeId;
      ++n_possible_edges;
    }
    ++nodeId;
  }
  std::cerr << "\tNode #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
  std::cerr << "\t" << n_edges << " edges added (out of " << n_possible_edges << " possible edges).\n";
}


void create_trans_arcs (Contigs &contigs, Graph &graph) {

  std::cerr << "Connecting distant nodes.\n";
  size_t nodeId1 = 0;
  size_t nodeId2 = 0;
  unsigned int n_edges = 0;
  for (size_t contigId1 = 0; contigId1 < contigs.size(); ++contigId1) {
    Contig &contig1 = contigs[contigId1];
    for (size_t contigPartId1 = 0; contigPartId1 < contig1.contigParts.size(); ++contigPartId1) {
      ContigPart &contigPart1 = contig1.contigParts[contigPartId1];
      assert(graph.nodes[nodeId1].contigId     == contigId1);
      assert(graph.nodes[nodeId1].contigPartId == contigPartId1);
      size_t contigPartId2 = contigPartId1;
      nodeId2 = nodeId1;
      for (size_t contigId2 = contigId1; contigId2 < contigs.size(); ++contigId2) {
        Contig &contig2 = contigs[contigId2];
        for (; contigPartId2 < contig2.contigParts.size(); ++contigPartId2) {
          ContigPart &contigPart2 = contig2.contigParts[contigPartId2];
          assert(graph.nodes[nodeId2].contigId     == contigId2);
          assert(graph.nodes[nodeId2].contigPartId == contigPartId2);
          // Use a strange trick to keep nodeId1 and nodeId2 synchronized
          if ((contigId1 != contigId2) || (contigPartId1 != contigPartId2)) {
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
          }
          ++nodeId2;
        }
        contigPartId2 = 0;
      }
      if (nodeId1 % 100 == 0) std::cerr << "\tNode #" << nodeId1 << " / " << graph.nodes.size() << ".\r" << std::flush;
      ++nodeId1;
    }
  }
  std::cerr << "\tNode #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
  std::cerr << "\t" << n_edges << " edges added (out of " << (2 * graph.nodes.size() * (graph.nodes.size() - 1)) << " possible edges).\n";
}

void create_arcs (Contigs &contigs, Graph &graph) {
  create_cis_arcs(contigs, graph);
  create_trans_arcs(contigs, graph);
}


void get_node_name (Contigs &contigs, Graph &graph, size_t nodeId, std::string &name) {

  Node       &node       = graph.nodes[nodeId];
  Contig     &contig     = contigs[node.contigId];
  ContigPart &contigPart = contig.contigParts[node.contigPartId];
  name = contig.name + "_" + std::to_string(contigPart.start) + "_" + std::to_string(contigPart.end);
}


void write_graph (Contigs &contigs, Graph &graph) {

  std::ofstream graph_file (Globals::graph_file_name, std::ofstream::out);
  std::string ctg1, ctg2;

  if (! graph_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::graph_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }

  graph_file << "H\tVN:Z:1.0\n";

  for (Contig &contig: contigs) {
    for (ContigPart &contigPart: contig.contigParts) {
      graph_file << "S" << "\t" << contig.name << "_" << contigPart.start << "_" << contigPart.end << "\t*\tLN:i:" << contigPart.getSize() << "\n";
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
