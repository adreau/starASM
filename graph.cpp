#include <cassert>
#include <fstream>

#include "constants.h"
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

// Transfer the content of the left contig part to the right contig part
void merge_contig_parts (ContigPart &left_contig_part, ContigPart &right_contig_part) {
  right_contig_part.barcodes_beg = std::move(left_contig_part.barcodes_beg);
  right_contig_part.start        = left_contig_part.start;
  right_contig_part.all_barcodes.insert(right_contig_part.all_barcodes.end(), left_contig_part.all_barcodes.begin(), left_contig_part.all_barcodes.end());
  std::sort(right_contig_part.all_barcodes.begin(), right_contig_part.all_barcodes.end());
  right_contig_part.all_barcodes.erase(std::unique(right_contig_part.all_barcodes.begin(), right_contig_part.all_barcodes.end()), right_contig_part.all_barcodes.end());
}

void create_cis_arcs (Contigs &contigs, Graph &graph) {
  std::ofstream cis_link_file;
  if (! Globals::cis_link_file_name.empty()) cis_link_file.open(Globals::cis_link_file_name, std::ofstream::out);
  std::cerr << TAB << "Reconnecting split contigs...\n";
  unsigned int n_edges          = 0;
  unsigned int n_possible_edges = 0;
  double jaccard;
  for (size_t nodeId = 1; nodeId < graph.nodes.size(); ++nodeId) {
    Node &prevNode = graph.nodes[nodeId - 1];
    Node &nextNode = graph.nodes[nodeId];
    if (prevNode.contigId == nextNode.contigId) {
      ++n_possible_edges;
      ContigPart &prevContigPart = get_contig_part(contigs, prevNode);
      ContigPart &nextContigPart = get_contig_part(contigs, nextNode);
      unsigned int n_inter;
      intersectMoleculesSize(prevContigPart.all_barcodes, nextContigPart.all_barcodes, n_inter, jaccard);
      bool merged = (n_inter >= Globals::min_n_reads);
      if (! Globals::cis_link_file_name.empty()) {
        cis_link_file << Globals::chrs[prevNode.contigId] << TAB << prevContigPart.start << TAB << prevContigPart.end << TAB <<
                         Globals::chrs[nextNode.contigId] << TAB << nextContigPart.start << TAB << nextContigPart.end << TAB <<
                         prevContigPart.all_barcodes.size() << TAB << nextContigPart.all_barcodes.size() << TAB << n_inter << TAB << jaccard << TAB << merged << "\n";
      }
      if (merged) {
        merge_contig_parts(prevContigPart, nextContigPart);
        prevNode.unset();
        //graph.add_edge(nodeId - 1, nodeId, Link_types::EB);
        ++n_edges;
      }
      if (nodeId % 100 == 0) std::cerr << TAB << TAB << "Split contig " << nodeId << "/" << graph.nodes.size() << "\r" << std::flush;
    }
  }
  if (! Globals::cis_link_file_name.empty()) cis_link_file.close();
  std::cerr << TAB << TAB << "Split contig " << graph.nodes.size() << "/" << graph.nodes.size() << "\n";
  std::cerr << TAB << TAB << n_edges << " edges added (out of " << n_possible_edges << " possible edges)\n";
}


void create_trans_arcs (Contigs &contigs, Graph &graph) {
  std::ofstream trans_link_file;
  if (! Globals::trans_link_file_name.empty()) trans_link_file.open(Globals::trans_link_file_name, std::ofstream::out);
  unsigned int edge_id = 0;
  unsigned int n_nodes = 0;
  unsigned int n_edges = 0;
  unsigned int total_n_edges;
  bool         kept;
  for (auto &node: graph.nodes) {
    if (node.is_set()) {
      ++n_nodes;
	}
  }
  std::cerr << TAB << "Connecting distant contigs from " << n_nodes << " contig parts...\n";
  total_n_edges = (n_nodes) * (n_nodes - 1) / 2;
  for (size_t nodeId1 = 0; nodeId1 < graph.nodes.size(); ++nodeId1) {
    Node &node1 = graph.nodes[nodeId1];
    if (node1.is_set()) {
      ContigPart &contigPart1 = get_contig_part(contigs, node1);
      for (size_t nodeId2 = nodeId1 + 1; nodeId2 < graph.nodes.size(); ++nodeId2) {
        Node &node2 = graph.nodes[nodeId2];
        if (node2.is_set()) {
          ContigPart &contigPart2 = get_contig_part(contigs, node2);
          unsigned int n_inter;
          double jaccard;
		  unsigned int n1 = contigPart1.barcodes_beg.size();
		  unsigned int n2 = contigPart2.barcodes_beg.size();
		  unsigned int n_min = std::min(n1, n2);
		  unsigned int n_max = std::max(n1, n2);
          intersectMoleculesSize(contigPart1.barcodes_beg, contigPart2.barcodes_beg, n_inter, jaccard);
          kept = false;
          if ((n_inter >= Globals::min_n_reads) && (jaccard >= Globals::jaccard)) {
            graph.add_edge(nodeId1, nodeId2, Link_types::BB, jaccard, n_inter, n_min, n_max);
            ++n_edges;
            kept = true;
          }
          if (! Globals::trans_link_file_name.empty()) {
            trans_link_file << Globals::chrs[node1.contigId] << TAB << contigPart1.start << TAB << contigPart1.end << TAB << 'B' << TAB <<
                               Globals::chrs[node2.contigId] << TAB << contigPart2.start << TAB << contigPart2.end << TAB << 'B' << TAB <<
                               contigPart1.barcodes_beg.size() << TAB << contigPart2.barcodes_beg.size() << TAB << n_inter << TAB << jaccard << TAB << kept << "\n";
          }
		  n1 = contigPart1.barcodes_beg.size();
		  n2 = contigPart2.barcodes_end.size();
		  n_min = std::min(n1, n2);
		  n_max = std::max(n1, n2);
          intersectMoleculesSize(contigPart1.barcodes_beg, contigPart2.barcodes_end, n_inter, jaccard);
          kept = false;
          if ((n_inter >= Globals::min_n_reads) && (jaccard >= Globals::jaccard)) {
            graph.add_edge(nodeId1, nodeId2, Link_types::BE, jaccard, n_inter, n_min, n_max);
            ++n_edges;
            kept = true;
          }
          if (! Globals::trans_link_file_name.empty()) {
            trans_link_file << Globals::chrs[node1.contigId] << TAB << contigPart1.start << TAB << contigPart1.end << TAB << 'B' << TAB <<
                               Globals::chrs[node2.contigId] << TAB << contigPart2.start << TAB << contigPart2.end << TAB << 'E' << TAB <<
                               contigPart1.barcodes_beg.size() << TAB << contigPart2.barcodes_end.size() << TAB << n_inter << TAB << jaccard << TAB << kept << "\n";
          }
		  n1 = contigPart1.barcodes_end.size();
		  n2 = contigPart2.barcodes_beg.size();
		  n_min = std::min(n1, n2);
		  n_max = std::max(n1, n2);
          intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_beg, n_inter, jaccard);
          kept = false;
          if ((n_inter >= Globals::min_n_reads) && (jaccard >= Globals::jaccard)) {
            graph.add_edge(nodeId1, nodeId2, Link_types::EB, jaccard, n_inter, n_min, n_max);
            ++n_edges;
            kept = true;
          }
          if (! Globals::trans_link_file_name.empty()) {
            trans_link_file << Globals::chrs[node1.contigId] << TAB << contigPart1.start << TAB << contigPart1.end << TAB << 'E' << TAB <<
                               Globals::chrs[node2.contigId] << TAB << contigPart2.start << TAB << contigPart2.end << TAB << 'B' << TAB <<
                               contigPart1.barcodes_end.size() << TAB << contigPart2.barcodes_beg.size() << TAB << n_inter << TAB << jaccard << TAB << kept << "\n";
          }
		  n1 = contigPart1.barcodes_end.size();
		  n2 = contigPart2.barcodes_end.size();
		  n_min = std::min(n1, n2);
		  n_max = std::max(n1, n2);
          intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_end, n_inter, jaccard);
          kept = false;
          if ((n_inter >= Globals::min_n_reads) && (jaccard >= Globals::jaccard)) {
            graph.add_edge(nodeId1, nodeId2, Link_types::EE, jaccard, n_inter, n_min, n_max);
            ++n_edges;
            kept = true;
          }
          if (! Globals::trans_link_file_name.empty()) {
            trans_link_file << Globals::chrs[node1.contigId] << TAB << contigPart1.start << TAB << contigPart1.end << TAB << 'E' << TAB <<
                               Globals::chrs[node2.contigId] << TAB << contigPart2.start << TAB << contigPart2.end << TAB << 'E' << TAB <<
                               contigPart1.barcodes_end.size() << TAB << contigPart2.barcodes_end.size() << TAB << n_inter << TAB << jaccard << TAB << kept << "\n";
          }
          if (edge_id % 1000000 == 0) std::cerr << TAB << TAB << "Inspecting edge " << edge_id << "/" << total_n_edges << "\r" << std::flush;
          ++edge_id;
        }
      }
    }
  }
  if (! Globals::trans_link_file_name.empty()) trans_link_file.close();
  std::cerr << TAB << TAB << "Inspecting edge " << total_n_edges << "/" << total_n_edges << "\n";
  std::cerr << TAB << TAB << n_edges << " edges added (out of " << total_n_edges << " possible edges)\n";
}

void create_arcs (Contigs &contigs, Graph &graph) {
  create_cis_arcs(contigs, graph);
  create_trans_arcs(contigs, graph);
}


void get_node_name (Contigs &contigs, Graph &graph, size_t nodeId, std::string &name) {
  Node       &node       = graph.nodes[nodeId];
  Contig     &contig     = contigs[node.contigId];
  ContigPart &contigPart = contig.contigParts[node.contigPartId];
  name = Globals::chrs[node.contigId] + "_" + std::to_string(contigPart.start + 1) + "_" + std::to_string(contigPart.end);
}


void write_graph (Contigs &contigs, Graph &graph) {
  std::ofstream graph_file (Globals::graph_file_name, std::ofstream::out);
  std::string ctg, ctg1, ctg2;
  if (! graph_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::graph_file_name << "'.\n";
      exit(EXIT_FAILURE);
  }
  graph_file << "H\tVN:Z:1.0\n";
  for (size_t nodeId = 0; nodeId < graph.nodes.size(); ++nodeId) {
    Node &node = graph.nodes[nodeId];
	if (node.is_set()) {
      get_node_name(contigs, graph, nodeId, ctg);
      graph_file << "S" << "\t" << ctg << "\t*\tLN:i:" << contigs[node.contigId].contigParts[node.contigPartId].getSize() << "\n";
    }
  }
  for (size_t nodeId1 = 0; nodeId1 < graph.nodes.size(); ++nodeId1) {
    Node &node = graph.nodes[nodeId1];
	if (node.is_set()) {
      get_node_name(contigs, graph, nodeId1, ctg1);
      for (Edge &edge: node.edges) {
        size_t nodeId2 = edge.nodeId;
        if (nodeId1 < nodeId2) {
          get_node_name(contigs, graph, nodeId2, ctg2);
          switch (edge.link_type) {
            case Link_types::BB:
              graph_file << "L\t" << ctg2 << "\t-\t" << ctg1 << "\t+\t*\tJA:f:" << edge.jaccard << " RC:i:" << edge.n_reads << " NI:i:" << edge.n_min << " NA:i:" << edge.n_max << "\n";
              break;
            case Link_types::BE:
              graph_file << "L\t" << ctg2 << "\t+\t" << ctg1 << "\t+\t*\tJA:f:" << edge.jaccard << " RC:i:" << edge.n_reads << " NI:i:" << edge.n_min << " NA:i:" << edge.n_max << "\n";
              break;
            case Link_types::EB:
              graph_file << "L\t" << ctg1 << "\t+\t" << ctg2 << "\t+\t*\tJA:f:" << edge.jaccard << " RC:i:" << edge.n_reads << " NI:i:" << edge.n_min << " NA:i:" << edge.n_max << "\n";
              break;
            case Link_types::EE:
              graph_file << "L\t" << ctg1 << "\t+\t" << ctg2 << "\t-\t*\tJA:f:" << edge.jaccard << " RC:i:" << edge.n_reads << " NI:i:" << edge.n_min << " NA:i:" << edge.n_max << "\n";
              break;
          }
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
  unsigned int n_nodes   = 0;
  for (size_t nodeId = 0; nodeId < graph.nodes.size(); ++nodeId) {
    Node &node = graph.nodes[nodeId];
	if (node.is_set()) {
      int n_begin = 0;
      int n_end   = 0;
      ++n_nodes;
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
  }
  std::cerr << TAB << "Removed " << n_removed << "/" << (2 * n_nodes) << " bifurcations\n";
}
