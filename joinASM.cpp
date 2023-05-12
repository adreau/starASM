#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "globals.h"
#include "contig.h"
#include "graph.h"
#include "scaffold.h"
#include "scaffolds_to_fasta.h"


static void show_usage(char *name) {

  std::cerr << "Usage: " << name << " <option(s)>\n"
    << "Options:\n"
    << "\t-h, --help                  Show this help message\n"
    << "\t-w, --window          INT   Window size for barcode consideration (default: " << Globals::window << ") \n"
    << "\t-a, --arcsCondition   INT   Condition used for connecting two contigs; values{1..8} (default: " << Globals::condition << ", lower is more strict) \n"
    << "\t-r, --nReads          INT   Min number of common barcodes to get a links (default: " << Globals::min_n_reads << ")\n"
    << "\t-b, --begRatio        FLOAT Ratio of the contig size that is considered as the beginning part (default: " << Globals::beginning_ratio << ", should be less than 0.5)\n"
    << "\t-v, --minOverlap      INT   Minimum overlap between a molecule and a contig (default: " << Globals::min_overlap << ")\n"
    << "\t-m, --maxContigDist   INT   Merge contigs if they are separated by not more that N bp (default: " << Globals::max_contig_distance << ")\n"
    << "\t-i, --molecule        FILE  Input molecule file\n"
    << "\t-p, --mapping         FILE  Log where the molecule map with respect to the contigs (optional)\n"
    << "\t-c, --joins           FILE  Contig bed file name (result of splitASM) \n"
    << "\t-f, --contigs         FILE  Input contigs in FASTA format\n"
    << "\t-s, --scaffolds       FILE  Output scaffolds file name\n"
    << "\t-g, --graph           FILE  Output gfa file name\n"
    << "\t-o, --scaffolds       FILE  Output FASTA file\n"
    << "\t-l, --fillerSize      INT   Size of the stretch of Ns between the sequences (default: 100)\n"
    << std::endl;
}

// Read/store a set of (split) contigs in a bed file.
void create_contigs(Contigs &contigs, std::unordered_map < std::string, size_t > &contig_ids){

  std::ifstream contig_file(Globals::joins_file_name.c_str());
  std::string contig_line, ctg;
  int pos_beg, pos_end;
  unsigned int n_contigs      = 0;
  unsigned int n_contig_parts = 0;

  if (! contig_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::joins_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }

  while (getline(contig_file, contig_line)){

    std::stringstream splitstream (contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;

    // BED format is 0-based on the start, and 1-based on the end
    ++pos_beg;
    // If the contig is split, the second (third, etc.) part should be appened to the contig
    if ((! contigs.empty()) && (ctg == contigs.back().name)) {
      contigs.back().addPart(pos_beg, pos_end);
      ++n_contig_parts;
    }
    else {
      auto pos = contig_ids.find(ctg);
      // We have a new contig
      if (pos == contig_ids.end()) {
        contig_ids[ctg] = contigs.size();
        contigs.emplace_back(ctg, pos_beg, pos_end);
        ++n_contigs;
        ++n_contig_parts;
      }
      // If the file is not ordered (but why?) the contig may be already stored
      else {
        contigs[pos->second].addPart(pos_beg, pos_end);
        ++n_contig_parts;
      }
    }
  }
  std::cout << n_contigs << " contigs, and " << n_contig_parts << " contig parts, seen.\n";
}

// Count the number of common barcodes between to sets
int intersectMoleculesSize(std::vector < unsigned long int > &b1, std::vector < unsigned long int > &b2){

  std::vector < unsigned long int > common_barcodes;

  set_intersection(b1.begin(), b1.end(),
          b2.begin(), b2.end(),
          std::back_inserter(common_barcodes));
  size_t s1 = b1.size();
  size_t s2 = b2.size();
  size_t sc = common_barcodes.size();
  if (sc == 0) return 0;

  switch (Globals::condition) {
      case 1:
          if ((sc >= s1 * 0.8) && (sc >= s2 * 0.8)) return sc;
          return 0;
      case 2:
          if ((sc >= s1 * 0.8) || (sc >= s2 * 0.8)) return sc;
          return 0;
      case 3:
          if ((sc >= s1 * 0.6) && (sc >= s2 * 0.6)) return sc;
          return 0;
      case 4:
          if ((sc >= s1 * 0.6) || (sc >= s2 * 0.6)) return sc;
          return 0;
      case 5:
          if ((sc >= s1 * 0.4) && (sc >= s2 * 0.4)) return sc;
          return 0;
      case 6:
          if ((sc >= s1 * 0.4) || (sc >= s2 * 0.4)) return sc;
          return 0;
      case 7:
          if ((sc >= s1 * 0.2) && (sc >= s2 * 0.2)) return sc;
          return 0;
      case 8:
          if ((sc >= s1 * 0.2) || (sc >= s2 * 0.2)) return sc;
          return 0;
  }
  std::cerr << "Error: arc condition should be between 1 and 8.\n";
  exit(EXIT_FAILURE);
  return 0;
}

// Find the number of common barcodes between contig ends
void add_molecules_to_contigs_extremites (Contigs &contigs, std::unordered_map < std::string, size_t > &contig_ids) {

  std::cerr << "Reading molecule file...\n";
  std::ifstream molecule_file (Globals::molecule_file_name.c_str());
  std::string molecule_line, barcode, ctg, prevCtg;
  unsigned long int beg_pos, end_pos, nReads, barcode_id;
  std::unordered_set < std::string > unseen_ctgs;
  std::unordered_map < std::string, unsigned long int > barcode_to_id;
  unsigned int n_barcodes_begin  = 0;
  unsigned int n_barcodes_end    = 0;
  unsigned int n_barcodes_other  = 0;
  unsigned int n_barcodes_unused = 0;
  size_t ctg_id, prevCtg_id = 0;

  if (! molecule_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::molecule_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }


  std::ofstream mapping_file;
  if (! Globals::mapping_file_name.empty()) {
    mapping_file.open(Globals::mapping_file_name);
  }

  long unsigned int n_lines;
  for (n_lines = 0; getline(molecule_file, molecule_line); ++n_lines){

    std::stringstream  splitstream(molecule_line);
    splitstream >> ctg >> beg_pos >> end_pos >> barcode >> nReads;
    auto pos = barcode_to_id.find(barcode);
    if (pos == barcode_to_id.end()) {
      barcode_id = barcode_to_id.size();
      barcode_to_id[barcode] = barcode_to_id.size();
    }
    else {
      barcode_id = pos->second;
    }

    if (ctg == prevCtg) {
      ctg_id = prevCtg_id;
    }
    else {
      auto pos = contig_ids.find(ctg);
      if (pos == contig_ids.end()) {
        unseen_ctgs.insert(ctg);
        ++n_barcodes_unused;
        continue;
      }
      ctg_id     = pos->second;
      prevCtg    = ctg;
      prevCtg_id = ctg_id;
    }

    Interval molecule_interval (beg_pos, end_pos);
    for (ContigPart &contigPart: contigs[ctg_id].contigParts) {
      if (contigPart.get_overlap(molecule_interval) >= Globals::min_overlap) {
        // This is the condition to set a barcode to the beginning of a contig part
        if ((beg_pos >= contigPart.start) &&
            (beg_pos <= contigPart.start + std::min(Globals::window, contigPart.getSize() / 2)) &&
            (end_pos <= contigPart.start + contigPart.getSize() * (1 - Globals::beginning_ratio))) {
          contigPart.add_beg_molecule(barcode_id);
          ++n_barcodes_begin;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tbegin\n";
          }
        }
        // This is the condition to set a barcode to the end of a contig part
        else if ((end_pos <= contigPart.end) &&
            (beg_pos >= contigPart.start + contigPart.getSize() * Globals::beginning_ratio) && 
            (end_pos >= contigPart.end - std::min(Globals::window, contigPart.getSize() / 2))) {
          contigPart.add_end_molecule(barcode_id);
          ++n_barcodes_end;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tend\n";
          }
        }
        else {
          contigPart.add_other_molecule(barcode_id);
          ++n_barcodes_other;
          if (! Globals::mapping_file_name.empty()) {
            mapping_file << "Barcode " << barcode << "\t" << ctg << "\t" << beg_pos << "\t" << end_pos << "\t" << nReads << " reads\t" << ctg << "\t" << contigPart.start << "\t" << contigPart.end << "\tother\n";
          }
        }
      }
      else {
        ++n_barcodes_unused;
      }
    }
    if (n_lines % 1000000 == 0) std::cout << n_lines << " lines read.\r" << std::flush;
  }
  std::cout << n_lines << " lines read.\n";
  std::cout << n_barcodes_begin << " anchored on the left part, " <<
    n_barcodes_end << " anchored on the right part, " <<
    (n_barcodes_begin + n_barcodes_end + n_barcodes_other) << " anchored in general, and " <<
    n_barcodes_unused << " unused.\n";
  std::cout << barcode_to_id.size() << " different barcodes seen.\n";

  if (! unseen_ctgs.empty()) {
    std::vector < std::string > sorted_unseen_ctgs (unseen_ctgs.begin(), unseen_ctgs.end()); 
    sort(sorted_unseen_ctgs.begin(), sorted_unseen_ctgs.end());
    std::cerr << "Warning, " << sorted_unseen_ctgs.size() << " contigs are seen in the molecules, but not in the contigs:\n";
    for (auto &unseen_ctg: sorted_unseen_ctgs) {
      std::cerr << "\t" << unseen_ctg << "\n";
    }
  }

  std::cout << "Sorting barcodes:\n";
  for (size_t i = 0; i < contigs.size(); ++i) {
    Contig &contig = contigs[i];
    contig.sort_barcodes();
    if (i % 100 == 0) std::cout << "\tContig #" << i << " / " << contigs.size() << ".\r" << std::flush;
  }
  std::cout << "\tContig #" << contigs.size() << " / " << contigs.size() << ".\n";
}

void create_nodes (Contigs &contigs, Graph &graph) {

  for (size_t contigId = 0; contigId < contigs.size(); ++contigId) {
    Contig &contig = contigs[contigId];
    for (size_t contigPartId = 0; contigPartId < contig.contigParts.size(); ++contigPartId) {
      graph.add_node(contigId, contigPartId);
    }
  }
}

void create_cis_arcs (Contigs &contigs, Graph &graph) {

  std::cout << "Reconnecting split contigs.\n";
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
      if (nodeId % 100 == 0) std::cout << "\tNode #" << nodeId << " / " << graph.nodes.size() << ".\r" << std::flush;
      ++nodeId;
      ++n_possible_edges;
    }
    ++nodeId;
  }
  std::cout << "\tNode #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
  std::cout << "\t" << n_edges << " edges added (out of " << n_possible_edges << " possible edges).\n";
}


void create_trans_arcs (Contigs &contigs, Graph &graph) {

  std::cout << "Connecting distant nodes.\n";
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
      if (nodeId1 % 100 == 0) std::cout << "\tNode #" << nodeId1 << " / " << graph.nodes.size() << ".\r" << std::flush;
      ++nodeId1;
    }
  }
  std::cout << "\tNode #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
  std::cout << "\t" << n_edges << " edges added (out of " << (2 * graph.nodes.size() * (graph.nodes.size() - 1)) << " possible edges).\n";
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

  std::ofstream graph_file (Globals::graph_file_name.c_str(), std::ofstream::out);
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
  std::cout << "Removed " << n_removed << " / " << (2 * graph.nodes.size()) << " bifurcations.\n";
}

// Find an element which has at most one neighbor
size_t find_first_scaffold_end (Graph &graph, size_t nodeIdStart, std::vector < bool > &seen_nodes) {

  for (size_t nodeId = nodeIdStart; nodeId < graph.nodes.size(); ++nodeId) {
    if ((! seen_nodes[nodeId]) && (graph.nodes[nodeId].get_n_edges() < 2)) {
      return nodeId;
    }
  }
  return graph.nodes.size();
}

// Find the list of scaffolds in the graph
void find_scaffolds (Graph &graph, Scaffolds &scaffolds) {

  size_t n_nodes = graph.nodes.size();
  std::vector < bool > seen_nodes (n_nodes, false);
  for (size_t nodeIdStart = find_first_scaffold_end(graph, 0, seen_nodes); nodeIdStart < n_nodes; nodeIdStart = find_first_scaffold_end(graph, ++nodeIdStart, seen_nodes)) {
    size_t prev_node_id = nodeIdStart;
    Node &node = graph.nodes[nodeIdStart];
    seen_nodes[nodeIdStart] = true;
    scaffolds.emplace_back();
    // First of a chain
    if (node.get_n_edges() > 0) {
      // Find whether we should go right, or reverse the first contig
      Edge &edge = node.get_first_edge();
      bool is_forward = ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE));
      scaffolds.back().emplace_back(graph.nodes[nodeIdStart], is_forward);
      while (true) {
        size_t nodeId = edge.nodeId;
        Node &node = graph.nodes[nodeId];
        // If the edge is forward, the next edge is forward iff the link starts with E.
        // If the edge is reverse, the next edge is reverse iff the link starts with B.
        is_forward = ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE));
        scaffolds.back().emplace_back(node, is_forward);
        seen_nodes[nodeId] = true;
        // We reached the end of the chain
        if (node.get_n_edges() == 1) {
          break;
        }
        edge = node.get_other_edge(prev_node_id, reverse_link_type[edge.link_type]);
        // This only should happen if the scaffold is circular
        if (seen_nodes[edge.nodeId]) {
          break;
        }
        prev_node_id = nodeId;
      }
    }
    // Lone node
    else {
      scaffolds.back().emplace_back(graph.nodes[nodeIdStart], true);
    }
  }
}


Contig &get_contig (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId];
}

ContigPart &get_contig_part (Contigs &contigs, NodeId &nodeId) {
  return contigs[nodeId.contigId].contigParts[nodeId.contigPartId];
}

// Merge two consecutive contigs if they belong to the same scaffold, have been previously split, and are not too distant from each other
void merge_close_contigs (RefIntervalsSet &refIntervalsSet) {
  unsigned int n_merges = 0;
  for (RefIntervals &refIntervals: refIntervalsSet) {
    size_t i = 0;
    for (size_t j = 1; j < refIntervals.size(); ++j) {
      if (refIntervals[i].can_merge(refIntervals[j])) {
        refIntervals[i].merge(refIntervals[j]);
        ++n_merges;
      }
      else {
        ++i;
        if (i != j) {
          refIntervals[i] = refIntervals[j];
        }
      }
    }
    if (i < refIntervals.size()) {
      refIntervals.erase(refIntervals.begin() + i + 1, refIntervals.end());
    }
  }
  std::cout << n_merges << " contig part merges.\n";
}


void scaffold_to_intervals (Scaffolds &scaffolds, Contigs &contigs, RefIntervalsSet &refIntervalsSet) {
  for (Scaffold &scaffold: scaffolds) {
    RefIntervals refIntervals;
    for (ScaffoldPart &scaffold_part: scaffold) {
      Contig     &contig      = get_contig(contigs, scaffold_part.nodeId);
      ContigPart &contig_part = get_contig_part(contigs, scaffold_part.nodeId);
      if (contig_part.is_set()) {
        refIntervals.emplace_back(contig.name, scaffold_part.is_forward, contig_part.start, contig_part.end);
      }
    }
    refIntervalsSet.push_back(refIntervals);
  }
}


void print_scaffold (RefIntervalsSet &refIntervalsSet) {
  std::ofstream scaffold_file (Globals::scaffold_file_name, std::ofstream::out);
  if (! scaffold_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::scaffold_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }
  for (auto &refIntervals: refIntervalsSet) {
    for (auto &refInterval: refIntervals) {
      scaffold_file << refInterval << ';';
    }
    scaffold_file << '\n';
  }
  scaffold_file.close();
}


// Global values
unsigned long int Globals::min_overlap         =    300;
float             Globals::beginning_ratio     =    0.4;
unsigned long int Globals::window              = 100000;
unsigned long int Globals::max_contig_distance =  20000;
int               Globals::condition           =      1;
int               Globals::min_n_reads         =      3;
std::string       Globals::joins_file_name     =     "";
std::string       Globals::graph_file_name     =     "";
std::string       Globals::molecule_file_name  =     "";
std::string       Globals::scaffold_file_name  =     "";
std::string       Globals::mapping_file_name   =     "";
std::string       Globals::contigs_file_name   =     "";
std::string       Globals::fasta_file_name     =     "";
unsigned int      Globals::filler_size         =    100;


int main (int argc, char* argv[]) {
  if (argc < 2) {
    show_usage(argv[0]);
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } else if ((arg == "-w") || (arg == "--window")) { //maximal distance between two read pairs of the same molecule
      Globals::window = std::stoi(argv[++i]);
    } else if ((arg == "-a") || (arg == "--arcsCondition")){
      Globals::condition = std::stoi(argv[++i]);
    } else if ((arg == "-r") || (arg == "--minReads")){
      Globals::min_n_reads = std::stoi(argv[++i]);
    } else if ((arg == "-b") || (arg == "--begRatio")){
      Globals::beginning_ratio = std::stof(argv[++i]);
    } else if ((arg == "-v") || (arg == "--minOverlap")){
      Globals::min_overlap = std::stoi(argv[++i]);
    } else if ((arg == "-l") || (arg == "--fillerSize")){
      Globals::filler_size = std::atoi(argv[++i]);
    } else if ((arg == "-m") || (arg == "--maxContigDist")){
      Globals::max_contig_distance = std::stoi(argv[++i]);
    } else if ((arg == "-c") || (arg == "--joins")){
      Globals::joins_file_name = argv[++i];
    } else if ((arg == "-g") || (arg == "--graph")){
      Globals::graph_file_name = argv[++i];
    } else if ((arg == "-s") || (arg == "--scaffolds")){
      Globals::scaffold_file_name = argv[++i];
    } else if ((arg == "-p") || (arg == "--mapping")){
      Globals::mapping_file_name = argv[++i];
    } else if ((arg == "-o") || (arg == "--scaffolds")){
      Globals::fasta_file_name = argv[++i];
    } else if ((arg == "-f") || (arg == "--contigs")){
      Globals::contigs_file_name = argv[++i];
    } else if ((arg == "-i") || (arg == "--molecules")){
      Globals::molecule_file_name = argv[++i];
    } else {
      std::cerr << "Error!  Parameter '" << argv[i] << "' is not understood.\nExiting.\n";
      exit(EXIT_FAILURE);
    }
  }

  if (Globals::molecule_file_name.empty()) {
    std::cerr << "Error!  Molecule file missing.\nExiting.\n";
    exit(EXIT_FAILURE);
  }
  if (Globals::contigs_file_name.empty()) {
    std::cerr << "Error!  Input contigs FASTA file missing.\nExiting.\n";
    exit(EXIT_FAILURE);
  }
  if (Globals::fasta_file_name.empty()) {
    std::cerr << "Error!  Output scaffold FASTA file missing.\nExiting.\n";
    exit(EXIT_FAILURE);
  }
  Contigs contigs;
  std::unordered_map < std::string, size_t > contig_ids;
  create_contigs(contigs, contig_ids);
  add_molecules_to_contigs_extremites(contigs, contig_ids);
  Graph graph;
  create_nodes(contigs, graph);
  create_arcs(contigs, graph);
  write_graph(contigs, graph);
  remove_bifurcations(graph);
  Scaffolds scaffolds;
  find_scaffolds(graph, scaffolds);
  RefIntervalsSet refIntervalsSet;
  scaffold_to_intervals(scaffolds, contigs, refIntervalsSet);
  merge_close_contigs(refIntervalsSet);
  print_scaffold(refIntervalsSet);
  scaffolds_to_fasta(refIntervalsSet);

  return EXIT_SUCCESS;
}
