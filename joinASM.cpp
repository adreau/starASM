#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "Globals.h"
#include "Contig.h"
#include "Graph.h"
#include "Scaffold.h"


static void show_usage(char *name) {

  std::cerr << "Usage: " << name << " <option(s)> MOLECULE FILE \n"
    << "Options:\n"
    << "\t-h, --help                  Show this help message\n"
    << "\t-w, --window          INT   Window size for barcode consideration (default: " << Globals::window << ") \n"
    << "\t-a, --arcsCondition   INT   Condition used for connecting two contigs; values{1..8} (default: " << Globals::condition << ", lower is more strict) \n"
    << "\t-r, --nReads          INT   Min number of common barcodes to get a links (default: " << Globals::min_n_reads << ")\n"
    << "\t-b, --begRatio        FLOAT Ratio of the contig size that is considered as the beginning part (default: " << Globals::beginning_ratio << ", should be less than 0.5)\n"
    << "\t-p, --pairReadsLength INT   ??? (default: " << Globals::pair_reads_length << ")\n"
    << "\t-m, --maxContigDist   INT   merge contigs if they are separated by not more that N bp (default: " << Globals::max_contig_distance << ")\n"
    << "\t-c, --contigs         FILE  Contig bed file name (result of splitASM) \n"
    << "\t-s, --scaffolds       FILE  Output scaffolds file name \n"
    << "\t-g, --graph           FILE  Output gfa file name \n"
    << std::endl;
}

// Read/store a set of (split) contigs in a bed file.
void create_contigs(Contigs &contigs, std::unordered_map < std::string, size_t > &contig_ids){

  std::ifstream contig_file(Globals::contig_file_name.c_str());
  std::string contig_line, ctg;
  int pos_beg, pos_end;

  if (! contig_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::contig_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }

  while (getline(contig_file, contig_line)){

    std::stringstream splitstream (contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;

    // If the contig is split, the second (third, etc.) part should be appened to the contig
    if ((! contigs.empty()) && (ctg == contigs.back().name)) {
      contigs.back().addPart(pos_beg, pos_end);
    }
    else {
      auto pos = contig_ids.find(ctg);
      // We have a new contig
      if (pos == contig_ids.end()) {
        contig_ids[ctg] = contigs.size();
        contigs.emplace_back(ctg, pos_beg, pos_end);
      }
      // If the file is not ordered (but why?) the contig may be already stored
      else {
        contigs[pos->second].addPart(pos_beg, pos_end);
      }
    }
  }
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
void add_molecules_to_contigs_extremites(Contigs &contigs, std::unordered_map < std::string, size_t > &contig_ids){

  std::ifstream molecule_file (Globals::molecule_file_name.c_str());
  std::string molecule_line, barcode, ctg, prevCtg;
  unsigned long int beg_pos, end_pos, nReads, barcode_id;
  std::unordered_set < std::string > unseen_ctgs;
  std::unordered_map < std::string, unsigned long int > barcode_to_id;
  size_t ctg_id, prevCtg_id = 0;

  if (! molecule_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::molecule_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
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
        continue;
      }
      ctg_id     = pos->second;
      prevCtg    = ctg;
      prevCtg_id = ctg_id;
    }

    for (ContigPart &contigPart: contigs[ctg_id].contigParts) {
      // This is the condition to set a barcode to the beginning of a contig part
      if ((beg_pos >= contigPart.begin) && (beg_pos <= contigPart.begin + std::min(Globals::window, contigPart.getSize() / 2)) &&
          (end_pos <= contigPart.begin + contigPart.getSize() * (1 - Globals::beginning_ratio))) {
        contigPart.add_beg_molecule(barcode_id);
      }
      // This is the condition to set a barcode to the end of a contig part
      else if ((beg_pos >= contigPart.begin + contigPart.getSize() * Globals::beginning_ratio) && (beg_pos <= contigPart.end - Globals::pair_reads_length) &&
          (end_pos >= contigPart.end - std::min(Globals::window, contigPart.getSize() / 2))) {
        contigPart.add_end_molecule(barcode_id);
      }
    }
    if (n_lines % 10000000 == 0) std::cout << n_lines << " lines read.\r" << std::flush;
  }
  std::cout << n_lines << " lines read.\n";

  if (! unseen_ctgs.empty()) {
    std::cerr << "Warning, several contigs are seen in the molecules, but not in the contigs:\n";
    for (auto &unseen_ctg: unseen_ctgs) {
      std::cerr << "\t" << unseen_ctg << "\n";
    }
  }

  for (Contig &contig: contigs) {
    contig.sort_barcodes();
  }
}

void create_nodes (Contigs &contigs, Graph &graph) {

  for (size_t contigId = 0; contigId < contigs.size(); ++contigId) {
    Contig &contig = contigs[contigId];
    for (size_t contigPartId = 0; contigPartId < contig.contigParts.size(); ++contigPartId) {
      graph.add_node(contigId, contigPartId);
    }
  }
}


void create_arcs (Contigs &contigs, Graph &graph) {

  size_t nodeId1 = 0;
  size_t nodeId2 = 0;
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
            }
            if (intersectMoleculesSize(contigPart1.barcodes_beg, contigPart2.barcodes_end) >= Globals::min_n_reads) {
              graph.add_edge(nodeId1, nodeId2, Link_types::BE);
            }
            if (intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_beg) >= Globals::min_n_reads) {
              graph.add_edge(nodeId1, nodeId2, Link_types::EB);
            }
            if (intersectMoleculesSize(contigPart1.barcodes_end, contigPart2.barcodes_end) >= Globals::min_n_reads) {
              graph.add_edge(nodeId1, nodeId2, Link_types::EE);
            }
          }
          ++nodeId2;
        }
        contigPartId2 = 0;
      }
      if (nodeId1 % 100 == 0) std::cout << "Node #" << nodeId1 << " / " << graph.nodes.size() << ".\r" << std::flush;
      ++nodeId1;
    }
  }
  std::cout << "Node #" << graph.nodes.size() << " / " << graph.nodes.size() << ".\n";
}


void get_node_name (Contigs &contigs, Graph &graph, size_t nodeId, std::string &name) {

  Node       &node       = graph.nodes[nodeId];
  Contig     &contig     = contigs[node.contigId];
  ContigPart &contigPart = contig.contigParts[node.contigPartId];
  name = contig.name + "_" + std::to_string(contigPart.begin) + "_" + std::to_string(contigPart.end);
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
      graph_file << "S" << "\t" << contig.name << "_" << contigPart.begin << "_" << contigPart.end << "\t*\tLN:i:" << contigPart.getSize() << "\n";
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

  for (size_t nodeId = 0; nodeId < graph.nodes.size(); ++nodeId) {
    Node &node  = graph.nodes[nodeId];
    int n_begin = 0;
    int n_end   = 0;
    for (Edge &edge: node.edges) {
      if (edge.nodeId != unset_value) {
        if ((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::BE)) {
          ++n_begin;
        }
        else {
          ++n_end;
        }
      }
      // Found a bifurcation at the beginning of the contig
      if (n_begin > 1) {
        for (Edge &edge: node.edges) {
          if ((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::BE)) {
            // Discard this edge
            edge.nodeId = unset_value;
            // Find the corresponding edge in the other node (each edge is present twice)
            graph.get_reciprocal_edge(nodeId, edge).nodeId = unset_value;
          }
        }
      }
      // Found a bifurcation at the end of the contig
      if (n_end > 1) {
        for (Edge &edge: node.edges) {
          if ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE)) {
            edge.nodeId = unset_value;
            graph.get_reciprocal_edge(nodeId, edge).nodeId = unset_value;
          }
        }
      }
    }
  }
}

// Find an element which has at most one neighbor
size_t find_first_scaffold_end (Graph &graph, size_t nodeIdStart, std::vector < bool > &seen_nodes) {

  for (size_t nodeId = nodeIdStart; nodeId < graph.nodes.size(); ++nodeId) {
    if (! seen_nodes[nodeId]) {
      Node        &node    = graph.nodes[nodeId];
      unsigned int n_nodes = 0;
      for (Edge &edge: node.edges) {
        if (edge.nodeId != unset_value) {
          ++n_nodes;
        }
        if (n_nodes < 2) {
          return nodeId;
        }
      }
    }
  }
  return graph.nodes.size();
}

// Find the list of scaffolds in the graph
void find_scaffolds (Graph &graph, Scaffolds &scaffolds) {

  size_t n_nodes = graph.nodes.size();
  std::vector < bool > seen_nodes (n_nodes, false);
  for (size_t nodeIdStart = 0; nodeIdStart < n_nodes; nodeIdStart = find_first_scaffold_end(graph, nodeIdStart, seen_nodes)) {
    Node &node = graph.nodes[nodeIdStart];
    seen_nodes[nodeIdStart] = true;
    scaffolds.emplace_back();
    // First of a chain
    if (node.get_n_edges() > 1) {
      // Find whether we should go right, or reverse the first contig
      Edge &edge = node.get_first_edge();
      bool is_forward = ((edge.link_type == Link_types::EB) || (edge.link_type == Link_types::EE));
      scaffolds.back().emplace_back(graph.nodes[nodeIdStart], is_forward);
      while (true) {
        Node &nextNode = graph.nodes[edge.nodeId];
        is_forward = (((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::EB)) == is_forward);
        scaffolds.back().emplace_back(graph.nodes[edge.nodeId], is_forward);
        seen_nodes[edge.nodeId] = true;
        if (nextNode.get_n_edges() == 1) {
          break;
        }
        edge = node.get_other_edge(edge);
        // This only should happen if the scaffold is circular
        if (seen_nodes[edge.nodeId]) {
          break;
        }
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

// Merge two consecutive contigs if the belong to the same scaffold, have been previously split, and are not too distant from each other
void merge_close_contigs (Graph &graph, Contigs &contigs, Scaffolds &scaffolds) {

  for (Scaffold &scaffold: scaffolds) {
    ScaffoldPart &prev_scaffold_part = scaffold.front();
    ContigPart   &prev_contig_part   = get_contig_part(contigs, prev_scaffold_part.nodeId);
    for (size_t scaffold_part_id = 1; scaffold_part_id < scaffold.size(); ++scaffold_part_id) {
      ScaffoldPart &next_scaffold_part = scaffold[scaffold_part_id];
      ContigPart   &next_contig_part   = get_contig_part(contigs, next_scaffold_part.nodeId);
      if ((prev_scaffold_part.is_forward == next_scaffold_part.is_forward) && (prev_contig_part.get_distance(next_contig_part) <= Globals::max_contig_distance)) {
        next_contig_part.merge(prev_contig_part);
        prev_contig_part.unset();
      }
      next_scaffold_part = prev_scaffold_part;
      next_contig_part   = prev_contig_part;
    }
  }
}


void print_scaffold (Scaffolds &scaffolds, Contigs &contigs) {

  std::ofstream scaffold_file (Globals::scaffold_file_name.c_str(), std::ofstream::out);

  if (! scaffold_file.is_open()){
      std::cerr << "Error!  Cannot open file '" << Globals::scaffold_file_name << "'" << std::endl;
      exit(EXIT_FAILURE);
  }

  for (Scaffold &scaffold: scaffolds) {
    for (ScaffoldPart &scaffold_part: scaffold) {
      Contig     &contig      = get_contig(contigs, scaffold_part.nodeId);
      ContigPart &contig_part = get_contig_part(contigs, scaffold_part.nodeId);
      if (contig_part.is_set()) {
        scaffold_file << (scaffold_part.is_forward? '+': '-') << contig.name << ' ' << contig_part.begin << ' ' << contig_part.end << ';';
      }
    }
    scaffold_file << '\n';
  }

  scaffold_file.close();
}



/*

void get_canonical_ctg(std::pair<std::string,std::string> &arc, std::string &ctg1, std::string &sign1, std::string &ctg2, std::string &sign2){

  std::string tmp_sign1, tmp_sign2;
  int pos = arc.first.find(":");
  ctg1 = arc.first.substr(0,pos);
  tmp_sign1 = arc.first.substr(pos);

  pos = arc.second.find(":");
  ctg2 = arc.second.substr(0,pos);
  tmp_sign2 = arc.second.substr(pos);

  if(tmp_sign1.compare(":beg")==0)
    sign1="-";
  else sign1="+";

  if(tmp_sign2.compare(":beg")==0)
    sign2="+";
  else sign2="-";

}

char opposing_sign(char sign){

  if(sign == '+') return '-';
  else return '+';

}

void print_scaffold( std::vector<std::string> &scaffold, std::vector<char> &signs, std::ofstream &scaffoldFile){

  std::cout << "Scaffold size: "<< scaffold.size()<<std::endl;
  if(signs[0] == '-') {
    for(int i=0; i<scaffold.size();i++)
      scaffoldFile << opposing_sign(signs[i]) <<scaffold[i]<<";";
  }else {
    for(int i=0; i<scaffold.size();i++)
      scaffoldFile << signs[i] <<scaffold[i]<<";";
  }
  scaffoldFile <<"\n";
}

void getNotVistedNeighbors(int ctg_beg_int , UndirectedGraph &undigraph, std::vector<bool> &visited, std::vector < int > &neighbors){

  typename graph_traits<UndirectedGraph>::adjacency_iterator vi, vi_end;

  for (boost::tie(vi, vi_end) = adjacent_vertices(ctg_beg_int, undigraph); vi != vi_end; ++vi){

    std::cout << "Neighbors all: "<<*vi<<std::endl;
    if(! visited[*vi] && out_degree(*vi,undigraph) <= 2)
      neighbors.push_back(*vi);

  }

}

bool extendScaffold (int & nextNode, std::vector<std::string> & scaffold, std::vector<char> &signs, UndirectedGraph &undigraph, std::vector<bool> &visited,
    std::map<int,std::string> &nodes_list_string,   std::map<std::string,int> &nodes_list_int){

  std::vector<int> neighbors;
  getNotVistedNeighbors(nextNode, undigraph, visited, neighbors);

  std::cout <<"Neighbors of the next node:" << neighbors.size() <<std::endl;

  if(neighbors.size() != 1){

    if(scaffold.size() == 0){//then it is an isolated contig

      visited[nextNode] = true;
      std::string node_sens = nodes_list_string.at(nextNode);
      std::string contig_name = node_sens.substr(0,node_sens.length()-4);
      std::string sens = node_sens.substr(node_sens.length()-3,3);
      std::cout << "Contig name: "<<contig_name<<std::endl;
      std::cout <<"Sense: "<<sens<<std::endl;
      if(sens.compare("beg")==0){
        std::string otherEx = contig_name.append(":end");

        std::cout <<"Other extr: "<< otherEx <<std::endl;
        int nodeOtherEx = nodes_list_int.at(otherEx);
        std::cout <<"Other extr in node: "<< nodeOtherEx <<std::endl;
        visited[nodeOtherEx] = true;

      }else {

        std::string otherEx = contig_name.append(":beg");
        std::cout <<"Other extr: "<< otherEx <<std::endl;
        int nodeOtherEx = nodes_list_int.at(otherEx);
        std::cout <<"Other extr in node: "<< nodeOtherEx <<std::endl;
        visited[nodeOtherEx] = true;


      }

      signs.push_back('+');
      scaffold.push_back(contig_name.substr(0,node_sens.length()-4));
    }

    return false;

  }else{

    int otherExtemity = neighbors[0];
    visited[nextNode] = true; //we start visiting but not yet sure if both extremities are ok

    getNotVistedNeighbors(otherExtemity , undigraph, visited, neighbors);


    if(neighbors.size() <= 1){ //we include the contig in the scaffold

      std::string node_sens = nodes_list_string.at(nextNode);
      std::cout << "1:"<< node_sens <<std::endl;
      std::string contig_name = node_sens.substr(0,node_sens.length()-4);
      std::cout << "2:"<< contig_name <<std::endl;
      std::string sens = node_sens.substr(node_sens.length()-3,3);
      std::cout << "3:"<< sens <<std::endl;

      visited[otherExtemity] = true;

      if(sens.compare("beg")==0){

        signs.push_back('+');

      }else {

        signs.push_back('-');
      }
      scaffold.push_back(contig_name);

      if(neighbors.size() == 0){ //the path has ended
        return false;
      }else {
        nextNode = neighbors[0];
        return true;
      }

    }else{//the contig is not visted in the end because the other extremity has many neighbors
      visited[nextNode] = false;
      return false;
    }
  }






}

*/

// Global values
int               Globals::pair_reads_length   =   300;
float             Globals::beginning_ratio     =   0.4;
unsigned long int Globals::window              = 10000;
unsigned long int Globals::max_contig_distance = 20000;
int               Globals::condition           =     1;
int               Globals::min_n_reads         =     3;
std::string       Globals::contig_file_name    =    "";
std::string       Globals::graph_file_name     =    "";
std::string       Globals::molecule_file_name  =    "";
std::string       Globals::scaffold_file_name  =    "";


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
    } else if ((arg == "-p") || (arg == "--pairReadsLength")){
      Globals::pair_reads_length = std::stoi(argv[++i]);
    } else if ((arg == "-m") || (arg == "--maxContigDist")){
      Globals::max_contig_distance = std::stoi(argv[++i]);
    } else if ((arg == "-c") || (arg == "--contigs")){
      Globals::contig_file_name = argv[++i];
    } else if ((arg == "-g") || (arg == "--graph")){
      Globals::graph_file_name = argv[++i];
    } else if ((arg == "-s") || (arg == "--scaffolds")){
      Globals::scaffold_file_name = argv[++i];
    } else {
      Globals::molecule_file_name = argv[i++];
    }
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
  merge_close_contigs(graph, contigs, scaffolds);
  print_scaffold(scaffolds, contigs);


/*
  std::cout << "number of contigs: "<< contigs_list.size() << std::endl;

  std::vector<std::pair<std::string, std::string>> arcs_list;
  std::vector<std::string> nodes_list;
  create_nodes(contigs_list,nodes_list);
  std::cout<<"nodes: "<<nodes_list.size()<<std::endl;

  std::vector<int> weight_list;
  create_arcs_with_size(contigs_list, arcs_list, weight_list);
  std::cout << "arcs: "<<arcs_list.size()<<std::endl;

  write_graph(contigs_list, arcs_list, weight_list);

  //traverse the graph


  std::ofstream scaffoldFile(Globals::scaffold_file.c_str(), std::ofstream::out);

  UndirectedGraph undigraph(arcs_list.size());


  std::map<std::string,int> nodes_list_int;
  std::map<int,std::string> nodes_list_string;
  std::vector<bool> visited;
  for(int i=0;i<nodes_list.size();i++){

    nodes_list_int.insert(std::pair<std::string,int>(nodes_list[i],i));
    nodes_list_string.insert(std::pair<int,std::string>(i, nodes_list[i]));
    visited.push_back(false);
  }



  for(int i=0; i<arcs_list.size(); i++){

    add_edge(nodes_list_int.at(arcs_list[i].first), nodes_list_int.at(arcs_list[i].second),undigraph);

  }

  std::cout<<"graph edges :" <<undigraph.m_edges.size() <<std::endl;
  std::cout<<"graph vertices :"<<undigraph.m_vertices.size() <<std::endl;

  std::vector<char> signs;
  std::vector<std::string> scaffold;

  std::vector <int> neighbors;

  int nextNode;
  bool unbranched;

  for(int i=0;i<contigs_list.size();i++){

    int ctg_beg_int = nodes_list_int.at(contigs_list[i].name+":beg");
    int ctg_end_int = nodes_list_int.at(contigs_list[i].name+":end");
    std::cout << contigs_list[i].name <<"\t" << ctg_beg_int << "\t"<< ctg_end_int <<std::endl;

    if (! visited[ctg_beg_int] ){ //if the contig was not yet included in a scaffold

      int degreeBeg = out_degree(ctg_beg_int, undigraph);
      int degreeEnd = out_degree(ctg_end_int, undigraph);

      std::cout <<"not visited"<<std::endl;

      if ((degreeBeg != 2) || (degreeEnd != 2)) {//if it's a unbranched path end

        std::cout <<"a branching path begins"<<degreeBeg <<"\t"<<degreeEnd<<std::endl;
        visited[ctg_beg_int] = true;
        visited[ctg_end_int] = true;

        if((degreeBeg == 1) ||
            (degreeBeg > 2) && (degreeEnd == 2) ||
            (degreeBeg > 2) && (degreeEnd > 2)) { //we start with the beginning

          signs.push_back('+');
        }else { //we start with the end

          signs.push_back('-');

        }

        scaffold.push_back(contigs_list[i].name);

        if ((degreeBeg == 1) && (degreeEnd == 1)){ //isolated contig
          std::cout <<"isolated contig" <<std::endl;
          print_scaffold(scaffold,signs,scaffoldFile);
          scaffold.clear();
          signs.clear();

        }else if ((degreeBeg >= 2) && (degreeEnd >=2)){ //the scaffold contains only one contig
          std::cout <<"contig with two extremites closed" <<std::endl;
          print_scaffold(scaffold,signs,scaffoldFile);
          scaffold.clear();
          signs.clear();

          getNotVistedNeighbors(ctg_beg_int , undigraph, visited, neighbors); //we start another scaffold from each neighbor of the begining node
          for(int i = 0; i < neighbors.size(); i++){

            nextNode = neighbors[i];
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          getNotVistedNeighbors(ctg_end_int , undigraph, visited, neighbors);
          for(int i = 0; i < neighbors.size(); i++){ //we start another scaffold from each neighbor of the ending node

            nextNode = neighbors[i];
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

        } else if(degreeBeg == 1) {
          std::cout <<"contig with begining extremity open" <<std::endl;
          if(degreeEnd > 2){ //the unbrancing path has ended with the ending node

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          getNotVistedNeighbors(ctg_end_int , undigraph, visited, neighbors);
          std::cout << "Neighbors size " <<neighbors.size() <<std::endl;
          for(int i = 0; i < neighbors.size(); i++){ //we start another scaffold from each neighbor of the ending node

            nextNode = neighbors[i];
            std::cout << "Neighbor: " << neighbors[i] <<std::endl;
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited, nodes_list_string, nodes_list_int);
            std::cout <<"continue? "<< unbranched <<std::endl;
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }
        } else if(degreeEnd == 1){

          std::cout << "contig with ending extremity open" <<std::endl;
          if(degreeBeg > 2){//the unbranching path has ended with the begining node

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          getNotVistedNeighbors(ctg_beg_int , undigraph, visited, neighbors);
          for(int i =0; i<neighbors.size();i++){  //we start another scaffold from each neighbor of the begining node

            nextNode = neighbors[i];
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

        }

      }

    }

  }

  scaffoldFile.close();
*/

  return EXIT_SUCCESS;
}
