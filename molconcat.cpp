#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cassert>

// #include <boost/config.hpp>
// #include <boost/graph/adjacency_list.hpp>

#include "Globals.h"
#include "Contig.h"
#include "Graph.h"

// using namespace boost;

const char tab='\t';

// typedef property < edge_weight_t, double >Weight;
// typedef adjacency_list < vecS, vecS, undirectedS, no_property, Weight > UndirectedGraph;


static void show_usage(char *name)
{
  std::cerr << "Usage: " << name << " <option(s)> MOLECULE FILE \n"
    << "Options:\n"
    << "\t-h, --help                  Show this help message\n"
    << "\t-w, --window          INT   Window size for barcode consideration (default: " << Globals::window << ") \n"
    << "\t-a, --arcsCondition   INT   Condition used for connecting two contigs; values{1..8} (default: " << Globals::condition << ", lower is more strict) \n"
    << "\t-r, --nReads          INT   Min number of common barcodes to get a links (default: " << Globals::min_n_reads << ")\n"
    << "\t-b, --begRatio        FLOAT Ratio of the contig size that is considered as the beginning part (default: " << Globals::beginning_ratio << ", should be less than 0.5)\n"
    << "\t-p, --pairReadsLength INT   ??? (default: " << Globals::pair_reads_length << ")\n"
    << "\t-c, --contigs         FILE  Contig bed file name \n"
    << "\t-s, --scaffolds       FILE  Scaffolds file name \n"
    << "\t-g, --graph           FILE  Output gfa file name \n"
    << std::endl;
}

// Read/store a set of (split) contigs in a bed file.
void create_contigs(std::vector < Contig > &contig_list, std::unordered_map < std::string, size_t > &contig_ids){

  std::ifstream contig_bed(Globals::contig_file_name.c_str());
  std::string contig_line, ctg;
  int pos_beg, pos_end;

  while (getline(contig_bed, contig_line)){

    std::stringstream splitstream (contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;

    // If the contig is split, the second (third, etc.) part should be appened to the contig
    if ((! contig_list.empty()) && (ctg == contig_list.back().name)) {
      contig_list.back().addPart(pos_beg, pos_end);
    }
    else {
      auto pos = contig_ids.find(ctg);
      // We have a new contig
      if (pos == contig_ids.end()) {
        contig_ids[ctg] = contig_list.size();
        contig_list.emplace_back(ctg, pos_beg, pos_end);
      }
      // If the file is not ordered (but why?) the contig may be already stored
      else {
        contig_list[pos->second].addPart(pos_beg, pos_end);
      }
    }
  }
}

int intersectMoleculesSize(std::vector < std::string > &b1, std::vector < std::string > &b2){

    std::vector < std::string > common_barcodes;

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

    return 0;
}


void add_molecules_to_contigs_extremites(std::vector < Contig > &contigs, std::unordered_map < std::string, size_t > &contig_ids){

  std::ifstream molecule_list(Globals::molecule_file_name.c_str());
  std::string molecule_line, barcode, ctg, prevCtg;
  unsigned long int beg_pos, end_pos, nReads;
  size_t ctg_id, prevCtg_id = 0;

  long unsigned int n_lines;
  for (n_lines = 0; getline(molecule_list, molecule_line); ++n_lines){

    std::stringstream  splitstream(molecule_line);
    splitstream >> ctg >> beg_pos >> end_pos >> barcode >> nReads;

    if (ctg == prevCtg) {
      ctg_id = prevCtg_id;
    }
    else {
      auto pos = contig_ids.find(ctg);
      if (pos == contig_ids.end()) {
        std::cerr << "Error, contig '" << ctg << "' (found in molecule file) is not present in contig file.\n";
        exit(EXIT_FAILURE);
      }
      ctg_id     = pos->second;
      prevCtg    = ctg;
      prevCtg_id = ctg_id;
    }

    for (ContigPart &contigPart: contigs[ctg_id].contigParts) {
      // This is the condition to set a barcode to the beginning of a contig part
      if ((beg_pos >= contigPart.begin) && (beg_pos <= contigPart.begin + std::min(Globals::window, contigPart.getSize() / 2)) &&
          (end_pos <= contigPart.begin + contigPart.getSize() * (1 - Globals::beginning_ratio))) {
        contigPart.add_beg_molecule(barcode);
      }
      // This is the condition to set a barcode to the end of a contig part
      else if ((beg_pos >= contigPart.begin + contigPart.getSize() * Globals::beginning_ratio) && (beg_pos <= contigPart.end - Globals::pair_reads_length) &&
          (end_pos >= contigPart.end - std::min(Globals::window, contigPart.getSize() / 2))) {
        contigPart.add_end_molecule(barcode);
      }
    }
    if (n_lines % 10000000 == 0) std::cout << n_lines << " lines read.\r" << std::flush;

  }//end while read molecules
  std::cout << n_lines << " lines read.\n";

  for (Contig &contig: contigs) {
    contig.sort_barcodes();
  }
}

void create_nodes (std::vector < Contig > &contigs, Graph &graph) {

  for (size_t contigId = 0; contigId < contigs.size(); ++contigId) {
    Contig &contig = contigs[contigId];
    for (size_t contigPartId = 0; contigPartId < contig.contigParts.size(); ++contigPartId) {
      graph.add_node(contigId, contigPartId);
    }
  }
}


void create_arcs (std::vector<Contig> &contigs, Graph &graph) {

  size_t nodeId1 = 0;
  size_t nodeId2 = 0;
  for (size_t contigId1 = 0; contigId1 < contigs.size(); ++contigId1) {
    Contig &contig1 = contigs[contigId1];
    for (size_t contigPartId1 = 0; contigPartId1 < contig1.contigParts.size(); ++contigPartId1) {
      ContigPart &contigPart1 = contig1.contigParts[contigPartId1];
      assert(graph.nodes[nodeId1].contigId     == contigId1);
      assert(graph.nodes[nodeId1].contigPartId == contigPartId1);
      size_t contigPartId2 = contigPartId1;
      for (size_t contigId2 = contigId1; contigId2 < contigs.size(); ++contigId2) {
        Contig &contig2 = contigs[contigId2];
        for (; contigPartId2 < contig2.contigParts.size(); ++contigPartId2) {
          ContigPart &contigPart2 = contig2.contigParts[contigPartId2];
          assert(graph.nodes[nodeId2].contigId     == contigId2);
          assert(graph.nodes[nodeId2].contigPartId == contigPartId2);
          // Use a strange trick to keep nodeId1 and nodeId2 synchronized
          nodeId2 = nodeId1;
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
      ++nodeId1;
    }
  }
}


void get_node_name (std::vector < Contig > &contigs, Graph &graph, size_t nodeId, std::string &name) {
  Node       &node       = graph.nodes[nodeId];
  Contig     &contig     = contigs[node.contigId];
  ContigPart &contigPart = contig.contigParts[node.contigPartId];
  name = contig.name + "_" + std::to_string(contigPart.begin) + "_" + std::to_string(contigPart.end);
}


void write_graph (std::vector < Contig > &contigs, Graph &graph) {

  std::ofstream graph_file (Globals::graph_file_name.c_str(), std::ofstream::out);
  std::string ctg1, ctg2;

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
            graph_file << "L" << "\t" << ctg2 << "\t-\t" << ctg1 << "\t+\n";
            break;
          case Link_types::BE:
            graph_file << "L" << "\t" << ctg2 << "\t+\t" << ctg1 << "\t+\n";
            break;
          case Link_types::EB:
            graph_file << "L" << "\t" << ctg1 << "\t+\t" << ctg2 << "\t+\n";
            break;
          case Link_types::EE:
            graph_file << "L" << "\t" << ctg1 << "\t+\t" << ctg2 << "\t-\n";
            break;
        }
      }
    }
  }
  graph_file.close();
}


// Remove arcs when bifurcations are observed
// Arcs are actually not removed, but replaced by a given value
void remove_bifurcations (Graph &graph) {

  size_t not_set = graph.nodes.size();
  for (size_t nodeId1 = 0; nodeId1 < graph.nodes.size(); ++nodeId1) {
    Node &node1 = graph.nodes[nodeId1];
    int n_begin = 0;
    int n_end   = 0;
    for (Edge &edge: node1.edges) {
      if (edge.nodeId != not_set) {
        if ((edge.link_type == Link_types::BB) || (edge.link_type == Link_types::BE)) {
          ++n_begin;
        }
        else {
          ++n_end;
        }
      }
      // Found a bifurcation at the beginning of the contig
      if (n_begin > 1) {
        for (Edge &edge1: node1.edges) {
          if ((edge1.link_type == Link_types::BB) || (edge1.link_type == Link_types::BE)) {
            // Discard this edge
            edge1.nodeId = not_set;
            // Find the corresponding edge in the other node (each edge is present twice)
            size_t nodeId2 = edge1.nodeId;
            for (Edge &edge2: graph.nodes[nodeId2].edges) {
              if ((edge2.nodeId == nodeId1) && ((edge2.link_type == Link_types::BB) || (edge2.link_type == Link_types::EB))) {
                edge2.nodeId = not_set;
              }
            }
          }
        }
      }
      // Found a bifurcation at the end of the contig
      if (n_end > 1) {
        for (Edge &edge1: node1.edges) {
          if ((edge1.link_type == Link_types::EB) || (edge1.link_type == Link_types::EE)) {
            edge1.nodeId = not_set;
            size_t nodeId2 = edge1.nodeId;
            for (Edge &edge2: graph.nodes[nodeId2].edges) {
              if ((edge2.nodeId == nodeId1) && ((edge2.link_type == Link_types::BE) || (edge2.link_type == Link_types::EE))) {
                edge2.nodeId = not_set;
              }
            }
          }
        }
      }
    }
  }
}


size_t find_first_scaffold_end (Graph &graph, size_t nodeIdStart, std::vector < bool > seen_nodes) {

  size_t not_set = graph.nodes.size();
  for (size_t nodeId = nodeIdStart; nodeId < graph.nodes.size(); ++nodeId) {
    if (! seen_nodes[nodeId]) {
      Node        &node    = graph.nodes[nodeId];
      unsigned int n_nodes = 0;
      for (Edge &edge: node.edges) {
        if (edge.nodeId != not_set) {
          ++n_nodes;
        }
        if (n_nodes < 2) {
          return nodeId;
        }
      }
    }
  }
  return not_set;
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
int               Globals::pair_reads_length  =   300;
float             Globals::beginning_ratio    =   0.4;
unsigned long int Globals::window             = 10000;
int               Globals::condition          =     1;
int               Globals::min_n_reads        =     3;
std::string       Globals::contig_file_name   =    "";
std::string       Globals::graph_file_name    =    "";
std::string       Globals::molecule_file_name =    "";
std::string       Globals::scaffold_file_name =    "";


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

  std::vector < Contig > contigs_list;
  std::unordered_map < std::string, size_t > contig_ids;
  create_contigs(contigs_list, contig_ids);

  add_molecules_to_contigs_extremites(contigs_list, contig_ids);
  Graph graph;
  create_nodes(contigs_list, graph);
  create_arcs(contigs_list, graph);
  write_graph(contigs_list, graph);


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
