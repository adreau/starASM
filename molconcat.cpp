#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <experimental/filesystem>
#include <vector>
#include <thread>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "Globals.h"
#include "Contig.h"

using namespace boost;

const char tab='\t';

typedef property < edge_weight_t, double >Weight;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, Weight > UndirectedGraph;


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
    << "\t-g, --graph           FILE  Output gfa file name \n"
    << std::endl;
}

void create_contigs(std::vector<Contig> &list_contigs){

  std::ifstream contig_bed(Globals::contigFile.c_str());
  std::string contig_line, ctg, ctg_name;
  int pos_beg, pos_end;

  while(getline(contig_bed,contig_line)){

    std::stringstream  splitstream(contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;

    ctg_name = ctg + "_";
    ctg_name = ctg_name + std::to_string(pos_beg);
    ctg_name = ctg_name + "_";
    ctg_name = ctg_name + std::to_string(pos_end);

    Contig contig(ctg_name, ctg, pos_beg, pos_end);

    list_contigs.push_back(contig);

  }
}

void add_molecules_to_contigs_extremites(std::vector<Contig> &nodes_list){

  std::ifstream molecule_list(Globals::molecule_file.c_str());
  std::string molecule_line, contig, barcode;
  int beg_pos, end_pos, noReads;

  int split_contig_count = 0;

  int split_ctg1_beg = nodes_list[split_contig_count].pos_beg;
  int split_ctg1_end = nodes_list[split_contig_count].pos_end;
  int length_ctg1 = split_ctg1_end - split_ctg1_beg;

  int split_ctg2_beg = nodes_list[split_contig_count+1].pos_beg;
  int split_ctg2_end = nodes_list[split_contig_count+1].pos_end;
  int length_ctg2 = split_ctg2_end - split_ctg2_beg;

  while(getline(molecule_list,molecule_line)){

    std::stringstream  splitstream(molecule_line);
    splitstream >> contig >> beg_pos >> end_pos >> barcode >> noReads;
    Molecule molecule (beg_pos,end_pos,barcode,noReads);


    if ( ( (beg_pos >= split_ctg1_end) && (contig.compare(nodes_list[split_contig_count].origin) == 0)) ||
        (contig.compare(nodes_list[split_contig_count].origin) > 0)){ //we reached the end of ctg1

      split_contig_count +=1;

      split_ctg1_beg = nodes_list[split_contig_count].pos_beg;
      split_ctg1_end = nodes_list[split_contig_count].pos_end;
      length_ctg1 = split_ctg1_end - split_ctg1_beg;

      split_ctg2_beg = nodes_list[split_contig_count+1].pos_beg;
      split_ctg2_end = nodes_list[split_contig_count+1].pos_end;
      length_ctg2 = split_ctg2_end - split_ctg2_beg;

    }

    if ( ( beg_pos >=  split_ctg1_beg ) && ( beg_pos <=  split_ctg1_beg + std::min(Globals::window,length_ctg1/2) ) &&
        ( end_pos <= split_ctg1_beg + length_ctg1 * (1 - Globals::beginning_ratio) ) ) { //condition begining

      nodes_list[split_contig_count].add_beg_molecule(molecule);

    }else if ( ( beg_pos >=  split_ctg1_beg + length_ctg1 * Globals::beginning_ratio)  && ( beg_pos <=  split_ctg1_end - Globals::pair_reads_length ) &&
        ( end_pos >=  split_ctg1_end - std::min(Globals::window,length_ctg1/2) )){//condition end


      if (contig.compare(nodes_list[split_contig_count+1].origin) != 0){ //ctg1 and ctg2 not the split of the same ctg

        nodes_list[split_contig_count].add_end_molecule(molecule);

      }else if(end_pos <= split_ctg2_beg +length_ctg2 * (1 - Globals::beginning_ratio) ){

        nodes_list[split_contig_count].add_end_molecule(molecule);

        if (end_pos >= split_ctg2_beg + Globals::pair_reads_length) { //condition begining of ctg2

          nodes_list[split_contig_count+1].add_beg_molecule(molecule);

        }
      }

    }

  }//end while read molecules

  for(int i=0; i<nodes_list.size();i++){
    sort(nodes_list[i].barcodes_beg.begin(), nodes_list[i].barcodes_beg.end());
    sort(nodes_list[i].barcodes_end.begin(), nodes_list[i].barcodes_end.end());
  }


}

void create_nodes(std::vector<Contig> &ctg_list, std::vector<std::string> &nodes){

  for(int i=0; i < ctg_list.size(); i++){
    nodes.push_back(ctg_list[i].name+":beg");
    nodes.push_back(ctg_list[i].name+":end");
  }

}


void create_arcs_with_size(std::vector<Contig> &ctg_list, std::vector<std::pair<std::string, std::string>> &arcs, std::vector<int> &weight){

  for(int i=0; i < ctg_list.size(); i++){ //construct the arcs between extremites of same contig
    arcs.push_back(std::pair<std::string, std::string>(ctg_list[i].name+":beg",ctg_list[i].name+":end"));
    weight.push_back(100000);
  }

  std::cout << "normal arcs: "<< arcs.size()<<std::endl;

  int diffCtgArc = 0;

  for(int i=0; i < ctg_list.size()-1; i++){

    for (int j = i; j < ctg_list.size(); j++){

      std::vector<int> connections;
      ctg_list[i].isNeighbourSize(ctg_list[j], connections);

      if(i==j) { //if the same extremity of a contig do not add arc
        connections[Link_types::BB] = 0;
        connections[Link_types::EE] = 0;
      }

      //add arc if there are three shared barcodes at least and the condition is satisfied
      if(connections[Link_types::BB] >= Globals::min_n_reads)  {
        arcs.push_back(std::pair<std::string, std::string>(ctg_list[i].name+":beg",ctg_list[j].name+":beg"));//bb
        weight.push_back(connections[Link_types::BB]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }
      if(connections[Link_types::BE] >= Globals::min_n_reads) {
        arcs.push_back(std::pair<std::string, std::string>(ctg_list[i].name+":beg",ctg_list[j].name+":end"));//be
        weight.push_back(connections[Link_types::BE]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }

      if(connections[Link_types::EB] >= Globals::min_n_reads) {
        arcs.push_back(std::pair<std::string, std::string>(ctg_list[i].name+":end",ctg_list[j].name+":beg"));//eb
        weight.push_back(connections[Link_types::EB]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }
      if(connections[Link_types::EE] >= Globals::min_n_reads) {
        arcs.push_back(std::pair<std::string, std::string>(ctg_list[i].name+":end",ctg_list[j].name+":end"));//ee
        weight.push_back(connections[Link_types::EE]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }


    }

    std::cout<< "progess:" << i<< " "<<arcs.size()<<std::endl;
  }

  std::cout << "Number of arcs between different original contigs: "<< diffCtgArc << std::endl;

}

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

void write_graph(std::vector<Contig> &contigs_list, std::vector<std::pair<std::string, std::string>> &arcs_list, std::vector<int> &weigth){

  std::ofstream graph_file(Globals::graph_file.c_str(), std::ofstream::out);

  graph_file << "H" << "\t" << "VN:Z:1.0"<< std::endl;

  int ctg_length;

  for(int i=0; i<contigs_list.size(); i++){

    ctg_length=contigs_list[i].pos_end-contigs_list[i].pos_beg;
    graph_file << "S" << "\t" << contigs_list[i].name <<"\t"<<"*"<< "\t" <<"LN:i:"<<ctg_length<<std::endl;

  }

  std::string ctg1,sign1,ctg2,sign2;

  for(int i=0; i<arcs_list.size();i++){

    get_canonical_ctg(arcs_list[i], ctg1,sign1, ctg2,sign2);

    if(ctg1.compare(ctg2)!=0)

      graph_file<< "L" <<"\t"<< ctg1 << "\t" << sign1 <<"\t"<<ctg2<<"\t"<<sign2<<"\t"<<"*"<<"\t"<<"BC:i:"<<weigth[i]<<std::endl;
  }


  graph_file.close();

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


// Global values
int         Globals::pair_reads_length =   300;
float       Globals::beginning_ratio   =   0.4;
int         Globals::window            = 10000;
int         Globals::condition         =     1;
int         Globals::min_n_reads       =     3;
std::string Globals::contigFile        =    "";
std::string Globals::graph_file        =    "";
std::string Globals::molecule_file     =    "";
std::string Globals::scaffold_file     =    "";


int main (int argc, char* argv[])
{
  if (argc < 2) {
    show_usage(argv[0]);
    return 1;
  }

  Globals globals;

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
      Globals::contigFile = argv[++i];
    } else if ((arg == "-g") || (arg == "--graph")){
      Globals::graph_file = argv[++i];
    } else if ((arg == "-s") || (arg == "--scaffolds")){
      Globals::scaffold_file = argv[++i];
    } else {
      Globals::molecule_file = argv[i++];
    }
  }

  std::vector<Contig> contigs_list;
  create_contigs(contigs_list);



  add_molecules_to_contigs_extremites(contigs_list);


  std::cout << "number of contigs: "<< contigs_list.size() << std::endl;

  std::cout <<"check molecules sets size" <<std::endl;

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

  return 0;
}
