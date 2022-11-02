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

#include "Contig.h"
#include "Graph.h"

using namespace std;
using namespace boost;

namespace fs = std::experimental::filesystem;
const char tab='\t';

typedef property < edge_weight_t, double >Weight;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, Weight > UndirectedGraph;

const int pair_reads_length = 300;


static void show_usage(string name)
{
    cerr << "Usage: " << name << " <option(s)> MOLECULE FILE \n"
              << "Options:\n"
              << "\t-h, --help\t\tShow this help message\n"
              << "\t-w, --window INT\tWindow size for barcode consideration (default 10kb) \n"
              << "\t-a, --arcsCondition INT\tcondition used for connecting two contigs; values{1..8} \n"
              << "\t-c, --contigs FILE\tContig bed file name \n"
              << "\t-g, --graph FILE\tOutput gfa file name \n"
              << endl;


}

vector<Contig> create_contigs(string contigFile){

  ifstream contig_bed(contigFile.c_str());
  string contig_line, ctg, ctg_name;
  int pos_beg, pos_end;

  vector<Contig> list_contigs;

  while(getline(contig_bed,contig_line)){

    stringstream  splitstream(contig_line);
    splitstream >> ctg >> pos_beg >> pos_end;

    ctg_name = ctg + "_";
    ctg_name = ctg_name + to_string(pos_beg);
    ctg_name = ctg_name + "_";
    ctg_name = ctg_name + to_string(pos_end);

    Contig contig(ctg_name, ctg, pos_beg, pos_end);

    list_contigs.push_back(contig);

  }

  return list_contigs;
}

void add_molecules_to_contigs_extremites(vector<Contig> &nodes_list, string molecule_file, int window){

  ifstream molecule_list(molecule_file.c_str());
  string molecule_line, contig, barcode;
  int beg_pos, end_pos, noReads;

  int split_contig_count = 0;

  int split_ctg1_beg = nodes_list[split_contig_count].pos_beg;
  int split_ctg1_end = nodes_list[split_contig_count].pos_end;
  int length_ctg1 = split_ctg1_end - split_ctg1_beg;

  int split_ctg2_beg = nodes_list[split_contig_count+1].pos_beg;
  int split_ctg2_end = nodes_list[split_contig_count+1].pos_end;
  int length_ctg2 = split_ctg2_end - split_ctg2_beg;

  while(getline(molecule_list,molecule_line)){

    stringstream  splitstream(molecule_line);
    splitstream >> contig >> beg_pos >> end_pos >> barcode >> noReads;


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

    if ( ( beg_pos >=  split_ctg1_beg ) && ( beg_pos <=  split_ctg1_beg + min(window,length_ctg1/2) ) &&
            ( end_pos <= split_ctg1_beg + length_ctg1*0.6 ) ) { //condition begining

        nodes_list[split_contig_count].add_beg_molecule(Molecule(beg_pos,end_pos,barcode,noReads));

    }else if ( ( beg_pos >=  split_ctg1_beg + length_ctg1*0.4 )  && ( beg_pos <=  split_ctg1_end - pair_reads_length ) &&
            ( end_pos >=  split_ctg1_end - min(window,length_ctg1/2) )){//condition end


        if (contig.compare(nodes_list[split_contig_count+1].origin) != 0){ //ctg1 and ctg2 not the split of the same ctg

            nodes_list[split_contig_count].add_end_molecule(Molecule(beg_pos,end_pos,barcode,noReads));

        }else if(end_pos <= split_ctg2_beg +length_ctg2*0.6){

            nodes_list[split_contig_count].add_end_molecule(Molecule(beg_pos,end_pos,barcode,noReads));

            if (end_pos >= split_ctg2_beg + pair_reads_length) { //condition begining of ctg2

                nodes_list[split_contig_count+1].add_beg_molecule(Molecule(beg_pos,end_pos,barcode,noReads));

            }
        }

    }

  }//end while read molecules

  for(int i=0; i<nodes_list.size();i++){
    sort(nodes_list[i].barcodes_beg.begin(), nodes_list[i].barcodes_beg.end());
    sort(nodes_list[i].barcodes_end.begin(), nodes_list[i].barcodes_end.end());
  }


}

void create_nodes(vector<Contig> ctg_list, vector<string>& nodes){

  for(int i=0; i < ctg_list.size(); i++){
     nodes.push_back(ctg_list[i].name+":beg");
     nodes.push_back(ctg_list[i].name+":end");
  }

}


void create_arcs_with_size(vector<Contig> ctg_list, vector<pair<string, string>>& arcs, vector<int>& weight, int condition){

  for(int i=0; i < ctg_list.size(); i++){ //construct the arcs between extremites of same contig
     arcs.push_back(pair<string, string>(ctg_list[i].name+":beg",ctg_list[i].name+":end"));
     weight.push_back(100000);
  }

  cout << "normal arcs: "<< arcs.size()<<endl;

  int diffCtgArc = 0;

  for(int i=0; i < ctg_list.size()-1; i++){

    for (int j = i; j < ctg_list.size(); j++){

      vector<int> connections = ctg_list[i].isNeighbourSize(ctg_list[j], condition);

      if(i==j) { //if the same extremity of a contig do not add arc
        connections[0] = 0;
        connections[3] = 0;
      }

      //add arc if there are three shared barcodes at least and the condition is satisfied
      if(connections[0]>=3)  {
        arcs.push_back(pair<string, string>(ctg_list[i].name+":beg",ctg_list[j].name+":beg"));//bb
        weight.push_back(connections[0]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }
      if(connections[1]>=3) {
        arcs.push_back(pair<string, string>(ctg_list[i].name+":beg",ctg_list[j].name+":end"));//be
        weight.push_back(connections[1]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }

      if(connections[2]>=3) {
        arcs.push_back(pair<string, string>(ctg_list[i].name+":end",ctg_list[j].name+":beg"));//eb
        weight.push_back(connections[2]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }
      if(connections[3]>=3) {
        arcs.push_back(pair<string, string>(ctg_list[i].name+":end",ctg_list[j].name+":end"));//ee
        weight.push_back(connections[3]);
        if(ctg_list[i].origin.compare(ctg_list[j].origin)!=0) diffCtgArc++;
      }


    }

    cout<< "progess:" << i<< " "<<arcs.size()<<endl;
  }

  cout << "Number of arcs between different original contigs: "<< diffCtgArc << endl;

}

void get_canonical_ctg(pair<string,string> arc, string& ctg1, string& sign1, string& ctg2, string& sign2){

  string tmp_sign1, tmp_sign2;
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

void write_graph(vector<Contig> contigs_list, vector<pair<string, string>> arcs_list, vector<int> weigth, string graphFile){

  ofstream graph_file(graphFile.c_str(), ofstream::out);

  graph_file << "H" << "\t" << "VN:Z:1.0"<< endl;

  int ctg_length;

  for(int i=0; i<contigs_list.size(); i++){

      ctg_length=contigs_list[i].pos_end-contigs_list[i].pos_beg;
      graph_file << "S" << "\t" << contigs_list[i].name <<"\t"<<"*"<< "\t" <<"LN:i:"<<ctg_length<<endl;

  }

  string ctg1,sign1,ctg2,sign2;

  for(int i=0; i<arcs_list.size();i++){

    get_canonical_ctg(arcs_list[i], ctg1,sign1, ctg2,sign2);

    if(ctg1.compare(ctg2)!=0)

      graph_file<< "L" <<"\t"<< ctg1 << "\t" << sign1 <<"\t"<<ctg2<<"\t"<<sign2<<"\t"<<"*"<<"\t"<<"BC:i:"<<weigth[i]<<endl;
  }


  graph_file.close();

}


char opposing_sign(char sign){

  if(sign == '+') return '-';
  else return '+';

}

void print_scaffold( vector<string> scaffold, vector<char> signs, ofstream& scaffoldFile){

  cout << "Scaffold size: "<< scaffold.size()<<endl;
  if(signs[0] == '-') {
    for(int i=0; i<scaffold.size();i++)
      scaffoldFile << opposing_sign(signs[i]) <<scaffold[i]<<";";
  }else {
    for(int i=0; i<scaffold.size();i++)
      scaffoldFile << signs[i] <<scaffold[i]<<";";
  }
  scaffoldFile <<"\n";
}

vector<int> getNotVistedNeighbors(int ctg_beg_int , UndirectedGraph undigraph, vector<bool> visited){

  vector<int> neighbors;

  typename graph_traits<UndirectedGraph>::adjacency_iterator vi, vi_end;

  for (boost::tie(vi, vi_end) = adjacent_vertices(ctg_beg_int, undigraph); vi != vi_end; ++vi){

    cout << "Neighbors all: "<<*vi<<endl;
    if(! visited[*vi] && out_degree(*vi,undigraph) <= 2)
      neighbors.push_back(*vi);

  }

  return neighbors;


}

bool extendScaffold (int & nextNode, vector<string> & scaffold, vector<char> & signs, UndirectedGraph undigraph, vector<bool> & visited,
                     map<int,string> nodes_list_string,   map<string,int> nodes_list_int){

  vector<int> neighbors = getNotVistedNeighbors(nextNode , undigraph, visited);

  cout <<"Neighbors of the next node:" << neighbors.size() <<endl;

  if(neighbors.size() != 1){

    if(scaffold.size() == 0){//then it is an isolated contig

      visited[nextNode] = true;
      string node_sens = nodes_list_string.at(nextNode);
      string contig_name = node_sens.substr(0,node_sens.length()-4);
      string sens = node_sens.substr(node_sens.length()-3,3);
      cout << "Contig name: "<<contig_name<<endl;
      cout <<"Sense: "<<sens<<endl;
      if(sens.compare("beg")==0){
        string otherEx = contig_name.append(":end");

        cout <<"Other extr: "<< otherEx <<endl;
        int nodeOtherEx = nodes_list_int.at(otherEx);
        cout <<"Other extr in node: "<< nodeOtherEx <<endl;
        visited[nodeOtherEx] = true;

      }else {

         string otherEx = contig_name.append(":beg");
        cout <<"Other extr: "<< otherEx <<endl;
        int nodeOtherEx = nodes_list_int.at(otherEx);
        cout <<"Other extr in node: "<< nodeOtherEx <<endl;
        visited[nodeOtherEx] = true;


      }

      signs.push_back('+');
      scaffold.push_back(contig_name.substr(0,node_sens.length()-4));
    }

    return false;

  }else{

    int otherExtemity = neighbors[0];
    visited[nextNode] = true; //we start visiting but not yet sure if both extremities are ok

    neighbors = getNotVistedNeighbors(otherExtemity , undigraph, visited);


    if(neighbors.size() <= 1){ //we include the contig in the scaffold

      string node_sens = nodes_list_string.at(nextNode);
      cout << "1:"<< node_sens <<endl;
      string contig_name = node_sens.substr(0,node_sens.length()-4);
      cout << "2:"<< contig_name <<endl;
      string sens = node_sens.substr(node_sens.length()-3,3);
      cout << "3:"<< sens <<endl;

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

int main (int argc, char* argv[])
{
  if (argc < 2) {
      show_usage(argv[0]);
      return 1;
  }

  int window = 10000;
  int condition = 1;
  string contigFile = "";
  string graph_file = "";
  string molecule_file = "";
  string scaffold_file = "";
  for (int i = 1; i < argc; ++i) {
      string arg = argv[i];
      if ((arg == "-h") || (arg == "--help")) {
        show_usage(argv[0]);
        return 0;
      } else if ((arg == "-w") || (arg == "--window")) { //maximal distance between two read pairs of the same molecule
        window = stoi(argv[++i]);
      } else if ((arg == "-a") || (arg == "--arcsCondition")){
        condition = stoi(argv[++i]);
      } else if ((arg == "-c") || (arg == "--contigs")){
        contigFile = argv[++i];
      } else if ((arg == "-g") || (arg == "--graph")){
        graph_file = argv[++i];
      } else if ((arg == "-s") || (arg == "--scaffolds")){
        scaffold_file = argv[++i];
      } else {
        molecule_file = argv[i++];
      }
  }

  vector<Contig> contigs_list = create_contigs(contigFile);



  add_molecules_to_contigs_extremites(contigs_list, molecule_file, window);


  cout << "number of contigs: "<< contigs_list.size() << endl;

  cout <<"check molecules sets size" <<endl;

  for (int i=0; i<contigs_list.size() ; i++){

    cout << contigs_list[i].name << " " <<contigs_list[i].molecules_beg.size()<< " "<<contigs_list[i].molecules_end.size()<<endl;
  }

  vector<pair<string, string>> arcs_list;
  vector<string> nodes_list;
  create_nodes(contigs_list,nodes_list);
  cout<<"nodes: "<<nodes_list.size()<<endl;

  vector<int> weight_list;
  create_arcs_with_size(contigs_list, arcs_list, weight_list, condition);
  cout << "arcs: "<<arcs_list.size()<<endl;

  write_graph(contigs_list, arcs_list, weight_list, graph_file);

  //traverse the graph


  ofstream scaffoldFile(scaffold_file.c_str(), ofstream::out);

  UndirectedGraph undigraph(arcs_list.size());


  map<string,int> nodes_list_int;
  map<int,string> nodes_list_string;
  vector<bool> visited;
  for(int i=0;i<nodes_list.size();i++){

    nodes_list_int.insert(pair<string,int>(nodes_list[i],i));
    nodes_list_string.insert(pair<int,string>(i, nodes_list[i]));
    visited.push_back(false);
  }



  for(int i=0; i<arcs_list.size(); i++){

    add_edge(nodes_list_int.at(arcs_list[i].first), nodes_list_int.at(arcs_list[i].second),undigraph);

  }

  cout<<"graph edges :" <<undigraph.m_edges.size() <<endl;
  cout<<"graph vertices :"<<undigraph.m_vertices.size() <<endl;

  vector<char> signs;
  vector<string> scaffold;

  vector <int> neighbors;

  int nextNode;
  bool unbranched;

  for(int i=0;i<contigs_list.size();i++){

    int ctg_beg_int = nodes_list_int.at(contigs_list[i].name+":beg");
    int ctg_end_int = nodes_list_int.at(contigs_list[i].name+":end");
    cout << contigs_list[i].name <<"\t" << ctg_beg_int << "\t"<< ctg_end_int <<endl;

    if (! visited[ctg_beg_int] ){ //if the contig was not yet included in a scaffold

      int degreeBeg = out_degree(ctg_beg_int, undigraph);
      int degreeEnd = out_degree(ctg_end_int, undigraph);

      cout <<"not visited"<<endl;

      if ((degreeBeg != 2) || (degreeEnd != 2)) {//if it's a unbrached path end

        cout <<"a branching path begins"<<degreeBeg <<"\t"<<degreeEnd<<endl;
        visited[ctg_beg_int] = true;
        visited[ctg_end_int] = true;

        if((degreeBeg == 1) ||
          (degreeBeg>2) && (degreeEnd == 2) ||
          (degreeBeg>2) && (degreeEnd >2)) { //we start with the begining

          signs.push_back('+');
        }else { //we start with the end

          signs.push_back('-');

        }

        scaffold.push_back(contigs_list[i].name);

        if ((degreeBeg == 1) && (degreeEnd == 1)){ //isolated contig
          cout <<"isolated contig" <<endl;
          print_scaffold(scaffold,signs,scaffoldFile);
          scaffold.clear();
          signs.clear();

        }else if ((degreeBeg >= 2) && (degreeEnd >=2)){ //the scaffold contains only one contig
          cout <<"contig with two extremites closed" <<endl;
          print_scaffold(scaffold,signs,scaffoldFile);
          scaffold.clear();
          signs.clear();

          neighbors = getNotVistedNeighbors(ctg_beg_int , undigraph, visited); //we start another scaffold from each neighbor of the begining node
          for(int i =0; i<neighbors.size();i++){

            nextNode = neighbors[i];
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          neighbors = getNotVistedNeighbors(ctg_end_int , undigraph, visited);
          for(int i =0; i<neighbors.size();i++){ //we start another scaffold from each neighbor of the ending node

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
          cout <<"contig with begining extremity open" <<endl;
          if(degreeEnd >2){ //the unbrancing path has ended with the ending node

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          neighbors = getNotVistedNeighbors(ctg_end_int , undigraph, visited);
          cout << "Neighbors size " <<neighbors.size() <<endl;
          for(int i =0; i<neighbors.size();i++){ //we start another scaffold from each neighbor of the ending node

            nextNode = neighbors[i];
            cout << "Neighbor: " <<neighbors[i] <<endl;
            unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited, nodes_list_string, nodes_list_int);
            cout <<"continue? "<<unbranched <<endl;
            while(unbranched){
              unbranched = extendScaffold (nextNode,scaffold,signs, undigraph,visited,nodes_list_string, nodes_list_int);
            }

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }
        } else if(degreeEnd == 1){

          cout <<"contig with ending extremity open" <<endl;
          if(degreeBeg >2){//the unbrancing path has ended with the begining node

            print_scaffold(scaffold,signs,scaffoldFile);
            scaffold.clear();
            signs.clear();

          }

          neighbors = getNotVistedNeighbors(ctg_beg_int , undigraph, visited);
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
