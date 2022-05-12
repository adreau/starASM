#ifndef DEF_Graph
#define DEF_Graph

#include <string>
#include <vector>

#include "Contig.h"

using namespace std;

class Graph {
  
 public:

  Graph(vector<Contig> &v, vector<std::pair<Contig, Contig>> &e)
      : v_(v), e_(e) {}

  vector<Contig> v_;

  vector<std::pair<Contig, Contig>> e_;


};


#endif
