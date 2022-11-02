#ifndef DEF_Graph
#define DEF_Graph

#include <string>
#include <vector>

#include "Contig.h"

class Graph {
  
 public:

  Graph(std::vector<Contig> &v, std::vector<std::pair<Contig, Contig>> &e)
      : v_(v), e_(e) {}

  std::vector<Contig> v_;

  std::vector<std::pair<Contig, Contig>> e_;


};


#endif
