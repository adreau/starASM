#include "globals.h"


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
