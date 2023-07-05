#include <algorithm>
#include "globals.h"

std::vector < std::string >                     Globals::chrs;
std::unordered_map < std::string, std::size_t > Globals::chrids;
std::vector < size_t >                          Globals::chr_sizes;
Sequences                                       Globals::sequences;


std::string   Globals::input_file_name;
unsigned int  Globals::max_read_distance      {              60000 };
unsigned int  Globals::min_mapq               {                  0 };
unsigned int  Globals::min_mapq_solid         {                 30 };
unsigned int  Globals::min_len_solid          {                120 };
unsigned int  Globals::min_n_reads_barcode    {                  2 };
std::string   Globals::output_molecules_file_name;
std::string   Globals::output_split_file_name;
unsigned long Globals::seed                   {                  0 };
float         Globals::threshold              {                0.0 };
unsigned long Globals::window                 {              10000 };
std::string   Globals::counts_file_name;
std::string   Globals::scores_file_name;
std::string   Globals::contigs_file_name;
std::string   Globals::input_molecules_file_name;
unsigned long Globals::n_sample               {              10000 };
int           Globals::min_ctg_size           {              20000 };


unsigned long Globals::min_overlap            {                300 };
float         Globals::beginning_ratio        {                0.4 };
unsigned long Globals::max_contig_distance    {              20000 };
float         Globals::jaccard                {                0.5 };
unsigned int  Globals::min_n_reads            {                  3 };
std::string   Globals::input_split_file_name;
std::string   Globals::graph_file_name;
std::string   Globals::molecule_file_name;
std::string   Globals::scaffold_file_name;
std::string   Globals::mapping_file_name;
std::string   Globals::cis_link_file_name;
std::string   Globals::trans_link_file_name;
std::string   Globals::output_file_name;
unsigned int  Globals::filler_size            {                100 };


bool starts_with(const std::string &s1, const std::string &s2) {
  return s1.rfind(s2, 0) == 0;
}

unsigned int count_different (std::vector < unsigned long > &names) {
  std::sort(names.begin(), names.end());
  auto last = std::unique(names.begin(), names.end());
  return std::distance(names.begin(), last);
}
