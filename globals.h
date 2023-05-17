#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstddef>
#include <unordered_map>
#include <vector>
#include <string>


static const int n_link_types = 4;

static const size_t unset_value = -1;

enum Link_types { BB = 0, BE = 1, EB = 2, EE = 3 };

// What happens if you swap first and second elements?
static const Link_types reverse_link_type[] = {Link_types::BB, Link_types::EB, Link_types::BE, Link_types::EE};

struct Globals {
  static std::vector < std::string >                     chrs;
  static std::unordered_map < std::string, std::size_t > chrids;
  static std::vector < std::size_t >                     chr_sizes;

  static unsigned int  max_read_distance;
  static unsigned int  min_mapq;
  static unsigned int  min_mapq_solid;
  static unsigned int  min_len_solid;
  static unsigned int  min_n_reads_barcode;
  static std::string   output_molecules_file_name;
  static std::string   output_split_file_name;
  static float         threshold;
  static unsigned long window;
  static std::string   counts_file_name;
  static std::string   scores_file_name;
  static std::string   contigs_file_name;
  static std::string   input_molecules_file_name;
  static unsigned long n_sample;
  static int           min_ctg_size;

  static unsigned long min_overlap;
  static float         beginning_ratio;
  static unsigned long max_contig_distance;
  static int           condition;
  static unsigned int  min_n_reads;
  static std::string   input_split_file_name;
  static std::string   graph_file_name;
  static std::string   molecule_file_name;
  static std::string   scaffold_file_name;
  static std::string   mapping_file_name;
  static std::string   output_file_name;
  static unsigned int  filler_size;

};


struct Molecule_stat {
  double coverage;
  double length;
  double read_density;
  double start;
  double end;
};


using Molecule_stats = std::vector < std::vector < Molecule_stat > >;


bool starts_with (const std::string &s1, const std::string &s2);


#endif
