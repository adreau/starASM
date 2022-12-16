#ifndef DEF_GLOBALS
#define DEF_GLOBALS

static const int n_link_types = 4;

static const size_t unset_value = -1;

enum Link_types { BB = 0, BE = 1, EB = 2, EE = 3 };

// What happens if you swap first and second elements?
static const Link_types reverse_link_type[] = {Link_types::BB, Link_types::EB, Link_types::BE, Link_types::EE};

struct Globals {
  static unsigned long int min_overlap;
  static float             beginning_ratio;
  static unsigned long int window;
  static unsigned long int max_contig_distance;
  static int               condition;
  static int               min_n_reads;
  static std::string       contig_file_name;
  static std::string       graph_file_name;
  static std::string       molecule_file_name;
  static std::string       scaffold_file_name;
  static std::string       mapping_file_name;
};

#endif
