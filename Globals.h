#ifndef DEF_GLOBALS
#define DEF_GLOBALS

static const int n_link_types = 4;

enum Link_types { BB = 0, BE = 1, EB = 2, EE = 3 };

// What happens if you swap first and second elements?
static const Link_types reverse_link_type[] = {Link_types::BB, Link_types::EB, Link_types::BE, Link_types::EE};

struct Globals {
  static int               pair_reads_length;
  static float             beginning_ratio;
  static unsigned long int window;
  static int               condition;
  static int               min_n_reads;
  static std::string       contig_file_name;
  static std::string       graph_file_name;
  static std::string       molecule_file_name;
  static std::string       scaffold_file_name;
};

#endif
