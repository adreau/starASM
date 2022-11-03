#ifndef DEF_GLOBALS
#define DEF_GLOBALS

static const int n_link_types = 4;

enum Link_types { BB = 0, BE = 1, EB = 2, EE = 3 };

struct Globals {
  static int         pair_reads_length;
  static float       beginning_ratio;
  static int         window;
  static int         condition;
  static int         min_n_reads;
  static std::string contigFile;
  static std::string graph_file;
  static std::string molecule_file;
  static std::string scaffold_file;
};

#endif
