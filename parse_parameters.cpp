#include <iostream>

#include "parse_parameters.h"


void show_usage(char *name) {
  std::cerr << "Usage: " << name << " <option(s)>\n"
    << "Options:\n"
    << "\t-h, --help                    Show this help message (parameters with stars '*' are compulsory)\n"
    << "=== Input ===\n"
    << "\t-f, --contigs         FILE  * Input contigs in FASTA format\n"
    << "\t-i, --molecule        FILE    Input molecule file\n"
    << "\t-c, --joins           FILE    Contig bed file name (result of splitASM) \n"
    << "\t-w, --window          INT     Window size for barcode consideration (default: "                             << Globals::window              << ") \n"
    << "=== Output ===\n"
    << "\t-o, --output          FILE  * Output FASTA file\n"
    << "\t-s, --scaffolds       FILE    Output scaffolds file name\n"
    << "\t-g, --graph           FILE    Output GFA file name\n"
    << "\t-p, --mapping         FILE    Log where the molecule map with respect to the contigs\n"
    << "=== JoinASM parameters ===\n"
    << "\t-a, --arcsCondition   INT     Condition used for connecting two contigs; values{1..8} (default: "           << Globals::condition           << ", lower is more strict) \n"
    << "\t-r, --nReads          INT     Min number of common barcodes to get a links (default: "                      << Globals::min_n_reads         << ")\n"
    << "\t-b, --begRatio        FLOAT   Ratio of the contig size that is considered as the beginning part (default: " << Globals::beginning_ratio     << ", should be less than 0.5)\n"
    << "\t-v, --minOverlap      INT     Minimum overlap between a molecule and a contig (default: "                   << Globals::min_overlap         << ")\n"
    << "\t-m, --maxContigDist   INT     Merge contigs if they are separated by not more that N bp (default: "         << Globals::max_contig_distance << ")\n"
    << "\t-l, --fillerSize      INT     Size of the stretch of Ns between the sequences (default: "                   << Globals::filler_size         << ")\n";
}


void parse_parameters (int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      exit(EXIT_SUCCESS);
    } else if ((arg == "-w") || (arg == "--window")) { //maximal distance between two read pairs of the same molecule
      Globals::window = std::stoi(argv[++i]);
    } else if ((arg == "-a") || (arg == "--arcsCondition")){
      Globals::condition = std::stoi(argv[++i]);
    } else if ((arg == "-r") || (arg == "--minReads")){
      Globals::min_n_reads = std::stoi(argv[++i]);
    } else if ((arg == "-b") || (arg == "--begRatio")){
      Globals::beginning_ratio = std::stof(argv[++i]);
    } else if ((arg == "-v") || (arg == "--minOverlap")){
      Globals::min_overlap = std::stoi(argv[++i]);
    } else if ((arg == "-l") || (arg == "--fillerSize")){
      Globals::filler_size = std::atoi(argv[++i]);
    } else if ((arg == "-m") || (arg == "--maxContigDist")){
      Globals::max_contig_distance = std::stoi(argv[++i]);
    } else if ((arg == "-c") || (arg == "--joins")){
      Globals::joins_file_name = argv[++i];
    } else if ((arg == "-g") || (arg == "--graph")){
      Globals::graph_file_name = argv[++i];
    } else if ((arg == "-s") || (arg == "--scaffolds")){
      Globals::scaffold_file_name = argv[++i];
    } else if ((arg == "-p") || (arg == "--mapping")){
      Globals::mapping_file_name = argv[++i];
    } else if ((arg == "-o") || (arg == "--output")){
      Globals::fasta_file_name = argv[++i];
    } else if ((arg == "-f") || (arg == "--contigs")){
      Globals::contigs_file_name = argv[++i];
    } else if ((arg == "-i") || (arg == "--molecules")){
      Globals::molecule_file_name = argv[++i];
    } else {
      std::cerr << "Error!  Parameter '" << argv[i] << "' is not understood.\nExiting.\n";
      exit(EXIT_FAILURE);
    }
  }
}
