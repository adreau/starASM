#include <iostream>

#include "parse_parameters.h"


void show_usage(char *name) {
  std::cerr << "Usage: " << name << " <option(s)>\n"
    << "Options:\n"
    << "    -h, --help                   Show this help message\n"
    << "=== Molecule creation step ===\n"
    << "  * Inputs\n"
    << "    -i, --input           FILE   Reads file in BAM format\n"
    << "    -f, --contigs         FILE   Contig file name\n"
    << "  * Output\n"
    << "    -m, --outputMolecule  FILE   Output file name\n"
    << "  * Other parameters\n"
    << "    -d, --distance        INT    Max. distance to cluster reads (default: "                                    << Globals::max_read_distance          << ")\n"
    << "    -q, --mapq            INT    Min. MAPQ (default: "                                                         << Globals::min_mapq                   << ")\n"
    << "    -Q, --mapqSolid       INT    Min. MAPQ for solid reads (default: "                                         << Globals::min_mapq_solid             << ")\n"
    << "    -l, --lengthSolid     INT    Min. mapping size for solid reads (default: "                                 << Globals::min_len_solid              << ")\n"
    << "    -n, --barcode         INT    Min. #reads per barcode (default: "                                           << Globals::min_n_reads_barcode        << ")\n"
    << "======= Split ASM step =======\n"
    << "  * Inputs\n"
    << "    -f, --contigs         FILE   Input contigs in FASTA format\n"
    << "    -M, --inputMolecule   FILE   Input molecule file (output of previous step)\n"
    << "    -w, --window          INT    Window size for barcode consideration (default: "                             << Globals::window                     << ") \n"
    << "  * Output\n"
    << "    -s, --outputSplit     FILE   Output split file in BED format\n"
    << "  * Other parameters\n"
    << "    -t, --threshold       FLOAT  Stringency threshold, higher is less stringent (default "                     << Globals::threshold                  << ")\n"
    << "    -e, --sampleSize      INT    Sample size for outlier detection (0 = input size (default: "                 << Globals::n_sample                   << ")\n"
    << "    -z, --minSize         INT    Minimum contig size (default: "                                               << Globals::min_ctg_size               << ")\n"
    << "    -E, --seed            INT    Seed (0 = time, default: "                                                    << Globals::seed                       << ")\n"
    << "  * Log files\n"
    << "    -a, --counts          FILE   Write raw counts to file\n"
    << "    -A, --scores          FILE   Write scores to file\n"
    << "======= Join ASM step ========\n"
    << "  * Inputs\n"
    << "    -f, --contigs         FILE   Input contigs in FASTA format\n"
    << "    -M, --inputMolecule   FILE   Input molecule file (output of previous step)\n"
    << "    -S, --inputSplit      FILE   Input split file in BED format (output of previous step)\n"
    << "    -w, --window          INT    Window size for barcode consideration (default: "                             << Globals::window                     << ")\n"
    << "  * Output\n"
    << "    -o, --output          FILE   Output scaffold file in FASTA format\n"
    << "  * Other parameters\n"
    << "    -c, --arcsCondition   INT    Condition used for connecting two contigs; values{1..8} (default: "           << Globals::condition                  << ", lower is more strict) \n"
    << "    -r, --nReads          INT    Min number of common barcodes to get a links (default: "                      << Globals::min_n_reads                << ")\n"
    << "    -R, --begRatio        FLOAT  Ratio of the contig size that is considered as the beginning part (default: " << Globals::beginning_ratio            << ", should be less than 0.5)\n"
    << "    -v, --minOverlap      INT    Minimum overlap between a molecule and a contig (default: "                   << Globals::min_overlap                << ")\n"
    << "    -D, --maxContigDist   INT    Merge contigs if they are separated by not more that N bp (default: "         << Globals::max_contig_distance        << ")\n"
    << "    -Z, --fillerSize      INT    Size of the stretch of Ns between the sequences (default: "                   << Globals::filler_size                << ")\n" 
    << "  * Log files\n"
    << "    -F, --scaffolds       FILE   Output scaffolds file name\n"
    << "    -g, --graph           FILE   Output in GFA format\n"
    << "    -p, --mapping         FILE   Output where the molecule map with respect to the contigs\n"
    << "    -C, --cisLinks        FILE   Output links between previous split contig parts\n"
    << "    -T, --transLinks      FILE   Output links between different contig parts\n";
}


void parse_parameters (int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      exit(EXIT_SUCCESS);
    } else if ((arg == "-i") || (arg == "--input")){
      Globals::input_file_name = argv[++i];
    } else if ((arg == "-f") || (arg == "--contigs")){
      Globals::contigs_file_name = argv[++i];
    } else if ((arg == "-m") || (arg == "--outputMolecule")){
      Globals::output_molecules_file_name = argv[++i];
    } else if ((arg == "-d") || (arg == "--distance")) {
      Globals::max_read_distance = std::stof(argv[++i]);
    } else if ((arg == "-q") || (arg == "--mapq")) {
      Globals::min_mapq = std::stoi(argv[++i]);
    } else if ((arg == "-Q") || (arg == "--mapqSolid")) {
      Globals::min_mapq_solid = std::stoi(argv[++i]);
    } else if ((arg == "-l") || (arg == "--lengthSolid")){
      Globals::min_len_solid = std::stoi(argv[++i]); 
    } else if ((arg == "-b") || (arg == "--barcodes")) {
      Globals::min_n_reads_barcode = std::stoi(argv[++i]);
    } else if ((arg == "-M") || (arg == "--inputMolecule")){
      Globals::input_molecules_file_name = argv[++i];
    } else if ((arg == "-w") || (arg == "--window")) {
      Globals::window = std::stoi(argv[++i]);
    } else if ((arg == "-s") || (arg == "--outputSplit")){
      Globals::output_split_file_name = argv[++i];
    } else if ((arg == "-t") || (arg == "--threshold")) {
      Globals::threshold = std::stof(argv[++i]);
    } else if ((arg == "-e") || (arg == "--sampleSize")){
      Globals::n_sample = std::stoi(argv[++i]);
    } else if ((arg == "-z") || (arg == "--minSize")){
      Globals::min_ctg_size = std::stoi(argv[++i]);
    } else if ((arg == "-E") || (arg == "--seed")){
      Globals::seed = std::stoi(argv[++i]);
    } else if ((arg == "-a") || (arg == "--counts")) {
      Globals::counts_file_name = argv[++i];
    } else if ((arg == "-A") || (arg == "--scores")) {
      Globals::scores_file_name = argv[++i];
    } else if ((arg == "-S") || (arg == "--inputSplit")){
      Globals::input_split_file_name = argv[++i];
    } else if ((arg == "-o") || (arg == "--output")){
      Globals::output_file_name = argv[++i];
    } else if ((arg == "-c") || (arg == "--arcsCondition")){
      Globals::condition = std::stoi(argv[++i]);
    } else if ((arg == "-r") || (arg == "--nReads")){
      Globals::min_n_reads = std::stoi(argv[++i]);
    } else if ((arg == "-R") || (arg == "--begRatio")){
      Globals::beginning_ratio = std::stof(argv[++i]);
    } else if ((arg == "-v") || (arg == "--minOverlap")){
      Globals::min_overlap = std::stoi(argv[++i]);
    } else if ((arg == "-Z") || (arg == "--fillerSize")){
      Globals::filler_size = std::atoi(argv[++i]);
    } else if ((arg == "-D") || (arg == "--maxContigDist")){
      Globals::max_contig_distance = std::stoi(argv[++i]);
    } else if ((arg == "-F") || (arg == "--scaffolds")){
      Globals::scaffold_file_name = argv[++i];
    } else if ((arg == "-g") || (arg == "--graph")){
      Globals::graph_file_name = argv[++i];
    } else if ((arg == "-p") || (arg == "--mapping")){
      Globals::mapping_file_name = argv[++i];
    } else if ((arg == "-C") || (arg == "--cisLinks")){
      Globals::cis_link_file_name = argv[++i];
    } else if ((arg == "-T") || (arg == "--transLinks")){
      Globals::trans_link_file_name = argv[++i];
    } else {
      std::cerr << "Error!  Parameter '" << argv[i] << "' is not understood.\nExiting.\n";
      exit(EXIT_FAILURE);
    }
  }
}
