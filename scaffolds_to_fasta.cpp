#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <stdexcept>      // std::out_of_range

const size_t line_len = 60;

static void show_usage(char *name) {

  std::cerr << "Usage: " << name << "\n"
    << "Options:\n"
    << "\t-h, --help            Show this help message\n"
    << "\t-c, --contigs   FILE  Input contigs FASTA file\n"
    << "\t-j, --joins     FILE  Input join file, provided by joinASM\n"
    << "\t-s, --scaffolds FILE  Output FASTA file\n"
    << "\t-s, --filerSize INT   Size of the stretch of Ns between the sequences (default: 100)\n";
}


void complement (std::string &DNAseq){

    reverse(DNAseq.begin(), DNAseq.end());
    for (size_t i = 0; i < DNAseq.length(); ++i){
        switch (DNAseq[i]) {
            case 'A':
                DNAseq[i] = 'T';
                break;    
            case 'C':
                DNAseq[i] = 'G';
                break;
            case 'G':
                DNAseq[i] = 'C';
                break;
            case 'T':
                DNAseq[i] = 'A';
                break;
            default:
                std::cerr << "Error! Unknown nucleotide '" << DNAseq[i] << "'.";
                exit(EXIT_FAILURE);
        }
    }
}

void read_fasta (std::map < std::string, std::string > &contigs, std::string &contig_fasta) {

    std::ifstream ctg_fasta(contig_fasta.c_str());
    if (! ctg_fasta.is_open()) {
        std::cerr << "Error!  Cannot open file '" << contig_fasta << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string ctg_line, ctg_name, ctg_seq;

    while (getline(ctg_fasta, ctg_line)) {
        if (! ctg_line.empty()) {
            if (ctg_line[0] == '>') {
                if (! ctg_seq.empty()) {
                    contigs[ctg_name] = ctg_seq;
                }
                ctg_name = ctg_line.substr(1);
                ctg_seq.clear();
            }
            else {
                ctg_seq += ctg_line;
            }
        }
    }
    if (! ctg_seq.empty()) {
        contigs[ctg_name] = ctg_seq;
    }

    ctg_fasta.close();
    std::cout << "Read " << contigs.size() << " contigs.\n";
}

void write_fasta_sequence (int id, std::string &sequence, std::ofstream &file) {
    file << ">scaff_" << id << "\n";
    for (size_t i = 0; i <= sequence.length() / line_len; ++i) {
        file << sequence.substr(i * line_len, line_len) << "\n";
    }
}

int main (int argc, char* argv[]){

    if (argc < 1) {
        show_usage(argv[0]);
        return EXIT_FAILURE;
    }

    std::string contigs_file_name;
    std::string joins_file_name;
    std::string scaffolds_file_name;
    size_t filler_size;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return EXIT_SUCCESS;
        } else if ((arg == "-c") || (arg == "--contigs")){
            contigs_file_name = argv[++i];
        } else if ((arg == "-j") || (arg == "--joins")){
            joins_file_name = argv[++i];
        } else if ((arg == "-s") || (arg == "--scaffolds")){
            scaffolds_file_name = argv[++i];
        } else if ((arg == "-f") || (arg == "--fillerSize")){
            filler_size = std::atoi(argv[++i]);
        } else {
            show_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    std::ifstream scaff(joins_file_name.c_str());
    std::ofstream scaff_fasta(scaffolds_file_name.c_str(), std::ofstream::out);	
    if (! scaff.is_open()) {
        std::cerr << "Error!  Cannot open file '" << joins_file_name << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (! scaff_fasta.is_open()) {
        std::cerr << "Error!  Cannot open file '" << scaffolds_file_name << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::map < std::string, std::string > contigs;
    read_fasta(contigs, contigs_file_name);

    std::string ctg_seq;
    std::string scaffold_composition;
    std::string scaff_seq;
    std::string scaff_gap('N', filler_size);

    unsigned int n_scaffolds;
    for (n_scaffolds = 1; getline(scaff, scaffold_composition); ++n_scaffolds) {

        std::istringstream ss(scaffold_composition);
        std::string contig_complete;
        char strand;
        std::string contig;
        unsigned long int start, end;

        while (getline(ss, contig_complete, ';')) {

            if (! contig_complete.empty()) {

                std::stringstream splitstream (contig_complete);
                splitstream >> strand >> contig >> start >> end;

                // C++ is 0-based
                --start;
                --end;

                try {
                    ctg_seq = contigs.at(contig);
                }
                catch (const std::out_of_range& oor) {
                    std::cerr << "Cannot find contig '" << contig << "'.\nContigs should be in (first ten)";
                    unsigned int i = 0;
                    for (auto &c: contigs) {
                        std::cerr << " '" << c.first << "'";
                        ++i;
                        if (i >= 10) break;
                    }
                    std::cerr << "...\nExiting now.\n";
                    exit(EXIT_FAILURE);
                }

                if (! scaff_seq.empty()) {
                    scaff_seq.append(scaff_gap);
                }

                ctg_seq.substr(start, end - start + 1);
                if (strand == '-') {
                    complement(ctg_seq);
                }
                scaff_seq.append(ctg_seq);
            }
        }

        write_fasta_sequence(n_scaffolds + 1, scaff_seq, scaff_fasta);
        scaff_seq.clear();

        if (n_scaffolds % 100 == 0) {
            std::cout << n_scaffolds << " lines read.\r" << std::flush;
        }
    }
    std::cout << n_scaffolds << " lines read, done.\n";

    scaff_fasta.close();
    scaff.close();

    return 0;
}
