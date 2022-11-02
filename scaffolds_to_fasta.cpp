#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <experimental/filesystem>
#include <map>
#include <stdexcept>      // std::out_of_range

const std::string scaff_gap = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";


std::string complement_fast(std::string &DNAseq){

	reverse(DNAseq.begin(), DNAseq.end());
	for (size_t i = 0; i < DNAseq.length(); ++i){
        switch (DNAseq[i]){
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
        }
    }

	return DNAseq;

}
std::string complement(std::string &dna){
    std::string rev="";

	for(int i=1;i<=dna.length();i++){

       char nucleotide=dna[dna.length()-i];

       switch (nucleotide){

		   case 'A': rev=rev+"T"; break;
		   case 'C': rev=rev+"G"; break;
		   case 'G': rev=rev+"C"; break;
		   case 'T': rev=rev+"A"; break;
       }

	   if(i%100000 == 0) std::cout <<i<<std::endl;

     }
	 return rev;
}

int main (int argc, char* argv[]){

    std::string contig_fasta = argv[1];
    std::string scaffolds = argv[2];
    std::string scaffolds_fasta = argv[3];
	
    std::ifstream ctg_fasta(contig_fasta.c_str());
    std::ifstream scaff(scaffolds.c_str());
    std::ofstream scaff_fasta(scaffolds_fasta.c_str(), std::ofstream::out);	

    std::map<std::string,std::string> contigs;

    std::string ctg_name, ctg_seq = "";
	while(getline(ctg_fasta,ctg_name)){

		getline(ctg_fasta,ctg_seq);
		contigs.insert(std::pair<std::string,std::string>(ctg_name.substr(1),ctg_seq));

	}

    std::string scaffold_composition;
	int scaff_number = 0;
    std::string scaff_seq = "";
	while(getline(scaff,scaffold_composition)){

		scaff_number++;

        std::istringstream ss(scaffold_composition);
        std::string contig_complete;
		char sens;
        std::string contig;

		while(getline(ss, contig_complete, ';')) {
			
			if(contig_complete.length() > 0){

				sens = contig_complete[0];
				contig = contig_complete.substr(1);
                std::cout << sens <<" "<< contig << std::endl;
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
          return 1;
        }

				if(scaff_seq.length() > 0){

					scaff_seq.append(scaff_gap);
				}

				if(sens == '-')
					scaff_seq.append(complement_fast(ctg_seq));
				else 
					scaff_seq.append(ctg_seq);

			}

		}



		scaff_fasta << ">scaff_" << scaff_number << std::endl;
		scaff_fasta << scaff_seq << std::endl;
		scaff_seq = "";



	}


	scaff_fasta.close();
	scaff.close();
	ctg_fasta.close();


	 return 0;
}
