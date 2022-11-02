#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <experimental/filesystem>
#include <map>
#include <stdexcept>      // std::out_of_range

using namespace std;
const string scaff_gap = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";


string complement_fast(string &DNAseq){

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
string complement(string &dna){
	string rev="";

	for(int i=1;i<=dna.length();i++){

       char nucleotide=dna[dna.length()-i];

       switch (nucleotide){

		   case 'A': rev=rev+"T"; break;
		   case 'C': rev=rev+"G"; break;
		   case 'G': rev=rev+"C"; break;
		   case 'T': rev=rev+"A"; break;
       }

	   if(i%100000 == 0) cout <<i<<endl;

     }
	 return rev;
}

int main (int argc, char* argv[]){

	string contig_fasta = argv[1];
	string scaffolds = argv[2];
	string scaffolds_fasta = argv[3];
	
	ifstream ctg_fasta(contig_fasta.c_str());
	ifstream scaff(scaffolds.c_str());
	ofstream scaff_fasta(scaffolds_fasta.c_str(), ofstream::out);	

	map<string,string> contigs;

	string ctg_name, ctg_seq = "";
	while(getline(ctg_fasta,ctg_name)){

		getline(ctg_fasta,ctg_seq);
		contigs.insert(pair<string,string>(ctg_name.substr(1),ctg_seq));

	}

	string scaffold_composition;
	int scaff_number = 0;
	string scaff_seq = "";
	while(getline(scaff,scaffold_composition)){

		scaff_number++;

		istringstream ss(scaffold_composition);
		string contig_complete;
		char sens;
		string contig;

		while(getline(ss, contig_complete, ';')) {
			
			if(contig_complete.length() > 0){

				sens = contig_complete[0];
				contig = contig_complete.substr(1);
				cout << sens <<" "<< contig << endl;
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



		scaff_fasta << ">scaff_" << scaff_number << endl;
		scaff_fasta << scaff_seq << endl;
		scaff_seq = "";



	}


	scaff_fasta.close();
	scaff.close();
	ctg_fasta.close();


	 return 0;
}
