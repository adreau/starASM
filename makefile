all: scaffolds_to_fasta molconcat

scaffolds_to_fasta: scaffolds_to_fasta.cpp
	g++ -std=c++11 -o scaffolds_to_fasta scaffolds_to_fasta.cpp

molconcat: molconcat.o Contig.o
	g++ -std=c++11 -o joinASM molconcat.o Contig.o

molconcat.o: molconcat.cpp Contig.h Globals.h
	g++ -std=c++11 -c molconcat.cpp

Contig.o: Contig.cpp Contig.h Molecule.h Globals.h
	g++ -std=c++11 -c Contig.cpp

Molecule.o: Molecule.cpp Molecule.h
	g++ -std=c++11 -c Molecule.cpp

clean:
	rm -f *~ *.o
