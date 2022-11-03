all: scaffolds_to_fasta molconcat

scaffolds_to_fasta: scaffolds_to_fasta.cpp
	g++ -std=c++17 -o scaffolds_to_fasta scaffolds_to_fasta.cpp -lstdc++fs -pthread

molconcat: molconcat.o Contig.o
	g++ -std=c++17 -o joinASM molconcat.o Contig.o -lstdc++fs -pthread

molconcat.o: molconcat.cpp Contig.h
	g++ -std=c++17 -c molconcat.cpp -lstdc++fs -pthread

Contig.o: Contig.cpp Contig.h Molecule.h
	g++ -std=c++17 -c Contig.cpp

Molecule.o: Molecule.cpp Molecule.h
	g++ -std=c++17 -c Molecule.cpp

clean:
	rm -f *~ *.o
