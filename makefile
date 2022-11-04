all: scaffolds_to_fasta molconcat

scaffolds_to_fasta: scaffolds_to_fasta.cpp
	g++ -std=c++11 -Wall -o scaffolds_to_fasta scaffolds_to_fasta.cpp

molconcat: molconcat.o Contig.o
	g++ -std=c++11 -Wall -o joinASM molconcat.o Contig.o

molconcat.o: molconcat.cpp Contig.h Globals.h
	g++ -std=c++11 -Wall -c molconcat.cpp

Contig.o: Contig.cpp Contig.h Globals.h
	g++ -std=c++11 -Wall -c Contig.cpp

clean:
	rm -f *~ *.o
