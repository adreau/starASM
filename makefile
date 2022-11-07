all: scaffolds_to_fasta joinASM

scaffolds_to_fasta: scaffolds_to_fasta.cpp
	g++ -std=c++11 -Wall -o scaffolds_to_fasta scaffolds_to_fasta.cpp

joinASM: joinASM.o
	g++ -std=c++11 -Wall -o joinASM joinASM.o

joinASM.o: joinASM.cpp Contig.h Globals.h Graph.h Scaffold.h
	g++ -std=c++11 -Wall -c joinASM.cpp

clean:
	rm -f *~ *.o
