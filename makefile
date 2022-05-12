molconcat: molconcat.o Graph.o Contig.o
	g++ -std=c++17 -o joinASM molconcat.o Graph.o Contig.o -lstdc++fs -pthread

molconcat.o: molconcat.cpp Graph.h Contig.h
	g++ -std=c++17 -c molconcat.cpp -lstdc++fs -pthread

Graph.o: Graph.cpp Graph.h Contig.h
	g++ -std=c++17 -c Graph.cpp

Contig.o: Contig.cpp Contig.h Molecule.h
	g++ -std=c++17 -c Contig.cpp

Molecule.o: Molecule.cpp Molecule.h
	g++ -std=c++17 -c Molecule.cpp

clean:
	rm -f *~ *.o
