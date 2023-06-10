SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

all: starASM

starASM: $(OBJS) ./htslib/libhts.so
	g++ -g -fsanitize=address -static-libasan -std=c++11 -Wall -pedantic -O3 -Wl,-rpath,"./htslib/" -o starASM $(OBJS) ./htslib/libhts.so

./htslib/libhts.so:
	cd htslib && make

%.o: %.cpp
	g++ -g -fsanitize=address -std=c++11 -Wall -pedantic -O3 -c $< -o $@

clean:
	rm -f *~ *.o
