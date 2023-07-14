SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

all: starASM

starASM: $(OBJS) ./htslib/libhts.so
	g++ -g -fsanitize=address -static-libasan -std=c++11 -Wall -pedantic -O2 -Wl,-rpath,$(shell pwd)/htslib -o starASM $(OBJS) ./htslib/libhts.so

./htslib/libhts.so:
	cd htslib && make

%.o: %.cpp
	g++ -g -fsanitize=address -std=c++11 -Wall -pedantic -O2 -c $< -o $@

clean:
	rm -f *~ *.o
