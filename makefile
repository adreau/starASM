SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

all: starASM

starASM: $(OBJS)
	g++ -g -fsanitize=address -std=c++11 -Wall -pedantic -O3 -lhts -o starASM $(OBJS)

%.o: %.cpp
	g++ -g -fsanitize=address -std=c++11 -Wall -pedantic -O3 -c $< -o $@

clean:
	rm -f *~ *.o
