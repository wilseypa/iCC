CC = g++
CFLAGS = -g -std=c++20 -fopenmp -O2

all: match

match: main.cpp bi_graph.cpp cycle_rm.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) main.cpp bi_graph.cpp cycle_rm.cpp -o match

clean: 
	rm -f match
