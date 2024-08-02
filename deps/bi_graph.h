#ifndef BI_GRAPH_H
#define BI_GRAPH_H

#include <vector>

struct Bi_Graph{
    private:
    
    void randomBiGraphGen();

    public:

    int u, v, maxdegree;
    std::vector<std::vector<int>> adj_list;
    std::vector<int> match;

    Bi_Graph();
    Bi_Graph(int leftnum, int rightnum, int degree, bool random);
    void printBiGraph();
};

#endif
