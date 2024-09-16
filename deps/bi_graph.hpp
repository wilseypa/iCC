#pragma once

#include <vector>
#include <atomic>

struct Bi_Graph
{
private:
    const int MAX_PRINT_SIZE = 20;
    void randomBiGraphGen();

public:
    int u, v, maxdegree;
    std::vector<std::vector<int>> adj_list;
    std::vector<int> match;

    Bi_Graph() = default;
    Bi_Graph(int, int, int, bool);
    Bi_Graph(int, int, int);
    void addEdge(int, int);
    void addEdge(int, std::vector<int> &);
    void printBiGraph();
    int parallelKarpSipserInit(int);
    int parallelDFSMatch(int);
    int dfsAugPath(int, std::vector<int> &, std::vector<int> &, std::vector<int> &);
    void findMatch(int, std::vector<std::atomic<int>> &, std::vector<std::atomic<int>> &);
};