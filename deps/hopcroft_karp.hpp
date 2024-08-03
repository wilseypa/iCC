#pragma once

#include <set>
#include <queue>
#include <climits>

class HKGraph
{
private:
    int m;
    int n;

    std::vector<std::vector<int>> adj;
    std::vector<int> pair_u;
    std::vector<int> pair_v;
    std::vector<int> dist;

    bool bfs();
    bool dfs(int u);
public:
    HKGraph() = default;
    HKGraph(int m, int n);
    void addEdge(int u, int v);
    std::vector<std::pair<int, int>> hopcroftKarpAlgorithm();
    std::pair<std::set<int>, std::set<int>> custom_hopcroftKarpAlgorithm();
};
