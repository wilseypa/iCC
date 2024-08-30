#pragma once

#include <memory>
#include <stack>
#include <queue>
#include <atomic>

#include "bi_graph.hpp"

// d simplex traversal
class Bi_Graph_Traversal
{
public:
    // associated bi-graph
    Bi_Graph *graphptr;

    std::vector<int> root_flag;

    Bi_Graph_Traversal() = default;
    Bi_Graph_Traversal(Bi_Graph *);

    int findRootIterative(int);

    std::vector<int> getParent(int uidx);

    std::vector<int> getChild(int uidx);

    std::vector<int> getAncestor(int uidx);

    bool isBackwardAcyclic(std::vector<int> &, std::vector<int> &, int);

    int lookAheadDFS(std::queue<int> &, std::vector<int> &, int);

    int serialBFS(std::queue<int> &, int);

    int serialCycleRemoval(int maxdegree);

    void threadBackwardDFS(int &, int);

    int parallelBFS(int, int);

    int parallelRathod(int maxdegree, int threadnum);
};
