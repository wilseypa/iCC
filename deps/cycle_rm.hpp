#pragma once

#include <deque>

#include "bi_graph.hpp"

//d simplex traversal
class Bi_Graph_Traversal {
public:
    //associated bi-graph
    Bi_Graph* graphptr;
    
    std::vector<int> root_flag;

    Bi_Graph_Traversal() = default;
    Bi_Graph_Traversal(Bi_Graph*);

    int findRootIterative(int);
    
    std::vector<int> getParent(int uidx);
    
    std::vector<int> getChild(int uidx);

    std::vector<int> getAncestor(int rootnum, int uidx);

    bool isBackwardAcyclic(std::vector<int>&, std::vector<int>&);

    int lookAheadDFS(std::deque<int>&, std::vector<int>&, std::vector<int>&, int, int);

    int serialBFS(int, int);

    int serialCycleRemoval(int maxdegree);

    int threadBackwardBFS(int rootnum, int root);

    int parallelBFS(int rootnum, int root);

    int parallelRathod(int maxdegree);
};
