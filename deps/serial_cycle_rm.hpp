#pragma once

#include <memory>
#include <stack>
#include <queue>

#include "bi_graph.hpp"

//d simplex traversal
class Bi_Graph_Traversal {
public:
    //associated bi-graph
    Bi_Graph* graphptr;
    
    std::vector<int> visit_flag;

    Bi_Graph_Traversal() = default;
    Bi_Graph_Traversal(Bi_Graph*);

    int findRootIterative(int);
    
    std::vector<int> getParent(int uidx);
    
    std::vector<int> getChild(int uidx);

    std::vector<int> getAncestor(int uidx);

    bool isBackwardAcyclic(std::vector<int>&, std::vector<int>&, int);

    int lookAheadDFS(std::queue<int>&, std::vector<int>&, int);

    int traversalBFS(std::queue<int>&, int root);

    int cycleRemoval(int);
};
