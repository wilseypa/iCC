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

    int findRoot();
    
    int getParent(std::vector<int>& parent_workspace, int uidx);
    
    int getChild(std::vector<int>& child_workspace, int uidx);

    void getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx);

    bool isBackwardAcyclic(std::vector<int>&, std::vector<int>&);

    int lookAheadDFS(std::deque<int>&, std::vector<int>&, std::vector<int>&, int, int);

    int serialBFS(std::vector<int>&, int, int);

    int serialCycleRemoval(int maxdegree);

    int threadBackwardBFS(int rootnum, int uparent, int uidx);

    int parallelBFS(int rootnum, int root);

    int parallelRathod(int maxdegree);
};
