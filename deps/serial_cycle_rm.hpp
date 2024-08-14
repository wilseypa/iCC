#pragma once

#include <memory>
#include <stack>
#include <queue>

#include "bi_graph.hpp"

//the d simplex on a non-forked forward path forms a super node
class Path_Node : public std::enable_shared_from_this<Path_Node>
{
public:
    //the d-simplex on the path
    std::vector<int> forward_simp_path;
    //init d-simplex index in the path node
    int uidx;
    //construction status of the path
    // bool finished = false;

    std::shared_ptr<Path_Node> parent_ptr;
    // std::vector<std::shared_ptr<Path_Node>> children_vec;

    Path_Node(int u): uidx(u) {
    }

    Path_Node(int u, std::shared_ptr<Path_Node>& par_ptr): uidx(u), parent_ptr(par_ptr) {
    }

    std::shared_ptr<Path_Node> getSharedPtr();
    bool acyclicBackward(Bi_Graph* bi_graph, int uidx);
};

//need to add these functions to bi_graph class
int findRoot(Bi_Graph* bi_graph, std::vector<int>& visit_flag, std::vector<int>& root_d_simp);

int pathDFS(Bi_Graph* bi_graph, std::stack<std::shared_ptr<Path_Node>>& dfs_stack, std::vector<int>& visit_flag, std::vector<int>& root_d_simp);

int serialDFSCycleRemove(Bi_Graph* bi_graph);

int pathBFS(Bi_Graph* bi_graph, std::queue<std::shared_ptr<Path_Node>>& bfs_queue, std::vector<int>& visit_flag, std::vector<int>& root_d_simp);

int serialBFSCycleRemove(Bi_Graph* bi_graph);
