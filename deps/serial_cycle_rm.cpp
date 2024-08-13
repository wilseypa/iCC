#include <memory>
#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <stack>
#include <iostream>

#include "bi_graph.hpp"
#include "serial_cycle_rm.hpp"

std::shared_ptr<Path_Node> Path_Node::getSharedPtr() {
    return shared_from_this();
}


//backward acyclicity check
bool Path_Node::acyclicBackward(Bi_Graph* bi_graph, int uidx) {
    //check acyclicity within path node
    for (auto& i: this->forward_simp_path) {
        int imate = bi_graph->match[i];
        //check if imate is part of uidx's facet edges
        auto iter = std::find(bi_graph->adj_list[uidx].begin(), bi_graph->adj_list[uidx].end(), imate);
        //found cycle
        if (iter != bi_graph->adj_list[uidx].end()) {
            return false;
        }
    }
    std::weak_ptr<Path_Node> weak_par_ptr = this->parent_ptr;
    while (!weak_par_ptr.expired()) {
        std::shared_ptr<Path_Node> shared_par_ptr = weak_par_ptr.lock();
        for (auto& i: shared_par_ptr->forward_simp_path) {
            int imate = bi_graph->match[i];
            //check if imate is part of uidx's facet edges
            auto iter = std::find(bi_graph->adj_list[uidx].begin(), bi_graph->adj_list[uidx].end(), imate);
            //found cycle
            if (iter != std::end(bi_graph->adj_list[uidx])) {
                return false;
            }
        }
        weak_par_ptr = shared_par_ptr->parent_ptr;
    }
    return true;
}


int findRoot(Bi_Graph* bi_graph, std::vector<int>& visit_flag) {
    //pick the smallest indexed u as the start (root)
    auto root_lambda = [&](const auto iter) {return (bi_graph->match[iter] > 0 && visit_flag[iter] == 0);};
    //range does not work on my machine at this time
    // auto view = std::ranges::iota_view(0, u);
    int u = visit_flag.size();
    std::vector<int> root_range(u);
    std::iota(root_range.begin(), root_range.end(), 0);
    auto iter = std::find_if(root_range.begin(), root_range.end(), root_lambda);
    if (iter != std::end(root_range)) {
        return *iter;
    }
    else {
        return -1;
    }
}

int pathDFS(Bi_Graph* bi_graph, std::stack<std::shared_ptr<Path_Node>>& dfs_stack, std::vector<int>& visit_flag, std::vector<int>& root_d_simp) {
    int revertedmatch = 0;
    std::vector<int> next_d_simp;
    while (!dfs_stack.empty()) {
        auto currentptr = dfs_stack.top();
        dfs_stack.pop();
        
        int branch = 1;

        //check number of branch
        next_d_simp.clear();
        next_d_simp.push_back(currentptr->uidx);
        //no fork
        while (branch == 1) {

            // std::cout << "current d simp in branch while loop = " << next_d_simp.back() << '\n';

            //if acyclic, add uidx(next_d_simp.back()) to forward_simp_path
            if (currentptr->acyclicBackward(bi_graph, next_d_simp.back())) {
                currentptr->forward_simp_path.push_back(next_d_simp.back());
            } else {
                //remove the match
                bi_graph->match[next_d_simp.back()] = -1;
                revertedmatch += 1;
            }
            visit_flag[next_d_simp.back()] = 1;

            //find next d simp, check the neighbor of next_d_simp.back()
            branch = 0;
            for (const auto& i: bi_graph->adj_list[next_d_simp.back()]) {
                int imate = bi_graph->match[i];
                if (imate != -1 && imate != next_d_simp.back()) {
                    //check if imate is an existing root node
                    if (std::find(root_d_simp.begin(), root_d_simp.end(), imate) == root_d_simp.end()) {
                        branch += 1;
                        next_d_simp.push_back(imate);
                        visit_flag[imate] = 1;
                    }
                }
            }
        }
        //branch == 0 or branch > 1, push the new branch to the queue
        for (int i = 0; i < branch; i++) {
            dfs_stack.push(std::make_shared<Path_Node>(next_d_simp.rbegin()[i], currentptr));
        }

        // std::cout << "branch in outer while iter = " << branch << '\n';
        // std::cout << "next d simp size = " << next_d_simp.size() << '\n';
        // std::cout << "stack size  = " << dfs_stack.size() << '\n';

    }
    return revertedmatch;
}


int serialCycleRemove(Bi_Graph* bi_graph) {
    //assume u is the d simplex
    int u = bi_graph->u;

    int revertedmatch = 0;
    
    std::vector<int> visit_flag(u, 0);
    std::vector<int> root_d_simp;

    int root = findRoot(bi_graph, visit_flag);
    std::cout << "first root = " << root << '\n';

    while (!(root < 0)) {
        
        std::cout << "root = " << root << '\n';

        root_d_simp.push_back(root);

        std::stack<std::shared_ptr<Path_Node>> dfs_stack;
        dfs_stack.push(std::make_shared<Path_Node>(root));
        visit_flag[root] = 1;
        
        revertedmatch += pathDFS(bi_graph, dfs_stack, visit_flag, root_d_simp);

        std::cout << "reverted match = " << revertedmatch << '\n';
        
        root = findRoot(bi_graph, visit_flag);
    }
    return revertedmatch;
}