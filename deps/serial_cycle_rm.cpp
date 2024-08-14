#include <memory>
#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <stack>
#include <queue>

#include <iostream>

#include "bi_graph.hpp"
#include "serial_cycle_rm.hpp"

std::shared_ptr<Path_Node> Path_Node::getSharedPtr() {
    return shared_from_this();
}


//backward acyclicity check
bool Path_Node::acyclicBackward(Bi_Graph* bi_graph, int uidx) {
    //check if uidx can reach one of the previous root
    //check acyclicity within path node
    for (auto& i: this->forward_simp_path) {
        int imate = bi_graph->match[i];
        //check if imate is part of uidx's facet edges
        auto iter = std::find(bi_graph->adj_list[uidx].begin(), bi_graph->adj_list[uidx].end(), imate);
        //found cycle
        if (iter != std::end(bi_graph->adj_list[uidx])) {
            return false;
        }
    }

    //check ancestors
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


int findRoot(Bi_Graph* bi_graph, std::vector<int>& visit_flag, std::vector<int>& root_d_simp) {
    int u = bi_graph->u;
    int v = bi_graph->v;
    //check the number of matched neighbors of d-1 simp
    std::vector<int> v_range{v};
    std::iota(v_range.begin(), v_range.end(), u);
    auto match_lambda = [&](const auto idx) {return(bi_graph->match[idx] != -1);};
    //pragma omp parallel for 
    for (auto& v: v_range) {
        int vmatch = std::count_if(std::begin(bi_graph->adj_list[v]), std::end(bi_graph->adj_list[v]), match_lambda);
        //keep the v with vmatch == 1 as the root
        if (vmatch == 1) {
            //pragma omp critical
            root_d_simp.push_back(bi_graph->match[v]);
        }
    }
    //if no proper root, pick the smallest matched d simp and cut the incoming path of that d simp
    if (root_d_simp.size() == 0) {
        int root = -1;
        for (auto& v: v_range) {
            if (bi_graph->match[v] != -1) {
                root = bi_graph->match[v];
                break;
            }
        }
        int rootmate = bi_graph->match[root];
        for (auto& u: bi_graph->adj_list[rootmate]) {
            if (bi_graph->match[u] != -1 && u != root) {
                int umate = bi_graph->match[u];
                bi_graph->match[u] = -1;
                bi_graph->match[umate] = -1;
            }
        }
        root_d_simp.push_back(root);
    }
    return root_d_simp.size();
}


int pathDFS(Bi_Graph* bi_graph, std::stack<std::shared_ptr<Path_Node>>& dfs_stack, std::vector<int>& visit_flag, std::vector<int>& root_d_simp) {
    int revertedmatch = 0;
    //work space for path node construction
    std::vector<int> next_d_simp;

    while (!dfs_stack.empty()) {
        auto currentptr = dfs_stack.top();
        dfs_stack.pop();

        //check number of branch
        next_d_simp.clear();
        next_d_simp.push_back(currentptr->uidx);
        std::cout << "the d simp before while true loop = " << next_d_simp.back() << '\n';
        int branch;
        //while no fork
        while (true) {

            std::cout << "current d simp in while true loop = " << next_d_simp.back() << '\n';

            //if acyclic, add uidx(next_d_simp.back()) to forward_simp_path
            if (currentptr->acyclicBackward(bi_graph, next_d_simp.back())) {
                currentptr->forward_simp_path.push_back(next_d_simp.back());
            } else {
                std::cout << "reverted simp = " << next_d_simp.back() << '\n';
                //remove the match
                int temp = bi_graph->match[next_d_simp.back()];
                bi_graph->match[next_d_simp.back()] = -1;
                bi_graph->match[temp] = -1;
                revertedmatch += 1;
            }
            visit_flag[next_d_simp.back()] = 1;

            //find next d simp, check the neighbor of next_d_simp.back()
            branch = 0;
            for (const auto& v: bi_graph->adj_list[next_d_simp.back()]) {
                int vmate = bi_graph->match[v];
                if (vmate != -1 && visit_flag[vmate] == 0) {
                    
                    branch += 1;
                    next_d_simp.push_back(vmate);
                    visit_flag[vmate] = 1;
                }
            }

            std::cout << "the branch # after current d simp in while true loop = " << branch << '\n';

            if (branch != 1) break;
        }
        //branch == 0 or branch > 1, push the new branch to the stack/queue
        for (int i = 0; i < branch; i++) {
            dfs_stack.push(std::make_shared<Path_Node>(next_d_simp.rbegin()[i], currentptr));
        }
    }
    return revertedmatch;
}


int serialDFSCycleRemove(Bi_Graph* bi_graph) {
    //assume u is the d simplex
    int u = bi_graph->u;

    int revertedmatch = 0;
    
    std::vector<int> visit_flag(u, 0);
    std::vector<int> root_d_simp;
    std::stack<std::shared_ptr<Path_Node>> dfs_stack;

    int rootnum = findRoot(bi_graph, visit_flag, root_d_simp);
    std::cout << "number of root = " << rootnum << '\n';

    while (!(root_d_simp.empty())) {

        int root = root_d_simp.back();
        
        std::cout << "current root = " << root << '\n';

        dfs_stack.push(std::make_shared<Path_Node>(root));
        visit_flag[root] = 1;
        
        revertedmatch += pathDFS(bi_graph, dfs_stack, visit_flag, root_d_simp);
        
        root_d_simp.pop_back();
    }
    return revertedmatch;
}

int pathBFS(Bi_Graph* bi_graph, std::queue<std::shared_ptr<Path_Node>>& bfs_queue, std::vector<int>& visit_flag, std::vector<int>& root_d_simp) {
    int revertedmatch = 0;
    //work space for path node construction
    std::vector<int> next_d_simp;

    while (!bfs_queue.empty()) {
        auto currentptr = bfs_queue.front();
        bfs_queue.pop();

        //check number of branch
        next_d_simp.clear();
        next_d_simp.push_back(currentptr->uidx);
        std::cout << "the d simp before while true loop = " << next_d_simp.back() << '\n';
        int branch;
        //while no fork
        while (true) {

            std::cout << "current d simp in while true loop = " << next_d_simp.back() << '\n';

            //if acyclic, add uidx(next_d_simp.back()) to forward_simp_path
            if (currentptr->acyclicBackward(bi_graph, next_d_simp.back())) {
                currentptr->forward_simp_path.push_back(next_d_simp.back());
            } else {
                std::cout << "reverted simp = " << next_d_simp.back() << '\n';
                //remove the match
                int temp = bi_graph->match[next_d_simp.back()];
                bi_graph->match[next_d_simp.back()] = -1;
                bi_graph->match[temp] = -1;
                revertedmatch += 1;
            }
            visit_flag[next_d_simp.back()] = 1;

            //find next d simp, check the neighbor of next_d_simp.back()
            branch = 0;
            for (const auto& v: bi_graph->adj_list[next_d_simp.back()]) {
                int vmate = bi_graph->match[v];
                if (vmate != -1 && visit_flag[vmate] == 0) {
                    
                    branch += 1;
                    next_d_simp.push_back(vmate);
                    visit_flag[vmate] = 1;
                }
            }

            std::cout << "the branch # after current d simp in while true loop = " << branch << '\n';

            if (branch != 1) break;
        }
        //branch == 0 or branch > 1, push the new branch to the stack/queue
        for (int i = 0; i < branch; i++) {
            bfs_queue.push(std::make_shared<Path_Node>(next_d_simp.rbegin()[i], currentptr));
        }
    }
    return revertedmatch;
}

int serialBFSCycleRemove(Bi_Graph* bi_graph) {
    //assume u is the d simplex
    int u = bi_graph->u;

    int revertedmatch = 0;
    
    std::vector<int> visit_flag(u, 0);
    std::vector<int> root_d_simp;
    std::queue<std::shared_ptr<Path_Node>> bfs_queue;

    int rootnum = findRoot(bi_graph, visit_flag, root_d_simp);
    std::cout << "number of root = " << rootnum << '\n';

    while (!(root_d_simp.empty())) {

        int root = root_d_simp.back();
        
        std::cout << "current root = " << root << '\n';

        bfs_queue.push(std::make_shared<Path_Node>(root));
        visit_flag[root] = 1;
        
        revertedmatch += pathBFS(bi_graph, bfs_queue, visit_flag, root_d_simp);
        
        root_d_simp.pop_back();
    }
    return revertedmatch;
}