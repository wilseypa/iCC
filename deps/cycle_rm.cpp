#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <stack>
#include <deque>
#include <queue>
#include <set>
#include <execution>
#include <atomic>
#include <thread>

#include <iostream>

#include "bi_graph.hpp"
#include "cycle_rm.hpp"

Bi_Graph_Traversal::Bi_Graph_Traversal(Bi_Graph* bi_graph) : graphptr(bi_graph) {
    root_flag.resize(graphptr->u, 0);
}


int Bi_Graph_Traversal::findRootIterative(int maxdegree) {
    int u = graphptr->u;

    //check the number of matched neighbors of d-1 simp
    std::vector<int> u_range(u);
    std::iota(u_range.begin(), u_range.end(), 0);
    auto match_lambda = [&](const auto& i) {return(graphptr->match[i] != -1);};

    int root = -1;
    int deg = 1;
    while (deg <= maxdegree) {
        for (auto& i: u_range) {
            if (root_flag[i] != 0) continue;
            //par?
            int matchnum = std::count_if(std::begin(graphptr->adj_list[i]), std::end(graphptr->adj_list[i]), match_lambda);
            if (matchnum <= deg) {
                root = i;
                // std::cout<<"returned root = "<<root<<'\n';
                return root;
            }
        }
        deg += 1;
    }
    return root;
}


//directly return vec. use rvo
std::vector<int> Bi_Graph_Traversal::getParent(int uidx) {
    std::vector<int> parent_d_simp;
    int vidx = graphptr->match[uidx];
    //uidx is a source node has no parent (umatched)
    if (vidx == -1) return parent_d_simp;
     //i is d simp
    for (auto& i: graphptr->adj_list[vidx]) {
        int imate = graphptr->match[i];
        if (imate != vidx) parent_d_simp.push_back(i);
    }
    return parent_d_simp;
}

//no need to check visit flag. if found child in searched part. it might be merging or cycle
//in later stages, need to check visit flag when pushing to bfs queue
std::vector<int> Bi_Graph_Traversal::getChild(int uidx) {
    std::vector<int> child_d_simp;
    //i is d-1 simp
    for (auto& i: graphptr->adj_list[uidx]) {
        int imate = graphptr->match[i];
        if (imate != -1 && imate != uidx) child_d_simp.push_back(imate); 
    }

    return child_d_simp;
}

//ancestor in bfs searched part
std::vector<int> Bi_Graph_Traversal::getAncestor(int rootnum, int uidx) {
    //use rootnum to confine the ancestor search in current "bfs component"
    std::vector<int> ancestor_d_simp;
    std::queue<int> bfs_queue;

    //temp change uidx's rootflag to avoid cycle infi loop
    int rootflag = root_flag[uidx];
    root_flag[uidx] = -1;
    
    std::vector<int> uidx_parent = getParent(uidx);
    for (auto& i: uidx_parent) {
        if(root_flag[i] == rootnum) bfs_queue.push(i);
    }

    // std::cout<<"get ancestor started \n";
    //cycle can be found during this process but need ancestors for look up
    while (!bfs_queue.empty()) {
        int front = bfs_queue.front();
        ancestor_d_simp.push_back(front);
        bfs_queue.pop();

        std::vector<int> front_parent = getParent(front);
        for(auto& i: front_parent) {
            if (root_flag[i] == rootnum && std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) == ancestor_d_simp.end()) {
                bfs_queue.push(i);
            }
        }
    }

    //recover rootflag
    root_flag[uidx] = rootflag;

    return ancestor_d_simp;
}


bool Bi_Graph_Traversal::isBackwardAcyclic(std::vector<int>& ancestor_d_simp, std::vector<int>& u_child) {
    //par
    for (auto& i: u_child) {
        if (root_flag[i] == 0) continue;
        if (std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) != ancestor_d_simp.end()) {
            return false;
        }
    }
    return true;
}

int Bi_Graph_Traversal::lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_d_simp, std::vector<int>& child_d_simp, int rootnum, int uidx) {
    if (child_d_simp.size() == 0) return 0;

    std::vector<int> lookahead_vec;

    // auto push_lambda = [&](const auto i) {return (std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end());};

    for (auto& i: child_d_simp) {
        std::vector<int> i_parent = getParent(i);
        //if uidx is the only parent of i (in the whole graph). use look ahead shortcut
        if (i_parent.size() == 1) {
            lookahead_vec.push_back(i);
        } else if (std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
            //if i has more than 1 par and has not been searched yet (or in current bfs component). push to global bfs queue
            graph_bfs_queue.push_back(i);
        }
    }

    //stat var
    // int maxsize = 0;

    //start look ahead op on children of uidx
    ancestor_d_simp.push_back(uidx);
    int reverted = 0;
    bool flag;    //working var
    //the nodes appeared in lookahead stack at the same time are independent 
    while (!lookahead_vec.empty()) {
        //for stat
        // if (lookahead_vec.size() > maxsize) maxsize = lookahead_vec.size();

        flag = false;

        int top = lookahead_vec.back();
        lookahead_vec.pop_back();
        root_flag[top] = rootnum;

        //std::cout<<"lookahead top = "<<top<<"  lookahead start = "<<uidx<<'\n';

        std::vector<int> top_child_d_simp = getChild(top);
        if (top_child_d_simp.size() == 0) continue;

        //child in unsearched part of current bfs component
        auto iter = std::remove_if(child_d_simp.begin(), child_d_simp.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});
        child_d_simp.erase(iter, child_d_simp.end());

        //if top is not acyclic. the lookahead of top is over
        //push its unsearched children or reencountered children to bfs queue
        if (!isBackwardAcyclic(ancestor_d_simp, top_child_d_simp)) {
            
            int temp = graphptr->match[top];
            graphptr->match[top] = -1;
            graphptr->match[temp] = -1;
            reverted += 1;

        } else {
        //check if the child of top can be pushed to look ahead stack
            for (auto& i: top_child_d_simp) {
                std::vector<int> i_parent = getParent(i);
                //if top is the only parent of i. keep look ahead op
                if (i_parent.size() == 1) {
                    lookahead_vec.push_back(i);
                    //raise flag for adding top to ancestor
                    flag = true;
                } else if (std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
                    graph_bfs_queue.push_back(i);
                } 
            }
        }
        //if any of top's child is pushed to look ahead stack, update ancestor
        if (flag) ancestor_d_simp.push_back(top);
        //not all nodes in the stack have ances-des relationship with top
        //but this does not bring extra cycle rm. simplfied implementation
    }

    // std::cout<<"max LA stack size = "<<maxsize<<'\n';

    return reverted;
 }


int Bi_Graph_Traversal::serialBFS(int rootnum, int root) {
    int reverted = 0;

    std::deque<int> bfs_queue;
    bfs_queue.push_back(root);

    // std::cout<<"started at root = "<<root<<'\n';

    while (!bfs_queue.empty()) {
        int front = bfs_queue.front();
        bfs_queue.pop_front();
        root_flag[front] = rootnum;
        
        std::vector<int> front_child = getChild(front);
        if (front_child.size() == 0) continue;

        //child in unsearched part of current bfs component
        auto iter = std::remove_if(front_child.begin(), front_child.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});
        front_child.erase(iter, front_child.end());

        std::vector<int> ancestor_d_simp = getAncestor(rootnum, front);
        
        if (!isBackwardAcyclic(ancestor_d_simp, front_child)) {
            int temp = graphptr->match[front];
            graphptr->match[front] = -1;
            graphptr->match[temp] = -1;
            reverted += 1;
        } else {
            //look ahead and push descendant to bfs queue
            reverted += lookAheadDFS(bfs_queue, ancestor_d_simp, front_child, rootnum, front);
        }

        // std::cout<<"queue front = "<<front<<'\n';

    }
    
    // std::cout<<"ended at root = "<<root<<'\n';


    return reverted;
}



int Bi_Graph_Traversal::serialCycleRemoval(int maxdegree) {
    //assume u is the d simplex of d interface
    // int u = graphptr->u;

    int reverted = 0;

    int root = findRootIterative(maxdegree);
    int rootnum = 0;

    while (root != -1) {

        std::cout << "current root = " << root << '\n';

        rootnum += 1;

        reverted += serialBFS(rootnum, root);

        root = findRootIterative(maxdegree);
    }
    return reverted;
}


int Bi_Graph_Traversal::threadBackwardBFS(int rootnum, int uidx) {

    int reverted = 0;

    //get ancestor  //backward bfs
    std::vector<int> ancestor_d_simp;
    std::queue<int> bfs_queue;
    
    std::vector<int> uidx_parent = getParent(uidx);
    for (auto& i: uidx_parent) {
        if(root_flag[i] == rootnum) bfs_queue.push(i);
    }

    while (!bfs_queue.empty()) {
        int front = bfs_queue.front();
        ancestor_d_simp.push_back(front);
        bfs_queue.pop();

        std::vector<int> front_parent = getParent(front);
        for(auto& i: front_parent) {
            //found cycle
            if (i == uidx) {
                int temp = graphptr->match[uidx];
                graphptr->match[uidx] = -1;
                graphptr->match[temp] = -1;
                reverted += 1;
                return reverted;
            }
            if (root_flag[i] == rootnum && std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) == ancestor_d_simp.end()) {
                bfs_queue.push(i);
            }
        }
    }

    return reverted;
}

//Rathod Alg 1 procedure BFSComponent
int Bi_Graph_Traversal::parallelBFS(int rootnum, int root) {

    int reverted = 0;

    std::deque<int> graph_bfs_queue;
    graph_bfs_queue.push_back(root);

    auto push_lambda = [&](const auto i) {return (graphptr->match[i] != -1 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end());};

    while (!graph_bfs_queue.empty()) {
        int front = graph_bfs_queue.front();
        graph_bfs_queue.pop_front();
        root_flag[front] = rootnum;

        //std::cout<<"root = " <<root<<"  front = "<<front<<'\n';
        
        //get "leading up-edges"
        std::vector<int> front_child = getChild(front);
        //keep the unvisited or reencountered
        //remove_if only logically remove elements, read the docs!!!

        auto iter = std::remove_if(front_child.begin(), front_child.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});

        front_child.erase(iter, front_child.end());

        //use -1 to mark in-search
        // std::for_each(front_child.begin(), front_child.end(), [&](const auto& i) {root_flag[i] = -1;});

        //dynamic distribute work
#pragma omp parallel for schedule(dynamic) reduction(+ : reverted)
        for (int i = 0; i < front_child.size(); i++)
        {
            // std::cout<<"front = "<<front<<"  front child i= "<<front_child[i]<<'\n';
            reverted += threadBackwardBFS(rootnum, front_child[i]);
        }
        
        //back to the main thread, mark up and push acyclic leading up-edges to the queue 
        std::copy_if(front_child.begin(), front_child.end(), std::back_inserter(graph_bfs_queue), push_lambda);
    }

    return reverted;
}

int Bi_Graph_Traversal::parallelRathod(int maxdegree) {
    //assume u is the d simplex of d interface
    int reverted = 0;
    int rootnum = 0;

    int root = findRootIterative(maxdegree);

    while (root != -1) {

        rootnum += 1;
        
        std::cout <<"root num = " <<rootnum << "  current root = " << root << '\n';

        reverted += parallelBFS(rootnum, root);

        root = findRootIterative(maxdegree);
    }
    return reverted;
}