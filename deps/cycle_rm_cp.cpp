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


int Bi_Graph_Traversal::findRoot() {
    int u = graphptr->u;
    int root = -1;

    for(int i = 0; i < u; i++) {
        if (root_flag[i] == 0) {
            root = i;
            break;
        }
    }
    return root;
}


//return the number of parent
int Bi_Graph_Traversal::getParent(std::vector<int>& parent_workspace, int uidx) {
    parent_workspace.clear();
    int vidx = graphptr->match[uidx];
    //uidx is a source node has no parent (umatched)
    if (vidx == -1) return 0;
     //i is d simp
    for (auto& i: graphptr->adj_list[vidx]) {
        int imate = graphptr->match[i];
        if (imate != vidx) parent_workspace.push_back(i);
    }
    return parent_workspace.size();
}

//no need to check visit flag. if found child in searched part. it might be merging or cycle
//in later stages, need to check visit flag when pushing to bfs queue
int Bi_Graph_Traversal::getChild(std::vector<int>& child_workspace, int uidx) {
    child_workspace.clear();
    //i is d-1 simp
    for (auto& i: graphptr->adj_list[uidx]) {
        int imate = graphptr->match[i];
        if (imate != -1 && imate != uidx) child_workspace.push_back(imate); 
    }

    return child_workspace.size();
}

//ancestor in bfs searched part
void Bi_Graph_Traversal::getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx) {
    //use rootnum to confine the ancestor search in current "bfs component"
    ancestor_workspace.clear();

    std::queue<int> bfs_queue;
    std::vector<int> parent_workspace;
    parent_workspace.reserve(graphptr->maxdegree);

    //temp change uidx's rootflag to avoid cycle infi loop
    int rootflag = root_flag[uidx];
    root_flag[uidx] = -1;
    
    getParent(parent_workspace, uidx);
    for (auto& i: parent_workspace) {
        if(root_flag[i] == rootnum) bfs_queue.push(i);
    }

    // std::cout<<"get ancestor started \n";
    //cycle can be found during this process but need ancestors for look up
    while (!bfs_queue.empty()) {
        int front = std::move(bfs_queue.front());
        bfs_queue.pop();
        ancestor_workspace.push_back(front);
        
        //parent of queue front
        getParent(parent_workspace, front);
        for(auto& i: parent_workspace) {
            if (root_flag[i] == rootnum && std::find(ancestor_workspace.begin(), ancestor_workspace.end(), i) == ancestor_workspace.end()) {
                bfs_queue.push(i);
            }
        }
    }

    //recover rootflag
    root_flag[uidx] = rootflag;

    return;
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

    std::vector<int> parent_workspace;
    parent_workspace.reserve(graphptr->maxdegree);

    std::vector<int> child_workspace;
    child_workspace.reserve(graphptr->maxdegree);

    for (auto& i: child_d_simp) {
        int parentnum = getParent(parent_workspace, i);
        //if uidx is the only parent of i (in the whole graph). use look ahead shortcut
        if (parentnum == 1) {
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

        int top = std::move(lookahead_vec.back());
        lookahead_vec.pop_back();
        root_flag[top] = rootnum;

        //child_workspace contains the children of top
        int childnum = getChild(child_workspace, top);
        if (childnum == 0) continue;

        //child in unsearched part of current bfs component
        auto iter = std::remove_if(child_workspace.begin(), child_workspace.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});
        child_workspace.erase(iter, child_workspace.end());

        //if top is not acyclic. the lookahead of top is over
        //push its unsearched children or reencountered children to bfs queue
        if (!isBackwardAcyclic(ancestor_d_simp, child_workspace)) {
            
            int temp = graphptr->match[top];
            graphptr->match[top] = -1;
            graphptr->match[temp] = -1;
            reverted += 1;

        } else {
        //check if the child of top can be pushed to look ahead stack
            for (auto& i: child_workspace) {
                if (root_flag[i] == 0 && getParent(parent_workspace, i) == 1) {
                    //if top is the only parent of i. keep look ahead op
                    lookahead_vec.push_back(i);
                    //raise flag for adding top to ancestor
                    flag = true;
                }else if (std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
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


int Bi_Graph_Traversal::serialBFS(std::vector<int>& ancestor_d_simp, int rootnum, int root) {
    int reverted = 0;

    std::deque<int> bfs_queue;
    bfs_queue.push_back(root);

    // std::cout<<"started at root = "<<root<<'\n';
    std::vector<int> front_child;
    front_child.reserve(graphptr->maxdegree);

    while (!bfs_queue.empty()) {
        int front = std::move(bfs_queue.front());
        bfs_queue.pop_front();
        root_flag[front] = rootnum;
        
        int childnum = getChild(front_child, front);
        if (childnum == 0) continue;

        //child in unsearched part of current bfs component
        auto iter = std::remove_if(front_child.begin(), front_child.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});
        front_child.erase(iter, front_child.end());

        getAncestor(ancestor_d_simp, rootnum, front);

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

    int root = findRoot();
    int rootnum = 0;

    //ancestor workspace. avoid multiple reallocation
    std::vector<int> ancestor_workspace;
    ancestor_workspace.reserve(graphptr->u);

    while (root != -1) {

        // std::cout << "current root = " << root << '\n';

        rootnum += 1;

        reverted += serialBFS(ancestor_workspace, rootnum, root);

        root = findRoot();
    }
    return reverted;
}


int Bi_Graph_Traversal::threadBackwardBFS(int rootnum, int uparent, int uidx) {

    int reverted = 0;

    std::vector<int> child_workspace;
    child_workspace.reserve(graphptr->maxdegree);

    getChild(child_workspace, uidx);
    auto iter = std::remove_if(child_workspace.begin(), child_workspace.end(), [&](const auto& i) {return (root_flag[i] != rootnum);});
    child_workspace.erase(iter, child_workspace.end());

    if (child_workspace.size() == 0) return reverted;

    //backward bfs
    std::queue<int> bfs_queue;
    std::vector<int> desendant_workspace;
    desendant_workspace.reserve(graphptr->u);
    
    for (auto& i: child_workspace) {
        bfs_queue.push(i);
    }

    while (!bfs_queue.empty()) {
        int front = std::move(bfs_queue.front());
        desendant_workspace.push_back(front);
        bfs_queue.pop();

        getChild(child_workspace, front);
        for(auto& i: child_workspace) {
            //found cycle
            if (i == uparent) {
                int temp = graphptr->match[uidx];
                graphptr->match[uidx] = -1;
                graphptr->match[temp] = -1;
                reverted += 1;
                return reverted;
            }
            if (std::find(desendant_workspace.begin(), desendant_workspace.end(), i) == desendant_workspace.end()) {
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

    std::vector<int> front_child;
    front_child.reserve(graphptr->maxdegree);

    auto push_lambda = [&](const auto i) {return (graphptr->match[i] != -1 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end());};

    while (!graph_bfs_queue.empty()) {
        int front = std::move(graph_bfs_queue.front());
        graph_bfs_queue.pop_front();
        root_flag[front] = rootnum;

        // std::cout<<"root = " <<root<<"  front = "<<front<<'\n';
        
        //get "leading up-edges"
        int childnum = getChild(front_child, front);
        if (childnum == 0) continue;
        //keep the unvisited or reencountered
        //remove_if only logically remove elements, read the docs!!!

        auto iter = std::remove_if(front_child.begin(), front_child.end(), [&](const auto& i) {return (root_flag[i] != 0 && root_flag[i] != rootnum);});

        front_child.erase(iter, front_child.end());

        //dynamic distribute work
#pragma omp parallel for schedule(dynamic) reduction(+ : reverted)
        for (int i = 0; i < front_child.size(); i++)
        {
            reverted += threadBackwardBFS(rootnum, front, front_child[i]);
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

    int root = findRoot();

    while (root != -1) {

        rootnum += 1;
        
        // std::cout <<"root num = " <<rootnum << "  current root = " << root << '\n';

        reverted += parallelBFS(rootnum, root);

        root = findRoot();
    }
    return reverted;
}