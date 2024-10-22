#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <deque>
#include <queue>
#include <set>
#include <execution>
#include "omp.h"

#include "biGraphMorseMatch.hpp"

#include <iostream>

void Bi_Graph_Match::parallelKarpSipserInit() {
    std::vector<int> node_deg(u, 0);
    std::vector<int> visit_flag(u + v, 0);
    // int nodecount = 0;

    omp_set_num_threads(maxthreadnum);

    std::transform(std::execution::par, adj_list.begin(), adj_list.begin() + u, node_deg.begin(), [](const auto& u_adj) { return u_adj.size(); });

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < u; i++) {
        if (node_deg[i] == 1)
            pairDegreeOne(i, visit_flag, node_deg);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < u; i++) {
        if (visit_flag[i] == 0 && node_deg[i] > 0) {
            pairUnmatched(i, visit_flag);
        }
    }

    return;
}

void Bi_Graph_Match::pairDegreeOne(int uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg) {
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }

    std::vector<int> dfs_stack;

    for (const auto& vidx : adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            match_list[uidx] = vidx;
            match_list[vidx] = uidx;
            //update degree of the neighbor of v
            for (const auto& index : adj_list[vidx]) {
                //found new node with degree == 1
                if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) dfs_stack.push_back(index);
            }
        }
    }
    while(!dfs_stack.empty()) {
        int top = std::move(dfs_stack.back());
        dfs_stack.pop_back();

        if (__sync_fetch_and_add(&(visit_flag[top]), 1) != 0) continue;

        for (const auto& vidx : adj_list[top]) {
            if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
                match_list[top] = vidx;
                match_list[vidx] = top;
                for (const auto& index : adj_list[vidx]) {
                    if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) dfs_stack.push_back(index);
                }
            }
        }
    }
    return;
}

void Bi_Graph_Match::pairUnmatched(int uidx, std::vector<int>& visit_flag) {
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }

    for (const auto& vidx : adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            match_list[uidx] = vidx;
            match_list[vidx] = uidx;
            return;
        }
    }
    return;
}


int Bi_Graph_Match::dfsAugPath(int startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid) {
    int topindex = -1;
    aug_path_tid[++topindex] = startnode;

    while (topindex >= 0) {
        int uidx = aug_path_tid[topindex];
        int endflag = 0;
        //look ahead, look for u's unmatched neighbor
        for (const auto& vidx : adj_list[uidx]) {
            if (__sync_fetch_and_add(&(look_ahead_flag[vidx]), 1) == 0) {
                if (match_list[vidx] < 0) {
                    __sync_fetch_and_add(&(dfs_flag[vidx]), 1);
                    aug_path_tid[++topindex] = vidx;
                    return topindex + 1;    //path length
                }
            }
        }
        //dfs
        for (const auto& vidx : adj_list[uidx]) {
            if (__sync_fetch_and_add(&(dfs_flag[vidx]), 1) == 0) {
                if (!(match_list[vidx] < 0)) {
                    aug_path_tid[++topindex] = vidx;
                    aug_path_tid[++topindex] = match_list[vidx];
                    endflag = 1;
                    break;
                }
            }
        }
        //dfs cannot augment, pop
        if (!endflag) {
            topindex -= 2;
        }
    }
    return topindex + 1;
}


int Bi_Graph_Match::findRoot() {
    auto root_lambda = [&](const auto i) {return (root_flag[i] == 0);};
    auto iota{std::views::iota(0, u)};
    auto it = std::ranges::find_if(iota, root_lambda);
    if (it != iota.end()) {
        return *it;
    } else {
        return -1;
    }
}

int Bi_Graph_Match::getParent(std::vector<int>& parent_workspace, int uidx) {
    parent_workspace.clear();
    int vidx = match_list[uidx];
    //uidx is a source node has no parent (umatched)
    if (vidx == -1) return 0;
     //i is d simp
    for (auto& i: adj_list[vidx]) {
        int imate = match_list[i];
        if (imate != vidx) parent_workspace.push_back(i);
    }
    return parent_workspace.size();
}

int Bi_Graph_Match::getChild(std::vector<int>& child_workspace, int uidx) {
    child_workspace.clear();
    //i is d-1 simp
    for (auto& i: adj_list[uidx]) {
        int imate = match_list[i];
        if (imate != -1 && imate != uidx) child_workspace.push_back(imate); 
    }

    return child_workspace.size();
}

void Bi_Graph_Match::getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx) {
    //use rootnum to confine the ancestor search in current "bfs component"
    ancestor_workspace.clear();

    std::queue<int> bfs_queue;
    std::vector<int> parent_workspace;
    parent_workspace.reserve(udegree);

    //temp change uidx's rootflag to avoid cycle infi loop
    int rootflag = root_flag[uidx];
    root_flag[uidx] = -1;
    
    getParent(parent_workspace, uidx);
    for (auto& i: parent_workspace) {
        if(root_flag[i] == rootnum) bfs_queue.push(i);
    }

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

bool Bi_Graph_Match::isBackwardAcyclic(std::vector<int>& ancestor_d_simp, std::vector<int>& u_child) {
    //par
    for (auto& i: u_child) {
        if (root_flag[i] == 0) continue;
        if (std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) != ancestor_d_simp.end()) {
            return false;
        }
    }
    return true;
}

int Bi_Graph_Match::lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_d_simp, std::vector<int>& child_d_simp, int rootnum, int uidx) {
    if (child_d_simp.size() == 0) return 0;

    std::vector<int> lookahead_vec;

    //par ch workspace do not overlap. can be reused
    std::vector<int> parent_workspace;
    parent_workspace.reserve(udegree);

    std::vector<int> child_workspace;
    child_workspace.reserve(udegree);

    for (auto& i: child_d_simp) {
        int parentnum = getParent(parent_workspace, i);
        //if uidx is the only parent of i (in the whole graph). use look ahead shortcut
        if (parentnum == 1) {
            lookahead_vec.push_back(i);
        } else if (root_flag[i] == 0 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
            //if i has more than 1 par and has not been searched yet (or not in current bfs component). push to global bfs queue
            graph_bfs_queue.push_back(i);
        }
    }

    //start look ahead op on children of uidx
    ancestor_d_simp.push_back(uidx);
    int reverted = 0;
    bool flag;    //working var
    //the nodes appeared in lookahead stack at the same time are independent 
    while (!lookahead_vec.empty()) {

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
            
            int temp = match_list[top];
            match_list[top] = -1;
            match_list[temp] = -1;
            reverted += 1;

        } else {
        //check if the child of top can be pushed to look ahead stack
            for (auto& i: child_workspace) {
                if (getParent(parent_workspace, i) == 1) {
                    //if top is the only parent of i. keep look ahead op
                    lookahead_vec.push_back(i);
                    //raise flag for adding top to ancestor
                    flag = true;
                }else if (root_flag[i] == 0 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
                    graph_bfs_queue.push_back(i);
                } 
            }
        }
        //if any of top's child is pushed to look ahead stack, update ancestor
        if (flag) ancestor_d_simp.push_back(top);
        //not all nodes in the stack have ances-des relationship with top
        //but this does not bring extra cycle rm. simplfied implementation
    }
    return reverted;
 }

 int Bi_Graph_Match::serialBFS(std::vector<int>& ancestor_d_simp, int rootnum, int root) {
    int reverted = 0;

    std::deque<int> bfs_queue;
    bfs_queue.push_back(root);

    std::vector<int> front_child;
    front_child.reserve(udegree);

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
            int temp = match_list[front];
            match_list[front] = -1;
            match_list[temp] = -1;
            reverted += 1;
        } else {
            //look ahead and push descendant to bfs queue
            reverted += lookAheadDFS(bfs_queue, ancestor_d_simp, front_child, rootnum, front);
        }
    }
    
    return reverted;
}


void Bi_Graph_Match::parallelDFSMatch() {

    std::vector<int> unmatched_u_init(u, -1);
    std::vector<int> unmatched_u_final(u, 0);
    std::vector<int> dfs_flag(u + v, 0);
    std::vector<int> look_ahead_flag(u + v, 0);

    omp_set_num_threads(maxthreadnum);

    int initialunmatched = 0;
    int finalunmatched = 0;
    //find unmatched left nodes from initialized matching
#pragma omp parallel
    {
        int ct = 0;
        std::vector<int> thread_buff;    //overflow? limit its size?
#pragma omp for
        for (int i = 0; i < u; i++) {
            if (match_list[i] < 0 && adj_list[i].size() > 0) {
                thread_buff.push_back(i);
                ct += 1;
            }
        }
        if (ct > 0) {
            int offset = __sync_fetch_and_add(&initialunmatched, ct);
            std::copy_n(thread_buff.begin(), ct, unmatched_u_init.begin() + offset);
        }
    }


    //storing aug path one path per thread
    std::vector<std::vector<int>> aug_path(maxthreadnum, std::vector<int>(u + v, 0));

    while (true) {
        //shared among threads
        finalunmatched = 0;

        std::fill(std::execution::par, dfs_flag.begin(), dfs_flag.end(), 0);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < initialunmatched; i++) {
            auto& aug_path_tid = aug_path[omp_get_thread_num()];

            int ustart = unmatched_u_init[i];
            int augpathlen = dfsAugPath(ustart, dfs_flag, look_ahead_flag, aug_path_tid);

            //augmentation
            for (int j = 0; j < augpathlen; j += 2) {
                match_list[aug_path_tid[j]] = aug_path_tid[j + 1];
                match_list[aug_path_tid[j + 1]] = aug_path_tid[j];
            }
            //store the unmatched node, need atomic op on the shared count var
            if (augpathlen <= 0)
                unmatched_u_final[__sync_fetch_and_add(&finalunmatched, 1)] = ustart;
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched)) {
            break;
        }
        //try more aug paths, let final_unmatched be the init_unmatched
        std::swap(unmatched_u_init, unmatched_u_final);

        initialunmatched = finalunmatched;
    }

    //non-isolated unmatched
    return;   
}

int Bi_Graph_Match::serialCycleRemoval() {
    //assume u is the d simplex of d interface
    root_flag.resize(u, 0);

    int reverted = 0;

    int root = findRoot();
    int rootnum = 0;

    //ancestor workspace. avoid multiple reallocation
    std::vector<int> ancestor_workspace;
    ancestor_workspace.reserve(u);

    while (root != -1) {

        // std::cout << "current root = " << root << '\n';

        rootnum += 1;

        reverted += serialBFS(ancestor_workspace, rootnum, root);

        root = findRoot();
    }
    return reverted;
}

// std::vector<int> Bi_Graph_Match::getTopDimCriticalIndex() 
// {
//     std::vector<int> critical_index;
//     for(int i = 0; i < u; i++) 
//     {
//         if (match_list[i] < 0) critical_index.push_back(i);
//     }
//     return critical_index;
// }


void Bi_Graph_Match::addEdge(int leftnode, int rightnode) {
    adj_list[leftnode].push_back(u + rightnode);
    adj_list[u + rightnode].push_back(leftnode);
    return;
}

Bi_Graph_Match::Bi_Graph_Match(int leftnum, int rightnum, int leftdim, int threadnum) : u(leftnum), v(rightnum), udegree(leftdim), maxthreadnum(threadnum) {
    adj_list.resize(u + v);
    match_list.resize(u + v, -1);
}

void Bi_Graph_Match::updateDimension(int newleftnum, int newrightnum) {
    u = newleftnum;
    v = newrightnum;
    udegree += 1;
    adj_list.clear();
    adj_list.resize(u + v);
    match_list.clear();
    match_list.resize(u + v, -1);
}

bool Bi_Graph_Match::isFacet(std::vector<int>& cofacet, std::vector<int>& facet)
{
    for(auto i: facet)
    {
        if (std::find(cofacet.begin(), cofacet.end(), i) == cofacet.end()) return false;
    }
    return true;
}

void Bi_Graph_Match::buildInterface(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& active_index) {    
    //openmp??????

    for (int i = 0; i < cofacet_bin.size(); i++) {
        for (auto j: active_index) {
            // if (std::includes(cofacet_bin[i].begin(), cofacet_bin[i].end(), simplex_bin[j].begin(), simplex_bin[j].end())) 
            if (isFacet(cofacet_bin[i], simplex_bin[j]))
            {   
                addEdge(i, j);
            }
            
        }
    }

    // for(int i = 0; i < 20; i++) std::cout<<adj_list[i].size()<<"  ";
    // std::cout<<'\n';
    // for(int i = 0; i < 20; i++) std::cout<<adj_list[u+i].size()<<"  ";
    // std::cout<<'\n';

}

std::vector<int> Bi_Graph_Match::getActiveIndex() {
    std::vector<int> active_index;
    //need to reserve the space?
    for (int i = 0; i < u; i++) {
        if (match_list[i] < 0) active_index.push_back(i);
    }
    return active_index;
}

std::vector<int> Bi_Graph_Match::getCriticalIndex(std::vector<int>& dim_active_index) 
{
    std::vector<int> critical_index;
    //need to reserve the space?
    for (int i = u; i < (u + v); i++) {
        //for kth simplex, its index in graph i == k + u
        if (match_list[i] < 0)
        {
            if (std::find(dim_active_index.begin(), dim_active_index.end(), i - u) != dim_active_index.end()) critical_index.push_back(i - u);
        }            
    }
    return critical_index;
}

//to do
//if epsmin != 0, pass empty vec as d-1 simp
//bi_graph_match.build_interface(d-simp, d-1-simp, dimension)

//outline
//get sorted edge
//get first edge bin
//for dim = 1 to n
//build interface 
//save d-1 cirtical
//build cc.cofacet. take d-1 simp's spot
//swap d d-1 simp