#include <iostream>
#include <random>
#include <set>
#include <algorithm>
#include <execution>
#include <ranges>
#include "omp.h"

#include "bi_graph.hpp"

Bi_Graph::Bi_Graph(int leftnum, int rightnum, int maxdeg, bool random) : u(leftnum), v(rightnum), maxdegree(maxdeg)
{
    if (this->maxdegree > this->v)
    {
        this->maxdegree = this->v;
    }
    this->adj_list.resize(this->u + this->v);
    this->match.resize(this->u + this->v, -1);
    if (random)
    {
        this->randomBiGraphGen();
    }
    else
    {
        std::cout << "only random graph at this time" << std::endl;
    }
}

Bi_Graph::Bi_Graph(int leftnum, int rightnum, int maxdeg = 0) : u(leftnum), v(rightnum), maxdegree(maxdeg)
{
    if (this->maxdegree == 0 || this->maxdegree > this->v)
    {
        this->maxdegree = this->v;
    }
    this->adj_list.resize(this->u + this->v);
    this->match.resize(this->u + this->v, -1);
}

void Bi_Graph::addEdge(int u, int v)
{
    this->adj_list[this->u].push_back(this->v);
    this->adj_list[this->v].push_back(this->u);
}

void Bi_Graph::addEdge(int u, std::vector<int> &v_list)
{
    this->adj_list[this->u] = v_list;
    for (auto v : v_list)
        this->adj_list[v].push_back(u);
}

void Bi_Graph::printBiGraph()
{
    if (this->adj_list.size() > this->MAX_PRINT_SIZE)
    {
        std::cout << "graph too large" << std::endl;
        return;
    }
    for (int i = 0; i < this->u; i++)
    {
        std::cout << "left node " << i << "->";
        for (auto j : this->adj_list[i])
        {
            std::cout << j << " ";
        }
        std::cout << '\n';
    }
    for (int i = this->u; i < (this->u + this->v); i++)
    {
        std::cout << "right node " << i << "->";
        for (auto j : this->adj_list[i])
        {
            std::cout << j << " ";
        }
        std::cout << '\n';
    }
}

void Bi_Graph::randomBiGraphGen()
{
    std::random_device rand_dev;
    std::mt19937 rand_int_gen(rand_dev());
    std::set<int> u_pair;
    for (int i = 0; i < this->u; i++)
    {
        std::uniform_int_distribution<> deg_dis(1, this->maxdegree);
        std::uniform_int_distribution<> pair_dis(0, this->v - 1);
        u_pair.clear();
        int udeg = deg_dis(rand_int_gen);
        for (int j = 0; j < udeg; j++)
        {
            int vidx = pair_dis(rand_int_gen) + this->u; // global index of nodes in this->v
            if (u_pair.insert(vidx).second)
            {
                this->adj_list[i].push_back(vidx);
                this->adj_list[vidx].push_back(i);
            }
        }
    }
}

void Bi_Graph::findMatch(int uidx, std::vector<std::atomic<int>> &visit_flag, std::vector<std::atomic<int>> &node_deg)
{
    for (const auto &vidx : this->adj_list[uidx])
    {
        if (visit_flag[vidx].fetch_add(1) == 0)
        {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            this->match[uidx] = vidx;
            this->match[vidx] = uidx;
            // update degree of the neighbor of this->v
            for (const auto &index : this->adj_list[vidx])
            {
                // found new node with degree == 1
                if (node_deg[index].fetch_sub(1) == 2 && visit_flag[index].fetch_add(1) == 0)
                {
                    this->findMatch(index, visit_flag, node_deg);
                }
            }
            break;
        }
    }
}

int Bi_Graph::parallelKarpSipserInit(int threadnum)
{
    std::vector<std::atomic<int>> node_deg(this->u);
    std::vector<std::atomic<int>> visit_flag(this->u + this->v);

    omp_set_num_threads(threadnum);

    std::transform(std::execution::par, this->adj_list.begin(), this->adj_list.begin() + this->u, node_deg.begin(), [](const auto &adj)
                   { return adj.size(); });

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < this->u; i++)
    {
        if (node_deg[i].load() == 1 && visit_flag[i].fetch_add(1)== 0)
        {
            findMatch(i, visit_flag, node_deg);
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < this->u; i++)
    {
        if (node_deg[i].load() > 1 && visit_flag[i].fetch_add(1) == 0)
        {
            findMatch(i, visit_flag, node_deg);
        }
    }

    return std::count_if(std::execution::par, this->match.begin(), this->match.end(), [](int value)
                         { return value < 0; });
}

int Bi_Graph::dfsAugPath(int startnode, std::vector<int> &dfs_flag, std::vector<int> &look_ahead_flag, std::vector<int> &aug_path_tid)
{
    int topindex = -1;
    aug_path_tid[++topindex] = startnode;

    while (topindex >= 0)
    {
        int uidx = aug_path_tid[topindex];
        int endflag = 0;
        // look ahead, look for this->u's unmatched neighbor
        for (const auto &vidx : this->adj_list[uidx])
        {
            if (__sync_fetch_and_add(&(look_ahead_flag[vidx]), 1) == 0)
            {
                if (this->match[vidx] < 0)
                {
                    __sync_fetch_and_add(&(dfs_flag[vidx]), 1);
                    aug_path_tid[++topindex] = vidx;
                    return topindex + 1; // path length
                }
            }
        }
        // dfs
        for (const auto &vidx : this->adj_list[uidx])
        {
            if (__sync_fetch_and_add(&(dfs_flag[vidx]), 1) == 0)
            {
                if (!(this->match[vidx] < 0))
                {
                    aug_path_tid[++topindex] = vidx;
                    aug_path_tid[++topindex] = this->match[vidx];
                    endflag = 1;
                    break;
                }
            }
        }
        // dfs cannot augment, pop
        if (!endflag)
        {
            topindex -= 2;
        }
    }
    return topindex + 1;
}

int Bi_Graph::parallelDFSMatch(int threadnum)
{
    std::vector<int> unmatched_u_init(this->u, -1);
    std::vector<int> unmatched_u_final(this->u, 0);
    std::vector<int> dfs_flag(this->u + this->v, 0);
    std::vector<int> look_ahead_flag(this->u + this->v, 0);

    omp_set_num_threads(threadnum);

    // find unmatched left nodes from initialized this->matching
    auto view = std::ranges::views::iota(0, this->u);
    std::transform(std::execution::par, view.begin(), view.end(), unmatched_u_init.begin(), [&](const auto idx)
                   { return (this->match[idx] < 0 && this->adj_list[idx].size() > 0) ? idx : -1; });
    unmatched_u_init.erase(std::remove(std::execution::par, unmatched_u_init.begin(), unmatched_u_init.end(), -1), unmatched_u_init.end());

    int initialunmatched = unmatched_u_init.size();
    int finalunmatched = 0;

    // storing aug path one path per thread
    std::vector<std::vector<int>> aug_path(threadnum, std::vector<int>(this->u + this->v, 0));

    while (true)
    {
        // shared among threads
        finalunmatched = 0;

        std::fill(std::execution::par, dfs_flag.begin(), dfs_flag.end(), 0);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < initialunmatched; i++)
        {
            auto &aug_path_tid = aug_path[omp_get_thread_num()];

            int ustart = unmatched_u_init[i];
            int augpathlen = dfsAugPath(ustart, dfs_flag, look_ahead_flag, aug_path_tid);

            // augmentation
            for (int j = 0; j < augpathlen; j += 2)
            {
                this->match[aug_path_tid[j]] = aug_path_tid[j + 1];
                this->match[aug_path_tid[j + 1]] = aug_path_tid[j];
            }
            // store the unmatched node, need atomic op on the shared count var
            if (augpathlen <= 0)
                unmatched_u_final[__sync_fetch_and_add(&finalunmatched, 1)] = ustart;
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched))
        {
            break;
        }
        // try more aug paths, let final_unmatched be the init_unmatched
        std::swap(unmatched_u_init, unmatched_u_final);

        initialunmatched = finalunmatched;
    }

    // non-isolated unmatched
    return finalunmatched;
}