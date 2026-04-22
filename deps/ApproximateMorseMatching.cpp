#include <algorithm>
#include <iostream>

#include "omp.h"

#include "ApproximateMorseMatching.hpp"

size_t ApproximateMorseMatching::match(MatchingContext &matching_context)
{
    parallelTwoPhaseInit(matching_context.graph);
    parallelFacetDFSMatch(matching_context.graph);
    dfsCycleRemoval(matching_context.graph);
    return 0;
}

void ApproximateMorseMatching::parallelTwoPhaseInit(BipartiteGraph &graph)
{
    // convert index
    size_t umin = 0;
    size_t umax = graph.unodes;
    size_t vmin = 0 + graph.unodes;
    size_t vmax = graph.unodes + graph.vnodes;

    std::vector<int> visit_flag_v(vmax - vmin, 0);

    omp_set_num_threads(threadnum_);

// mark max facet
#pragma omp parallel for schedule(dynamic)
    for (size_t i = vmin; i < vmax; i++)
    {
        for (const auto &uidx : graph.adj_list[i]) // uidx in v adj is in ascending order
        {
            auto vit = graph.adj_list[uidx].begin(); // vidx in u adj is in descending order
            if (i == *vit)                           // i == max facet of uidx
            {
                visit_flag_v[i - vmin] += 1; // no race here
                break;
            }
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (size_t i = vmin; i < vmax; i++)
    {
        if (visit_flag_v[i - vmin] == 0)
            continue;

        int64_t min2ndcofacet = -1;
        int64_t max2ndcofacet = -1;
        int64_t min2ndfacet = i;
        int64_t max2ndfacet = 0;

        for (const auto &uidx : graph.adj_list[i]) // uidx in v adj is in ascending order
        {
            auto vit = graph.adj_list[uidx].begin(); // vidx in u adj is in descending order
            if (i == *vit)                           // i == max facet of uidx
            {
                // i is the largest facet of uidx
                // no race for uidx here
                // check the 2nd largest facet of uidx
                auto vit2nd = graph.adj_list[uidx].begin() + 1;

                if (vit2nd != graph.adj_list[uidx].end())
                {
                    // std::cout<<max2ndcofacet<<"  "<<min2ndcofacet<<"  "<<i<<'\n';

                    if (*vit2nd < min2ndfacet)
                    {
                        // std::cout<<"min 2nd cofacet pair. 2nd min facet = "<<*vit2nd<<"  "<<i<<'\n';
                        min2ndfacet = *vit2nd;
                        min2ndcofacet = uidx;
                    }
                    if (*vit2nd > max2ndfacet)
                    {
                        // std::cout<<"max 2nd cofacet pair"<<'\n';
                        max2ndfacet = *vit2nd;
                        max2ndcofacet = uidx;
                    }
                }
                else
                {
                    min2ndcofacet = uidx;
                    min2ndfacet = -1; // no other 2nd facet less than -1
                }
            }
        }

        if (min2ndcofacet >= 0)
        {
            // no race for min2ndcofacet here
            graph.match_list[min2ndcofacet] = i;
            graph.match_list[i] = min2ndcofacet;
        }

        if (max2ndcofacet >= 0 && max2ndcofacet != min2ndcofacet)
        {
            // no race for max2ndcofacet here
            if (__sync_fetch_and_add(&(visit_flag_v[max2ndfacet - vmin]), 1) == 0)
            {
                graph.match_list[max2ndcofacet] = max2ndfacet;
                graph.match_list[max2ndfacet] = max2ndcofacet;
            }
        }
    }

    return;
}

void ApproximateMorseMatching::parallelFacetDFSMatch(BipartiteGraph &graph)
{
    auto u = graph.unodes;
    auto v = graph.vnodes;

    std::vector<size_t> unmatched_v_init(v, 0);
    std::vector<size_t> unmatched_v_final(v, 0);
    std::vector<int> dfs_flag(u + v, 0);
    std::vector<int> look_ahead_flag(u, 0);

    omp_set_num_threads(threadnum_);

    size_t initialunmatched = 0;
    size_t finalunmatched = 0;
    // find unmatched left nodes from initialized matching
#pragma omp parallel
    {
        size_t ct = 0;
        std::vector<size_t> thread_buff; // overflow? limit its size?
#pragma omp for
        for (size_t i = u; i < u + v; i++)
        {
            if (graph.match_list[i] < 0 && graph.adj_list[i].size() > 0)
            {
                thread_buff.push_back(i);
                ct += 1;
            }
        }
        if (ct > 0)
        {
            size_t offset = __sync_fetch_and_add(&initialunmatched, ct);
            std::copy_n(thread_buff.begin(), ct, unmatched_v_init.begin() + offset);
        }
    }

    // storing aug path one path per thread
    std::vector<std::vector<size_t>> aug_path(threadnum_, std::vector<size_t>(u + v, 0));

    while (true)
    {
        // shared among threads
        finalunmatched = 0;

        std::fill(dfs_flag.begin(), dfs_flag.end(), 0);

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < initialunmatched; i++)
        {
            auto &aug_path_tid = aug_path[omp_get_thread_num()];

            size_t vstart = unmatched_v_init[i];
            int64_t augpathlen = facetDfsAugPath(graph, vstart, dfs_flag, look_ahead_flag, aug_path_tid);

            // augmentation
            for (int64_t j = 0; j < augpathlen; j += 2)
            {
                graph.match_list[aug_path_tid[j]] = aug_path_tid[j + 1];
                graph.match_list[aug_path_tid[j + 1]] = aug_path_tid[j];
            }
            // store the unmatched node, need atomic op on the shared count var
            if (augpathlen <= 0)
                unmatched_v_final[__sync_fetch_and_add(&finalunmatched, 1)] = vstart;
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched))
        {
            break;
        }
        // try more aug paths, let final_unmatched be the init_unmatched
        std::swap(unmatched_v_init, unmatched_v_final);

        initialunmatched = finalunmatched;
    }

    // non-isolated unmatched
    return;
}

int64_t ApproximateMorseMatching::facetDfsAugPath(BipartiteGraph &graph, const size_t startnode, std::vector<int> &dfs_flag, std::vector<int> &look_ahead_flag, std::vector<size_t> &aug_path_tid)
{
    int64_t topindex = -1;
    aug_path_tid[++topindex] = startnode;

    while (topindex >= 0)
    {
        size_t vidx = aug_path_tid[topindex];
        int endflag = 0;
        // look ahead, look for v's unmatched neighbor
        for (const auto &uidx : graph.adj_list[vidx])
        {
            if (__sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0 && graph.match_list[uidx] < 0)
            {
                // if (vidx != adj_list[uidx].back())
                {
                    __sync_fetch_and_add(&(dfs_flag[uidx]), 1);
                    aug_path_tid[++topindex] = uidx;
                    return topindex + 1; // path length
                }
            }
        }
        // dfs
        for (const auto &uidx : graph.adj_list[vidx])
        {
            if (__sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
            {
                if (graph.match_list[uidx] > 0 && graph.match_list[uidx] > vidx) // enforce direction
                {
                    aug_path_tid[++topindex] = uidx;
                    aug_path_tid[++topindex] = graph.match_list[uidx];
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

void ApproximateMorseMatching::dfsCycleRemoval(BipartiteGraph &graph)
{
    int reverted = 0;

    // 0: not visited, 1: visiting, 2: visited
    std::vector<uint8_t> state_flag(graph.unodes, 0);

    std::vector<size_t> dfs_stack;
    dfs_stack.reserve(graph.vnodes);

    int maxdegree = 5; // assume max degree is 5, can be changed
    std::vector<size_t> child_workspace;
    child_workspace.reserve(maxdegree);

    for (size_t i = 0; i < graph.unodes; i++)
    {
        if (state_flag[i] != 0)
            continue;

        dfs_stack.push_back(i);

        while (!dfs_stack.empty())
        {
            size_t top = dfs_stack.back();
            dfs_stack.pop_back();

            if (state_flag[top] == 0)
            {
                state_flag[top] = 1;
                dfs_stack.push_back(top);

                child_workspace.clear();
                // i is d-1 simp
                for (auto &i : graph.adj_list[top])
                {
                    int64_t imate = graph.match_list[i];
                    if (imate != -1 && imate != top)
                        child_workspace.push_back(imate);
                }

                for (size_t child : child_workspace)
                {

                    if (state_flag[child] == 0)
                    {
                        dfs_stack.push_back(child);
                    }
                    else if (state_flag[child] == 1)
                    {
                        // found back edge
                        int64_t temp = graph.match_list[child];

                        graph.match_list[child] = -1;
                        graph.match_list[temp] = -1;
                        reverted += 1;
                    }
                }
            }
            else if (state_flag[top] == 1)
                state_flag[top] = 2;
        }
    }

    std::cout << "Reverted cycles: " << reverted << '\n';

    return;
}
