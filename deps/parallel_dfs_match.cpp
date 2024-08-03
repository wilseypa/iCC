#include <omp.h>
#include <algorithm>
#include <ranges>
#include <execution>

#include "parallel_dfs_match.hpp"

int dfsAugPath(Bi_Graph* bi_graph, int startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid) {
    int topindex = -1;
    aug_path_tid[++topindex] = startnode;

    while (topindex >= 0) {
        int uidx = aug_path_tid[topindex];
        //look ahead, look for u's unmatched neighbor
        for (const auto& vidx : bi_graph->adj_list[uidx]) {
            if (__sync_fetch_and_add(&(look_ahead_flag[vidx]), 1) == 0) {
                if (bi_graph->match[vidx] < 0) {
                    __sync_fetch_and_add(&(dfs_flag[vidx]), 1);
                    aug_path_tid[++topindex] = vidx;
                    return topindex + 1;    //path length
                }
            }
        }
        //dfs
        for (const auto& vidx : bi_graph->adj_list[uidx]) {
            if (__sync_fetch_and_add(&(dfs_flag[vidx]), 1) == 0) {
                if (!(bi_graph->match[vidx] < 0)) {
                    aug_path_tid[++topindex] = vidx;
                    aug_path_tid[++topindex] = bi_graph->match[vidx];
                    break;
                }
            }
        }
        //dfs cannot augment, pop
        if (!bi_graph->adj_list[uidx].empty()) {
            topindex -= 2;
        }
    }
    return topindex + 1;
}

int parallelDFSMatch(Bi_Graph* bi_graph, int threadnum) {
    int u = bi_graph->u;
    int v = bi_graph->v;

    std::vector<int> unmatched_u_init(u, -1);
    std::vector<int> unmatched_u_final(u, 0);
    std::vector<int> dfs_flag(u + v, 0);
    std::vector<int> look_ahead_flag(u + v, 0);

    omp_set_num_threads(threadnum);

    //find unmatched left nodes from initialized matching
    auto view = std::ranges::views::iota(0, u);
    std::transform(std::execution::par, view.begin(), view.end(), unmatched_u_init.begin(), [&](const auto idx) {return (bi_graph->match[idx] < 0 && bi_graph->adj_list[idx].size() > 0) ? idx : -1;});
    unmatched_u_init.erase(std::remove(std::execution::par, unmatched_u_init.begin(), unmatched_u_init.end(), -1), unmatched_u_init.end());

    int initialunmatched = unmatched_u_init.size();
    int finalunmatched = 0;

    //storing aug path one path per thread
    std::vector<std::vector<int>> aug_path(threadnum, std::vector<int>(u + v, 0));

    while (true) {
        //shared among threads
        finalunmatched = 0;

        std::fill(std::execution::par, dfs_flag.begin(), dfs_flag.end(), 0);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < initialunmatched; i++) {
            auto& aug_path_tid = aug_path[omp_get_thread_num()];

            int ustart = unmatched_u_init[i];
            int augpathlen = dfsAugPath(bi_graph, ustart, dfs_flag, look_ahead_flag, aug_path_tid);

            //augmentation
            for (int j = 0; j < augpathlen; j += 2) {
                bi_graph->match[aug_path_tid[j]] = aug_path_tid[j + 1];
                bi_graph->match[aug_path_tid[j + 1]] = aug_path_tid[j];
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
    return finalunmatched;
}