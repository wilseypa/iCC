#include <iostream>
#include <vector>
#include <omp.h>
#include <algorithm>

#include "bi_graph.h"
#include "parallel_karp_sipser_init.h"

int dfsAugPath(Bi_Graph* bi_graph, int startnode, int* dfs_flag, int* look_ahead_flag, int* aug_path_tid) {
    int topindex = -1;
    aug_path_tid[++topindex] = startnode;

    while(topindex >= 0) {
        int uidx = aug_path_tid[topindex];
        int endflag = 0;
        //look ahead, look for u's unmatched neighbor
        for (const auto& vidx: bi_graph->adj_list[uidx]) {
            if (__sync_fetch_and_add(&(look_ahead_flag[vidx]), 1) == 0) {
                if (bi_graph->match[vidx] < 0) {
                    __sync_fetch_and_add(&(dfs_flag[vidx]), 1);
                    aug_path_tid[++topindex] = vidx;
                    return topindex + 1;    //path length
                }
            }
        }
        //dfs
        for (const auto& vidx: bi_graph->adj_list[uidx]) {
            if (__sync_fetch_and_add(&(dfs_flag[vidx]), 1) == 0) {
                if (!(bi_graph->match[vidx] < 0)) {
                    aug_path_tid[++topindex] = vidx;
                    aug_path_tid[++topindex] = bi_graph->match[vidx];
                    break;
                }
            }
            //searched all u's neighbor
            endflag = 1;
        }
        //dfs cannot augment, pop
        if (endflag) {
            topindex -= 2;
        }
    }
    return topindex + 1;
}

int parallelDFSMatch(Bi_Graph *bi_graph, int threadnum) {
    int u = bi_graph->u;
    int v = bi_graph->v;

    int* unmatched_u_init = new int[u];
    int* unmatched_u_final = new int[u];
    int* dfs_flag = new int[u + v];
    int* look_ahead_flag = new int[u + v];

    int initialunmatched = 0;
    int finalunmatched = 0;
    int iterationcount = 0;

    omp_set_num_threads(threadnum);

    //find unmatched left nodes from initialized matching
#pragma omp parallel
    {
        int ct = 0;
        std::vector<int> thread_buff;    //overflow? limit its size?
#pragma omp for
        for (int i = 0; i < u; i++) {
            if (bi_graph->match[i] < 0 && bi_graph->adj_list[i].size() > 0) {
                thread_buff.push_back(i);
                ct += 1;
            }
        }
        if (ct > 0) {
            int offset = __sync_fetch_and_add(&initialunmatched, ct);
            std::copy_n(thread_buff.begin(), ct, unmatched_u_init + offset);
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < (u + v); i++) {
        look_ahead_flag[i] = 0;
    }

    //storing aug path
    //one path per thread
    int** aug_path = new int* [threadnum];
    for (int i = 0; i < threadnum; i++) {
        aug_path[i] = new int [u + v];
    }

    while(true) {
        //shared among threads
        iterationcount += 1;
        finalunmatched = 0;

#pragma omp parallel for schedule(static)
        for (int i = 0; i < (u + v); i++) {
            dfs_flag[i] = 0;
            //no need to reset look ahead flag each time
        }

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < initialunmatched; i++) {
            auto tid = omp_get_thread_num();
            int* aug_path_tid = aug_path[tid];

            int ustart = unmatched_u_init[i];
            int augpathlen = -1;

            augpathlen = dfsAugPath(bi_graph, ustart, dfs_flag, look_ahead_flag, aug_path_tid);

            //augmentation
            if (augpathlen > 0) {
                for (int j = 0; j < augpathlen; j += 2) {
                    bi_graph->match[aug_path_tid[j]] = aug_path_tid[j + 1];
                    bi_graph->match[aug_path_tid[j + 1]] = aug_path_tid[j];
                }
            }
            //store the unmatched node, need atomic op on the shared count var
            else {
                unmatched_u_final[__sync_fetch_and_add(&finalunmatched, 1)] = ustart;
            }
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched)) {
            break;
        }
        //try more aug paths, let final_unmatched be the init_unmatched
        int* tempptr = unmatched_u_final;
        unmatched_u_final = unmatched_u_init;
        unmatched_u_init = tempptr;

        initialunmatched = finalunmatched;
    }

    //deallocate
    delete[] unmatched_u_init;
    delete[] unmatched_u_final;
    delete[] dfs_flag;
    delete[] look_ahead_flag;
    for (int i = 0; i < threadnum; i++){
        delete[] aug_path[i];
    }
    delete[] aug_path;

    //non-isolated unmatched
    return finalunmatched;
}

int main() {
    Bi_Graph test_rand_graph = Bi_Graph(4, 4, 2, true);
    test_rand_graph.printBiGraph();

    int unmatchedinit = parallelKarpSipserInit(&test_rand_graph, 2);

    std::cout<<"unmatched after init "<<unmatchedinit<<std::endl;

    int unmatchedfinal = parallelDFSMatch(&test_rand_graph, 2);

    std::cout<<"unmatched after dfs "<<unmatchedfinal<<std::endl;

    for (int i = 0; i < 5; i++) {
        std::cout<<i<<" match "<<test_rand_graph.match[i]<<'\n';
    }

    return 0;
}