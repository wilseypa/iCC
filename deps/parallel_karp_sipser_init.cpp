#include <iostream>
#include <vector>
#include <omp.h>
#include <algorithm>

#include "bi_graph.h"
#include "parallel_karp_sipser_init.h"

void findMatch(Bi_Graph* bi_graph, int uidx, int* visit_flag, int* node_deg) {
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }

    for (const auto& vidx: bi_graph->adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            bi_graph->match[uidx] = vidx;
            bi_graph->match[vidx] = uidx;
            //update degree of the neighbor of v
            for (const auto& index: bi_graph->adj_list[vidx]) {
                //found new node with degree == 1
                if (__sync_fetch_and_add(&(node_deg[index]), -1) == 2) {
                    findMatch(bi_graph, index, visit_flag, node_deg);
                }
            }
            break;
        }
    }
}

int parallelKarpSipserInit(Bi_Graph* bi_graph, int threadnum) {
    int u = bi_graph->u;
    int v = bi_graph->v;

    int* node_deg = new int[u];
    int* deg_one_node = new int[u];
    int* visit_flag = new int[u + v];
    int nodecount = 0;
    int leftunmatched = 0;

    std::fill_n(visit_flag, u + v, 0);    //omp para as well?

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < u; i++) {
        node_deg[i] = bi_graph->adj_list[i].size();
        // std::cout<<"node deg "<<node_deg[i]<<'\n';
        if (node_deg[i] == 1) {
            deg_one_node[__sync_fetch_and_add(&nodecount, 1)] = i;
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nodecount; i++) {
        findMatch(bi_graph, deg_one_node[i], visit_flag, node_deg);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < u; i++) {
        if (visit_flag[i] > -1 && node_deg[i] > 0) {
            findMatch(bi_graph, i, visit_flag, node_deg);
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < u; i++) {
        if (bi_graph->match[i] < 0) {
            __sync_fetch_and_add(&leftunmatched, 1);
        }
    }

delete[] node_deg;
delete[] deg_one_node;
delete[] visit_flag;

return leftunmatched;
}


// int main() {
//     Bi_Graph test_rand_graph = Bi_Graph(7, 7, 3, true);
//     test_rand_graph.printBiGraph();

//     int unmatchedinit = parallelKarpSipserInit(&test_rand_graph, 2);

//     std::cout<<"unmatched after init "<<unmatchedinit<<std::endl;

//     return 0;
// }