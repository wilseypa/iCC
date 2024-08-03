#include <omp.h>
#include <execution>
#include <algorithm>

#include "parallel_karp_sipser_init.hpp"

void findMatch(Bi_Graph* bi_graph, int uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg) {
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }

    for (const auto& vidx : bi_graph->adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            bi_graph->match[uidx] = vidx;
            bi_graph->match[vidx] = uidx;
            //update degree of the neighbor of v
            for (const auto& index : bi_graph->adj_list[vidx]) {
                //found new node with degree == 1
                if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) {
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

    std::vector<int> node_deg(u, 0);
    std::vector<int> visit_flag(u + v, 0);
    int nodecount = 0;

    omp_set_num_threads(threadnum);

    std::transform(std::execution::par, bi_graph->adj_list.begin(), bi_graph->adj_list.begin() + u, node_deg.begin(), [](const auto& adj) { return adj.size(); });

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < u; i++) {
        if (node_deg[i] == 1)
            findMatch(bi_graph, i, visit_flag, node_deg);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < u; i++) {
        if (visit_flag[i] > -1 && node_deg[i] > 0) {
            findMatch(bi_graph, i, visit_flag, node_deg);
        }
    }

    return std::count_if(std::execution::par, bi_graph->match.begin(), bi_graph->match.end(), [](int value) { return value < 0; });
}


// int main() {
//     Bi_Graph test_rand_graph = Bi_Graph(7, 7, 3, true);
//     test_rand_graph.printBiGraph();

//     int unmatchedinit = parallelKarpSipserInit(&test_rand_graph, 2);

//     std::cout<<"unmatched after init "<<unmatchedinit<<std::endl;

//     return 0;
// }