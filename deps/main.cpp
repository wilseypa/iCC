#include <iostream>
#include <chrono>

#include "bi_graph.hpp"
#include "cycle_rm.hpp"

int main()
{
    Bi_Graph rand_graph = Bi_Graph(1000, 800, 4, true);
    // test_rand_graph.printBiGraph();

    int unmatchedinit = rand_graph.parallelKarpSipserInit(4);

    std::cout << "unmatched after init " << unmatchedinit << '\n';

    int unmatchedfinal = rand_graph.parallelDFSMatch(4);

    std::cout << "unmatched after dfs " << unmatchedfinal << '\n';

    // for (int i = 0; i < rand_graph.u; i++) {
    //     int imate = rand_graph.match[i];
    //     if (imate != -1 && rand_graph.match[imate] != i) {
    //         std::cout<<"matching error !!!  i mate = "<<rand_graph.match[i]<<"  imate match = "<<rand_graph.match[imate]<<'\n';
    //     }
    // }

    Bi_Graph_Traversal bi_graph_trav0 = Bi_Graph_Traversal(&rand_graph);
    // Bi_Graph rand_graph_cp(rand_graph);
    // Bi_Graph_Traversal bi_graph_trav1 = Bi_Graph_Traversal(&rand_graph_cp);

    auto st0 = std::chrono::high_resolution_clock::now();
    int revertedmatch0 = bi_graph_trav0.serialCycleRemoval(4);
    auto st1 = std::chrono::high_resolution_clock::now();
    auto st_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);

    // auto pt0 = std::chrono::high_resolution_clock::now();
    // int revertedmatch1 = bi_graph_trav1.parallelRathod(4);
    // auto pt1 = std::chrono::high_resolution_clock::now();
    // auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pt1 - pt0);

    std::cout << "serial reverted = " << revertedmatch0 << '\n';
    // std::cout << "para reverted = " << revertedmatch1 << '\n';

    // auto diff = st_ms.count() - pt_ms.count();

    // std::cout << "time difference (ser - para) in ms = " << diff << std::endl;

    return 0;
}