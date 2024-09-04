#include <iostream>

#include "bi_graph.hpp"
#include "parallel_karp_sipser_init.hpp"
#include "parallel_dfs_match.hpp"
#include "cycle_rm.hpp"

int main() {
    Bi_Graph rand_graph = Bi_Graph(100, 80, 4, true);
    // test_rand_graph.printBiGraph();

    int unmatchedinit = parallelKarpSipserInit(&rand_graph, 4);

    std::cout<<"unmatched after init "<<unmatchedinit<<'\n';

    int unmatchedfinal = parallelDFSMatch(&rand_graph, 4);

    std::cout << "unmatched after dfs " << unmatchedfinal << '\n';


    Bi_Graph_Traversal bi_graph_trav = Bi_Graph_Traversal(&rand_graph);

    // int revertedmatch = bi_graph_trav.serialCycleRemoval(4);
    int revertedmatch = bi_graph_trav.parallelRathod(4);


    std::cout << "reverted match " << revertedmatch << std::endl;

    return 0;
}
