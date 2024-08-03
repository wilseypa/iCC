#include "bi_graph.hpp"
#include "parallel_karp_sipser_init.hpp"
#include "parallel_dfs_match.hpp"
#include <iostream>

int main() {
    Bi_Graph test_rand_graph = Bi_Graph(4, 4, 2, true);
    test_rand_graph.printBiGraph();

    int unmatchedinit = parallelKarpSipserInit(&test_rand_graph, 2);

    std::cout << "unmatched after init " << unmatchedinit << std::endl;

    int unmatchedfinal = parallelDFSMatch(&test_rand_graph, 2);

    std::cout << "unmatched after dfs " << unmatchedfinal << std::endl;

    for (int i = 0; i < 5; i++) {
        std::cout << i << " match " << test_rand_graph.match[i] << '\n';
    }

    return 0;
}