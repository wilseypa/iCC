#include "omp.h"

#include "QuotientAndExpand.hpp"

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getVirtualVertexIndices(const size_t maxdim, const double initeps, const int threadnumber)
{
    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);
    
    auto sorted_simplex = simplex_enumerator.getSortedVREdges(initeps);

    auto active_index_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binomial_table_, sorted_simplex, npts_);

    auto sorted_cofacet = simplex_enumerator.getSortedVRCofacets(binomial_table_, sorted_simplex, 1, initeps, threadnumber);

    std::vector<std::unordered_set<size_t>> virtual_vertex_indices;

    BipartiteGraph bi_graph(1, 1, 1);

    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        bi_graph.updateDimension(sorted_cofacet.size(), sorted_simplex.size());
        buildInterface(bi_graph, sorted_cofacet, active_index_hash_table, binomial_table_, dim);

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_simplex, sorted_cofacet);

        MaximumMorseMatching morse_matching(threadnumber);

        std::vector<std::pair<double, double>> dim_persistent_pair;

        if (dim != maxdim)
        {
            dim_persistent_pair.clear();
            auto critsimpnum = morse_matching.matchAndPersistence(matching_context, dim_persistent_pair);
            std::cout << "critical simplex number: " << critsimpnum << std::endl;
        }
        else
        {
            dim_persistent_pair.clear();
            auto critidx = morse_matching.matchAndReturnMinCriticalIndex(matching_context, dim_persistent_pair);

            if (critidx > 0) auto minfacetindex = critidx - bi_graph.unodes;

            auto minfacetweight = sorted_simplex[minfacetindex].second;

            //get gradient paths

        }
    }

}