#include <unordered_map>

#include "omp.h"

#include "SimplexUtility.hpp"
#include "SimplexEnumerator.hpp"
#include "MaximumMorseMatching.hpp"
#include "QuotientAndExpand.hpp"

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getVirtualVertexIndices(const size_t maxdim, const double initeps, const int threadnumber)
{
    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);
    
    auto sorted_simplex = simplex_enumerator.getSortedVREdges(initeps);

    auto active_index_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binomial_table_, sorted_simplex, npts_);

    auto sorted_cofacet = simplex_enumerator.getSortedVRCofacets(binomial_table_, sorted_simplex, 1, initeps, threadnumber);

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

            if (minfacetweight < 0.0) minfacetweight = initeps;

            //get gradient paths
            auto gradient_paths = extractGradientPahts(matching_context, minfacetweight);
            //get gradient path vertex indices
            auto virtual_vertex_indices = getGradientPathVertexSets(matching_context, gradient_paths, dim);
        }

        if (dim != maxdim)
        {
            active_index_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binomial_table_, sorted_cofacet, npts_);

            sorted_simplex = simplex_enumerator.getSortedVRCofacets(binomial_table_, sorted_cofacet, dim, initeps, threadnumber);
            std::swap(sorted_simplex, sorted_cofacet);
        }

    }

    return virtual_vertex_indices;
}


template <typename DistMatType>
std::vector< std::vector<size_t> > QuotientAndExpand<DistMatType>::extractGradientPahts(const MatchingContext& matching_context, const double minfacetweight)
{
    auto bi_graph = matching_context.graph;
    auto sorted_simplex = matching_context.sorted_facets;
    auto sorted_cofacet = matching_context.sorted_cofacets;

    const size_t u = bi_graph.unodes;
    const size_t v = bi_graph.vnodes;

    auto facetcutoff = std::lower_bound(sorted_simplex.begin(), sorted_simplex.end(), minfacetweight,
        [](const std::pair<int64_t, double>& pair, double weight) { return pair.second < weight; });

    size_t maxfacetindex = std::distance(sorted_simplex.begin(), facetcutoff);

    auto cofacetcutoff = std::lower_bound(sorted_cofacet.begin(), sorted_cofacet.end(), minfacetweight,
        [](const std::pair<int64_t, double>& pair, double weight) { return pair.second < weight; });

    size_t maxcofacetindex = std::distance(sorted_cofacet.begin(), cofacetcutoff);

    GradientPathUnionFind union_find(u);

    std::vector<bool> is_in_gradient_path(u, false);

    for (size_t facetidx = 0; facetidx < maxfacetindex; ++facetidx)
    {
        size_t cofacetidx = bi_graph.match_list[facetidx + u];

        if (cofacetidx < 0 || cofacetidx >= maxcofacetindex) continue;

        is_in_gradient_path[cofacetidx] = true;

        for (auto vidx : bi_graph.adj_list[cofacetidx])
        {
            size_t nextfacetidx = vidx - u;

            if (nextfacetidx == facetidx || nextfacetidx >= maxfacetindex) continue;

            size_t nextcofacetidx = bi_graph.match_list[vidx];

            if (nextcofacetidx < 0 || nextcofacetidx >= maxcofacetindex) continue;

            is_in_gradient_path[nextcofacetidx] = true;
            union_find.unionSets(cofacetidx, nextcofacetidx);
        }
    }

    std::unordered_map<size_t, std::vector<size_t>> component_map;

    for (size_t cofacetidx = 0; cofacetidx < maxcofacetindex; ++cofacetidx)
    {
        if (is_in_gradient_path[cofacetidx])
        {
            size_t root = union_find.unionFind(cofacetidx);
            component_map[root].push_back(cofacetidx);
        }
    }

    std::vector<std::vector<size_t>> gradient_paths;
    gradient_paths.reserve(component_map.size());

    for (const auto& [root, path] : component_map)    //structured binding requires C++17
    {
        if (!path.empty())
        {
            gradient_paths.push_back(std::move(path));
        }
    }

    return gradient_paths;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getGradientPathVertexSets(const MatchingContext& matching_context, const std::vector< std::vector<size_t> >& gradient_paths, const size_t dim)
{
    std::vector<std::unordered_set<size_t>> gradient_path_vertex_sets;
    gradient_path_vertex_sets.reserve(gradient_paths.size());

    auto npt = matching_context.graph.unodes + matching_context.graph.vnodes;

    for (const auto& grad_path: gradient_paths)
    {
        std::unordered_set<size_t> vertex_set;

        for (auto cofacetidx : grad_path)
        {
            auto bindex = matching_context.sorted_cofacets[cofacetidx].first;
            auto simplex_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, bindex, npt, dim);
            vertex_set.insert(simplex_vertices.begin(), simplex_vertices.end());
        }
        gradient_path_vertex_sets.push_back(std::move(vertex_set));
    }

    return gradient_path_vertex_sets;
}