#pragma once

#include "BipartiteGraph.hpp"

#include "robin_hood.h"

#include <vector>
#include <cstdint>

struct MatchingContext
{
    BipartiteGraph& graph;

    // References to data needed by the algorithms
    std::vector<std::vector<int64_t>>& binomial_table;
    std::vector<std::pair<int64_t, double>>& sorted_facets;
    std::vector<std::pair<int64_t, double>>& sorted_cofacets;

    robin_hood::unordered_map<int64_t, size_t> active_index_hash_table;

    MatchingContext (BipartiteGraph& bi_graph,
                     std::vector<std::vector<int64_t>>& binomial_table,
                     std::vector<std::pair<int64_t, double>>& sorted_facets,
                     std::vector<std::pair<int64_t, double>>& sorted_cofacets,
                     robin_hood::unordered_map<int64_t, size_t>& active_index_hash_table)
        : graph(bi_graph),
          binomial_table(binomial_table),
          sorted_facets(sorted_facets),
          sorted_cofacets(sorted_cofacets),
          active_index_hash_table(active_index_hash_table) {}

};