#pragma once

#include "BipartiteGraph.hpp"

#include "robin_hood.h"

#include <vector>
#include <cstdint>

namespace    //for explicit/legacy
{
  inline const robin_hood::unordered_map<int64_t, size_t>& emptyMap() 
  {
    static const robin_hood::unordered_map<int64_t, size_t> empty_map;
    return empty_map;
  }
}

struct MatchingContext
{
    //default values for explicit/legacy implementation
    const size_t npts = 0;    //total number of vertices (original + virtual)
    const size_t dim = 0;    //current dim of interface (cofacet dim)

    BipartiteGraph& graph;

    // References to data needed by the algorithms
    std::vector<std::vector<int64_t>>& binomial_table;
    std::vector<std::pair<int64_t, double>>& sorted_facets;
    std::vector<std::pair<int64_t, double>>& sorted_cofacets;

    // Hash tables required for implicit graph lookups
    const robin_hood::unordered_map<int64_t, size_t>& facet_bindex_to_list_index;
    const robin_hood::unordered_map<int64_t, size_t>& cofacet_bindex_to_list_index;

    MatchingContext (BipartiteGraph& bi_graph,
                     std::vector<std::vector<int64_t>>& binomial_table,
                     std::vector<std::pair<int64_t, double>>& sorted_facets,
                     std::vector<std::pair<int64_t, double>>& sorted_cofacets,
                     const robin_hood::unordered_map<int64_t, size_t>& facet_hash = emptyMap(),
                     const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash = emptyMap())
    : graph(bi_graph),
      binomial_table(binomial_table),
      sorted_facets(sorted_facets),
      sorted_cofacets(sorted_cofacets),
      facet_bindex_to_list_index(facet_hash),
      cofacet_bindex_to_list_index(cofacet_hash)
      {}

    MatchingContext (BipartiteGraph& bi_graph,
                     std::vector<std::vector<int64_t>>& binomial_table,
                     std::vector<std::pair<int64_t, double>>& sorted_facets,
                     std::vector<std::pair<int64_t, double>>& sorted_cofacets,
                     const robin_hood::unordered_map<int64_t, size_t>& facet_hash,
                     const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash,
                     const size_t num_points,
                     const size_t current_dim)
    : graph(bi_graph),
      binomial_table(binomial_table),
      sorted_facets(sorted_facets),
      sorted_cofacets(sorted_cofacets),
      facet_bindex_to_list_index(facet_hash),
      cofacet_bindex_to_list_index(cofacet_hash),
      npts(num_points),
      dim(current_dim)
      {}
    
};