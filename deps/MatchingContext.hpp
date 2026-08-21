#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "BipartiteGraph.hpp"
#include "CofacetIndex.hpp"
#include "robin_hood.h"

struct NormalDistMat;

struct ApparentPairSearchConfig
{
    enum class Mode : uint8_t
    {
        Disabled,
        NormalVR,
        QuotientVR
    };

    Mode mode = Mode::Disabled;
    const NormalDistMat* dist_mat = nullptr;

    // Quotient-only data. Original labels are [0, original_vertex_count).
    const robin_hood::unordered_map<uint64_t, double>* pv_label_distance_hash = nullptr;
    const std::vector<uint8_t>* active_label_mask = nullptr;
    size_t original_vertex_count = 0;
};

namespace    //for explicit/legacy
{
  inline const robin_hood::unordered_map<int64_t, size_t>& emptyMap() 
  {
    static const robin_hood::unordered_map<int64_t, size_t> empty_map;
    return empty_map;
  }

  inline const CofacetIndex& emptyCofacetIndex()
  {
    static const CofacetIndex empty_index;
    return empty_index;
  }
}

struct MatchingContext
{
    //default values for explicit/legacy implementation
    const size_t total_label_count = 0;    //total labels: original vertices + virtual/PV labels
    const size_t dim = 0;    //current dim of interface (cofacet dim)

    // Disabled by default so explicit, approximate, and other callers retain the
    // original full-cofacet matching path unless they opt in explicitly.
    ApparentPairSearchConfig apparent_pair_search{};

    BipartiteGraph& graph;

    // References to data needed by the algorithms
    std::vector<std::vector<int64_t>>& binomial_table;
    std::vector<std::pair<int64_t, double>>& sorted_facets;
    std::vector<std::pair<int64_t, double>>& sorted_cofacets;

    // Lookup structures required for implicit graph queries
    const robin_hood::unordered_map<int64_t, size_t>& facet_bindex_to_list_index;
    const CofacetIndex& cofacet_bindex_to_list_index;

    MatchingContext (BipartiteGraph& bi_graph,
                     std::vector<std::vector<int64_t>>& binomial_table,
                     std::vector<std::pair<int64_t, double>>& sorted_facets,
                     std::vector<std::pair<int64_t, double>>& sorted_cofacets,
                     const robin_hood::unordered_map<int64_t, size_t>& facet_hash = emptyMap(),
                     const CofacetIndex& cofacet_index = emptyCofacetIndex())
    : total_label_count(0),
      dim(0),
      graph(bi_graph),
      binomial_table(binomial_table),
      sorted_facets(sorted_facets),
      sorted_cofacets(sorted_cofacets),
      facet_bindex_to_list_index(facet_hash),
      cofacet_bindex_to_list_index(cofacet_index)
      {}

    MatchingContext (BipartiteGraph& bi_graph,
                     std::vector<std::vector<int64_t>>& binomial_table,
                     std::vector<std::pair<int64_t, double>>& sorted_facets,
                     std::vector<std::pair<int64_t, double>>& sorted_cofacets,
                     const robin_hood::unordered_map<int64_t, size_t>& facet_hash,
                     const CofacetIndex& cofacet_index,
                     const size_t total_label_count_in,
                     const size_t current_dim)
    : total_label_count(total_label_count_in),
      dim(current_dim),
      graph(bi_graph),
      binomial_table(binomial_table),
      sorted_facets(sorted_facets),
      sorted_cofacets(sorted_cofacets),
      facet_bindex_to_list_index(facet_hash),
      cofacet_bindex_to_list_index(cofacet_index)
      {}
    
};
