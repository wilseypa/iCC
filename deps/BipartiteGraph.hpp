#pragma once

#include <vector>

struct ImplicitConstructionTag
{
}; // tag for implicit bipartite graph

/**
 * @struct BipartiteGraph
 * @brief A simple data structure to represent the bipartite graph interface.
 * This struct holds the adjacency list and the matching results.
 */
struct BipartiteGraph
{
  size_t unodes; // Number of nodes on the left (e.g., cofacets)
  size_t vnodes; // Number of nodes on the right (e.g., facets)

  // Adjacency list for the graph. adj_list[0...u-1] are left nodes,
  // adj_list[u...u+v-1] are right nodes.
  std::vector<std::vector<size_t>> adj_list;

  // Stores the matching result. match_list[i] = j means i is matched with j.
  std::vector<int64_t> match_list;

  explicit BipartiteGraph(size_t u, size_t v)
      : unodes(u),
        vnodes(v),
        adj_list(u + v),
        match_list(u + v, -1) {}

  BipartiteGraph(size_t u, size_t v, ImplicitConstructionTag /* tag */)
      : unodes(u),
        vnodes(v),
        match_list(u + v, -1) {}

  void updateDimension(size_t cofacetnum, size_t facetnum)
  {
    unodes = cofacetnum;
    vnodes = facetnum;
    adj_list.assign(unodes + vnodes, std::vector<size_t>());
    match_list.assign(unodes + vnodes, -1);
  }
};