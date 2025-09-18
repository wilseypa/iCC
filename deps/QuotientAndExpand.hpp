#pragma once

#include <unordered_set>

#include <algorithm>

#include "robin_hood.h"

#include "DistanceMatrix.hpp"
#include "BipartiteGraph.hpp"
#include "MatchingContext.hpp"


template <typename DistMatType>
class QuotientAndExpand
{
public:
    
    QuotientAndExpand(DistMatType& dist_mat, std::vector<std::vector<int64_t>>& binomial_table, const size_t originalvertexnumber):
    dist_mat_(dist_mat), binomial_table_(binomial_table) {}

    std::vector<std::unordered_set<size_t>> runQuotient(const size_t maxdim, const double initeps, const int threadnumber);

    void runExpand(const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices, const size_t maxdim, const double maxeps, const int threadnumber);

private:

    DistMatType& dist_mat_;
    std::vector<std::vector<int64_t>>& binomial_table_;

    struct UnionFind
    {
        std::vector<size_t> parent;
        std::vector<size_t> rank;

        UnionFind(size_t n) : parent(n), rank(n, 0)
        {
            std::iota(parent.begin(), parent.end(), 0);
        }

        size_t unionFind(size_t x)
        {
            if (parent[x] != x)
            {
                parent[x] = unionFind(parent[x]); // Path compression
            }
            return parent[x];
        }

        bool unionSets(size_t x, size_t y)
        {
            size_t px = unionFind(x);
            size_t py = unionFind(y);
            if (px == py) return false; // Already in the same set

            if (rank[px] < rank[py])
            {
                parent[px] = py;
            }
            else if (rank[px] > rank[py])
            {
                parent[py] = px;
            }
            else
            {
                parent[py] = px;
                rank[px]++;
            }
            return true; // Union successful
        }
    };

    struct EdgeRecord    // edge = (idx0, idx1). edge record for clique finding
    {
        double weight;
        size_t virtualidx0, virtualidx1;
        size_t localidx0, localidx1;

        bool operator < (const EdgeRecord& edge) const
        {
            return weight < edge.weight;
        }
    };


    std::vector<std::unordered_set<size_t>> getVirtualVertexIndices(const size_t maxdim, const double initeps, const int threadnumber);

    std::vector< std::vector<size_t> > extractGradientPaths(const MatchingContext& matching_context, const double minfacetweight);

    std::vector<std::unordered_set<size_t>> getGradientPathVertexSets(const MatchingContext& matching_context, const std::vector< std::vector<size_t> >& gradient_paths, const size_t dim);

    std::vector<std::unordered_set<size_t>> trimVertexSets(std::vector<std::unordered_set<size_t>>& gradient_path_vertex_sets);

    std::vector<size_t> getActiveVertexIndices(const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);    //list of active indices to construct edges and cofaces

    double computeVirtualDistance(const size_t i, const size_t j, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);

    robin_hood::unordered_map<uint64_t, double> getVirtualDistanceHashTable(const std::vector<size_t>& active_vertices, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices, int threadnum);

    double getVirtualDistance(size_t i, size_t j, const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table);

    std::vector<std::pair<int64_t, double>> getVirtualSortedEdge(const std::vector<size_t>& active_vertices,
                                                                 const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table, const double maxeps, int threadnum);

    robin_hood::unordered_map<int64_t, size_t> getVirtualActiveEdgeIndexHashTable(const std::vector<std::pair<int64_t, double>>& sorted_virtual_edge, const size_t virtualvtnum);

    std::vector< std::pair<int64_t, double> > getVirtualCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex_list,
                                                                          const std::vector<size_t>& active_vertices, const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table,
                                                                          const size_t dim, const double maxeps, int threadnum);

    bool findCliqueRecursive(const std::vector< std::vector< std::vector<uint64_t> > >& adj_mask,
                             std::vector<uint64_t>& candidate_mask,
                             std::vector<size_t>& current_clique_local_indices, size_t depth);

    double getGeometricVirtualSimplexWeight(const std::vector<size_t>& simplex_vertices, 
                                   const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices,
                                   size_t dim);
    
    std::vector< std::pair<int64_t, double> > getGeometricVirtualCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex_list,
                                                                             const std::vector<size_t>& active_vertices, 
                                                                             const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices,
                                                                             const size_t dim, const double maxeps, int threadnum);                                                               
    
    void buildInterface(BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& sorted_cofacet_list, 
                        const robin_hood::unordered_map<int64_t, size_t>& active_facet_index_hash_table, const size_t dim);


};


extern template class QuotientAndExpand<NormalDistMat>;
// extern template class QuotientAndExpand<SparseDistMat>;

