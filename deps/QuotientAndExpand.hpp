#pragma once

#include <unordered_set>

#include <algorithm>

#include "DistanceMatrix.hpp"
#include "BipartiteGraph.hpp"
#include "MatchingContext.hpp"


template <typename DistMatType>
class QuotientAndExpand
{
public:
    
    QuotientAndExpand(DistMatType& dist_mat, const std::vector<std::vector<int64_t>>& binomial_table, const size_t npts):
    dist_mat_(dist_mat), binomial_table_(binomial_table), npts_(npts) {}

    std::vector<std::unordered_set<size_t>> getVirtualVertexIndices(const size_t maxdim, const double initeps, const int threadnumber);

private:

    size_t npts_;
    DistMatType& dist_mat_;
    std::vector<std::vector<int64_t>>& binomial_table_;

    std::vector< std::vector<size_t> > extractGradientPahts(const MatchingContext& matching_context, const double minfacetweight);

    std::vector<std::unordered_set<size_t>> getGradientPathVertexSets(const MatchingContext& matching_context, const std::vector< std::vector<size_t> >& gradient_paths, const size_t dim);

    struct GradientPathUnionFind
    {
        std::vector<size_t> parent;
        std::vector<size_t> rank;

        GradientPathUnionFind(size_t n) : parent(n), rank(n, 0)
        {
            std::iota(parent.begin(), parent.end(), 0);
        }

        size_t unionFind(size_t x)
        {
            if (parent[x] != x)
            {
                parent[x] = find(parent[x]); // Path compression
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

};