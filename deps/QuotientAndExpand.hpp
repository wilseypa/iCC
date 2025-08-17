#pragma once

#include <unordered_set>

#include "DistanceMatrix.hpp"
#include "SimplexEnumerator.hpp"
#include "BipartiteGraph.hpp"
#include "MatchingContext.hpp"
#include "MaximumMorseMatching.hpp"


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
    std::vector<std::vector<int64_t>> binomial_table_;
};