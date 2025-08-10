#pragma once

#include <type_traits>
#include <vector>
#include <utility>

#include "SimplexUtility.hpp"
#include "DistanceMatrix.hpp"

template <typename DistMatType>
class SimplexEnumerator
{
public:
    SimplexEnumerator(const DistMatType& dist_mat, const std::vector<std::vector<int64_t>>& binomial_table)
        : dist_matrix(dist_mat), binomial_table(binomial_table) {}


    std::vector<std::pair<int64_t, double>> getSortedEdges(const double maxeps);

    template <typename ComplexType>
    std::vector<std::pair<int64_t, double>> getSortedCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum);

    
private:
    const NormallDistMat& dist_matrix;
    const std::vector<std::vector<int64_t>>& binomial_table;

    std::vector<std::pair<int64_t, double>> getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum);

    std::vector<std::pair<int64_t, double>> getSortedAlphaCofacets();
};

template <typename DistMatType>
template <typename ComplexType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum)
{
    if constexpr (std::is_same_v<ComplexType, VR>)
    {
        return this->getSortedVRCofacets(sorted_simplex, dim, maxeps, threadnum);
    }
    else if constexpr (std::is_same_v<ComplexType, Alpha>)
    {
        return this->getSortedAlphaCofacets();
    }
    else
    {
        throw std::invalid_argument("Unsupported ComplexType for SimplexEnumerator.");
    }
}