#include "omp.h"

#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVREdges(const double maxeps)
{
    std::vector<std::pair<int64_t, double>> sorted_edge;

    size_t npt = dist_mat_.dist_matrix.size();

    for (size_t i = 0; i < npt - 1; i++)
    {
        for (size_t j = i + 1; j < npt; j++)
        {
            double weight = dist_mat_.getDistance(i, j);
            if (weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(binomial_table, j, i);
                sorted_edge.push_back(std::make_pair(bindex, weight));
            }
        }
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_edge);

    return sorted_edge;
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum)
{
    //dim == simplex dimension == cofacet dimension - 1
    std::vector<std::pair<int64_t, double>> cofacet_list;

    size_t npt = dist_mat_.dist_matrix.size();

    std::vector< std::vector< std::pair<size_t, double> > > thread_workspace(sorted_simplex.size(), std::vector<std::pair<size_t, double>>());

    omp_set_num_threads(threadnum);

#pragma omp parallel for
    for (size_t i = 0; i < sorted_simplex.size(); i++)
    {
        int64_t bindex = sorted_simplex[i].first;
        double originalweight = sorted_simplex[i].second;
        std::vector<size_t> simplex_vt = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npt, dim);

        //cofacet of {i, j} = {i, j, ...}. i > j
        // auto minidx = *std::min_element(simplex_vt.begin(), simplex_vt.end());
        auto minidx = simplex_vt.back();
        for (size_t j = 0; j < minidx; j++)
        {
            auto idx = *std::max_element(simplex_vt.begin(), simplex_vt.end(), [this, j](auto first, auto second)
                                                                                { return this->dist_matrix[j][first] < this->dist_matrix[j][second]; });
            double weight = dist_matrix[j][idx];    //new facet weight
            double cofacetweight = (weight > originalweight) ? weight : originalweight; 

            if (cofacetweight < maxeps) thread_workspace[i].emplace_back(j, cofacetweight);
        }
    }

    for (size_t i = 0; i < sorted_simplex.size(); i++)
    {
        int64_t bindex = sorted_simplex[i].first;
        std::vector<size_t> simplex_vt = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npt, dim);
        //left shift bindex. left shift one position for each simplex vt
        bindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vt, 1);
        for (auto& idx_weight_pair: thread_workspace[i])
        {
            auto j = idx_weight_pair.first;
            auto w = idx_weight_pair.second;
            //push new bindex = shifted bin
            cofacet_list.emplace_back(bindex + j, w);
        }
    }

    SimplexUtility::sortSimplexByWeightThenIndex(cofacet_list);

    return cofacet_list;
}


template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedAlphaCofacets()
{
    // This function is a placeholder for the Alpha cofacet enumeration.
    throw std::runtime_error("Alpha complex has not implemented yet.");
}


