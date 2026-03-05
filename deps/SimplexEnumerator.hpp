#pragma once

#include <vector>
#include <utility>

#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include "SimplexUtility.hpp"
#include "DistanceMatrix.hpp"

// Forward-declare complex type tags
struct VR;
struct Alpha;

template <typename DistMatType>
class SimplexEnumerator
{
public:
    SimplexEnumerator(const DistMatType &dist_mat, const std::vector<std::vector<int64_t>> &binomial_table)
        : dist_mat_(dist_mat), binomial_table_(binomial_table) {}

    std::vector<std::pair<int64_t, double>> getSortedVREdges(const double maxeps);

    std::vector<std::pair<int64_t, double>> getSortedVRCofacets(const std::vector<std::pair<int64_t, double>> &sorted_simplex, const size_t dim, const double maxeps, const int threadnum);

    std::vector<std::pair<int64_t, double>> getSortedAlphaCells(const std::vector<std::vector<int64_t>> &binomial_table,
                                                                std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t> &vertex_handle_index,
                                                                CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> &delaunay_d, const size_t dim, double maxeps);

    // Geometric enumeration for complex with pseudo-vertices (virtual vertices).
    // Each pv represents a set of original vertices, and simplex weight is the smallest
    // clique realization weight among representative choices.
    std::vector<std::pair<int64_t, double>> getGeometricVirtualCofacetList(const std::vector<std::pair<int64_t, double>> &sorted_virtual_simplex_list,
                                                                           const std::vector<size_t> &active_vertices,
                                                                           const std::vector<std::unordered_set<size_t>> &virtual_vertex_indices,
                                                                           const size_t dim, const double maxeps, const int threadnum);


private:
    const DistMatType &dist_mat_;
    const std::vector<std::vector<int64_t>> &binomial_table_;

    double getAlphaSimplexWeight(const std::vector<size_t> &alpha_simplex);

    struct EdgeRecord
    {
        double weight;
        uint8_t virtualidx0, virtualidx1;
        uint8_t localidx0, localidx1;

        bool operator<(const EdgeRecord &edge) const { return weight < edge.weight; }
    };

    bool findCliqueRecursive(const std::vector<std::vector<std::vector<uint64_t>>>& adj_mask,
                             std::vector<uint64_t>& candidate_mask,
                             std::vector<size_t>& current_clique_local_indices, size_t depth);

    double getGeometricVirtualSimplexWeight(const std::vector<size_t>& simplex_vertices,
                                            const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices, size_t dim);
};

// Explicit template instantiations
// loss of generality, but avoid code bloat in header
extern template class SimplexEnumerator<NormalDistMat>;
// extern template class SimplexEnumerator<SparseDistMat>;    //to do

// template <typename DistMatType>
// template <typename ComplexType>
// std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum)
// {
//     if constexpr (std::is_same_v<ComplexType, VR>)
//     {
//         return this->getSortedVRCofacets(sorted_simplex, dim, maxeps, threadnum);
//     }
//     else if constexpr (std::is_same_v<ComplexType, Alpha>)
//     {
//         return this->getSortedAlphaCofacets();
//     }
//     else
//     {
//         throw std::invalid_argument("Unsupported ComplexType for SimplexEnumerator.");
//     }
// }

// template <typename DistMatType>
// template <typename ComplexType>
// std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedEdges(const double maxeps)
// {
//     if constexpr (std::is_same_v<ComplexType, VR>)
//     {
//         return this->getSortedVREdges(maxeps);
//     }
//     else if constexpr (std::is_same_v<ComplexType, Alpha>)
//     {
//         throw std::runtime_error("Alpha complex edges enumeration is not implemented.");
//     }
//     else
//     {
//         throw std::invalid_argument("Unsupported ComplexType for SimplexEnumerator.");
//     }
// }