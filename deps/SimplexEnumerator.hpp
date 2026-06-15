#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "robin_hood.h"

#ifdef BUILD_ALPHA_COMPLEX
#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#endif

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

#ifdef BUILD_ALPHA_COMPLEX
    std::vector<std::pair<int64_t, double>> getSortedAlphaCells(const std::vector<std::vector<int64_t>> &binomial_table,
                                                                std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t> &vertex_handle_index,
                                                                CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> &delaunay_d, const size_t dim, double maxeps);
#endif

    // Geometric enumeration for complex with pseudo-vertices (virtual vertices).
    // Each pv represents a set of original vertices, and simplex weight is the smallest
    // clique realization weight among representative choices.
    std::vector<std::pair<int64_t, double>> getGeometricVirtualCofacetList(const std::vector<std::pair<int64_t, double>> &sorted_virtual_simplex_list,
                                                                           const std::vector<size_t> &active_vertices,
                                                                           const std::vector<std::unordered_set<size_t>> &virtual_vertex_indices,
                                                                           const robin_hood::unordered_map<uint64_t, double> &virtual_distance_hash,
                                                                           const size_t dim, const double maxeps, const int threadnum);


private:
    const DistMatType &dist_mat_;
    const std::vector<std::vector<int64_t>> &binomial_table_;

#ifdef BUILD_ALPHA_COMPLEX
    double getAlphaSimplexWeight(const std::vector<size_t> &alpha_simplex);
#endif

    static constexpr size_t REP_STRIDE_ = 64;
    static constexpr size_t UNCHOSEN_ = std::numeric_limits<size_t>::max();

    struct EdgeRecord
    {
        double weight;
        uint8_t virtualidx0, virtualidx1;
        uint8_t localidx0, localidx1;

        bool operator<(const EdgeRecord &edge) const { return weight < edge.weight; }
    };

    struct WitnessWorkspace
    {
        std::vector<size_t> facet_vertices;
        std::vector<size_t> cofacet_vertices;
        std::vector<std::vector<size_t>> singleton_slots;
        std::vector<const std::vector<size_t> *> rep_ptrs;
        std::vector<EdgeRecord> facet_edges;
        std::vector<EdgeRecord> covt_edges;
        std::vector<uint64_t> adj_slab;
        std::vector<uint64_t> candidate_mask;
        std::vector<std::vector<uint64_t>> candidate_stack;
        std::vector<size_t> chosen;
    };

    void prepareFacetWitnessContext(WitnessWorkspace &ws,
                                    const std::vector<size_t> &facet_labels,
                                    const std::vector<std::vector<size_t>> &pv_rep_lists,
                                    const size_t originalvtnum, const double maxeps) const;

    void prepareCovtWitnessGroup(WitnessWorkspace &ws, const size_t covt, const size_t facet_label_count,
                                 const std::vector<std::vector<size_t>> &pv_rep_lists,
                                 const size_t originalvtnum, const double maxeps) const;

    double getGeometricVirtualSimplexWeight(WitnessWorkspace &ws, const size_t groups,
                                            const double lower_bound, const double maxeps) const;

    bool findCliqueRecursive(const uint64_t *adj_slab, const size_t groups, WitnessWorkspace &ws,
                             const std::vector<uint64_t> &candidate, const size_t chosen_count) const;
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
