#pragma once

#include <cstdint>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <vector>

#include "robin_hood.h"

#include "DistanceMatrix.hpp"
#include "BipartiteGraph.hpp"
#include "MatchingContext.hpp"
#include "MaximumMorseMatching.hpp"

template <typename DistMatType>
class QuotientAndExpand
{
public:
    QuotientAndExpand(DistMatType& dist_mat, std::vector<std::vector<int64_t>>& binomial_table, const size_t originalvertexnumber) : dist_mat_(dist_mat), binomial_table_(binomial_table) {}

    void runPiecewisePH(const std::vector<double>& eps_breaks, const size_t maxdim, const int thread_number,
                        const double pv_cap_scale, const double pv_min_separation = 0.0,
                        const bool verbose = false);

    //legacy QE

    std::vector<std::unordered_set<size_t>> runQuotient(const size_t maxdim, const double initeps, const int thread_number);

    void runExpand(const std::vector<std::unordered_set<size_t>>& pv_index_sets, const size_t maxdim, const double maxeps, const int thread_number);

private:
    static constexpr int MAX_SIZE_ = 64; // PV size cap
    static constexpr size_t MAX_FFI_PACKED_LABELS_ = sizeof(uint64_t);

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
            if (px == py)
                return false; // Already in the same set

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

    struct WindowState
    {
        size_t original_vertex_number = 0;

        std::vector<std::unordered_set<size_t>> pv_flat_index_set_list;
        // Kept in the same order as pv_flat_index_set_list.
        std::vector<double> pv_diameter_list;

        // active labels in the current window:
        // original vertices: 0 .. original_vertex_number-1
        // PV labels: original_vertex_number + i
        std::vector<size_t> active_label_list;

        WindowState() = default;

        explicit WindowState(size_t npts): original_vertex_number(npts), pv_flat_index_set_list(), pv_diameter_list(), active_label_list(npts)
        {
            std::iota(active_label_list.begin(), active_label_list.end(), 0);
        }
    };

    struct SelectedPV
    {
        // Flattened to original vertex indices and its max pairwise distance.
        std::unordered_set<size_t> flat_index_set;
        double diameter = 0.0;
    };

    struct QuotientEdgeData
    {
        robin_hood::unordered_map<uint64_t, double> label_distance_hash;
        std::vector<std::pair<int64_t, double>> sorted_edges;
    };

    std::vector<std::unordered_set<size_t>> runWindow(const WindowState& win_state, const size_t maxdim, const double eps_lo, const double eps_hi,
                                                      const int thread_number, const bool collect_pv, const bool verbose);

    std::vector<SelectedPV> trimPVCandidates(const WindowState& win_state, const std::vector<std::unordered_set<size_t>>& raw_label_sets,
                                             const double eps_hi, const double pv_cap_scale, const double min_separation);

    std::unordered_set<size_t> flattenLabelSet(const WindowState& win_state, const std::unordered_set<size_t>& raw_label_set);

    double getMaxPairwiseDistance(const std::unordered_set<size_t>& index_set) const;

    void rebuildWindowState(WindowState& win_state, std::vector<SelectedPV>&& new_pv_list);

    void collectProtectedIndices(const MatchingContext& matching_context,
                                 const MaximumMorseMatching::MatchSupportInfo& match_support_info,
                                 std::unordered_set<size_t>& protected_indices);

    std::vector<std::unordered_set<size_t>> getNonMergingPVSupport(const MatchingContext& matching_context,
                                                                   const std::vector<std::vector<size_t>>& raw_pv_support_cofacet_indices,
                                                                   const std::unordered_set<size_t>& protected_indices,
                                                                   const size_t origin_vt_num,
                                                                   const bool verbose);

    double computeLabelDistance(const size_t i, const size_t j, const std::vector<std::unordered_set<size_t>>& pv_index_sets);

    QuotientEdgeData buildQuotientEdges(const std::vector<size_t>& active_labels,
                                         const std::vector<std::unordered_set<size_t>>& pv_index_sets,
                                         const double maxeps, int threadnum);

    // Verbose-only top-interface diagnostics. Realizations are compared exactly; the
    // multi-coordinate carrier-clique test is sufficient, but not necessary, for a filler.
    void reportFalseFacetIdentificationStats(const WindowState& win_state,
                                             const std::vector<std::pair<int64_t, double>>& facet_list,
                                             const std::vector<std::pair<int64_t, double>>& cofacet_list,
                                             const robin_hood::unordered_map<int64_t, uint64_t>& facet_pv_realizations,
                                             const robin_hood::unordered_map<int64_t, uint64_t>& cofacet_pv_realizations,
                                             const size_t interface_dim, const double eps_lo, const double eps_hi,
                                             const int threadnum);

    robin_hood::unordered_map<int64_t, size_t> getQuotientActiveEdgeIndexHashTable(const std::vector<std::pair<int64_t, double>>& sorted_quotient_edge, const size_t pvnum);

    //legacy QE
    std::vector<std::unordered_set<size_t>> getPVIndexSets(const size_t maxdim, const double initeps, const int thread_number);

    std::vector<std::unordered_set<size_t>> trimIndexSets(std::vector<std::unordered_set<size_t>>& pv_support_vertex_sets, const double initeps);

    std::vector<size_t> getActiveLabelIndices(const std::vector<std::unordered_set<size_t>>& pv_index_sets); // labels used to construct quotient edges and cofacets
};

extern template class QuotientAndExpand<NormalDistMat>;
// extern template class QuotientAndExpand<SparseDistMat>;
