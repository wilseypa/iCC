#include <unordered_map>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>

#include "omp.h"

#include "SimplexUtility.hpp"
#include "SimplexEnumerator.hpp"
#include "MaximumMorseMatching.hpp"
#include "QuotientAndExpand.hpp"

namespace
{
    void printSimplexLabels(const std::vector<std::vector<int64_t>>& binomial_table,
                            const char* prefix,
                            const int64_t simplex_bindex,
                            const size_t total_label_count,
                            const size_t simplex_dim)
    {
        std::cout << "    " << prefix << " labels: ";
        if (simplex_bindex < 0)
        {
            std::cout << "(none)\n";
            return;
        }

        const auto simplex_labels = SimplexUtility::getSimplexVertices(
            binomial_table, simplex_bindex, total_label_count, simplex_dim);

        if (simplex_labels.empty())
        {
            std::cout << "(empty)";
        }
        else
        {
            for (size_t i = 0; i < simplex_labels.size(); ++i)
            {
                if (i > 0) std::cout << ' ';
                std::cout << simplex_labels[i];
            }
        }
        std::cout << '\n';
    }

    void printPersistentPairs(
        const std::vector<MaximumMorseMatching::PersistentPairInfo>& persistent_pairs,
        const std::vector<std::vector<int64_t>>& binomial_table,
        const double eps_lo,
        const size_t dim,
        const size_t maxdim,
        const size_t total_label_count,
        const bool verbose)
    {
        bool printed_any = false;
        for (const auto& pair_info : persistent_pairs)
        {
            const double facetweight = pair_info.facet_weight;
            const double cofacetweight = pair_info.cofacet_weight;
            if (!(cofacetweight >= eps_lo || cofacetweight < 0)) continue;

            if (verbose)
            {
                std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;

                if (dim == maxdim)
                {
                    printSimplexLabels(binomial_table, "birth facet", pair_info.facet_bindex,
                                       total_label_count, dim - 1);
                    printSimplexLabels(binomial_table, "death cofacet", pair_info.cofacet_bindex,
                                       total_label_count, dim);
                }
            }
            else
            {
                std::cout << dim << ", " << facetweight << ", " << cofacetweight << '\n';
            }

            printed_any = true;
        }

        if (verbose && !printed_any)
            std::cout << "  (no new interval or surviving interval from previous eps range)" << std::endl;
    }
}

template <typename DistMatType>
double QuotientAndExpand<DistMatType>::getMaxPairwiseDistance(const std::unordered_set<size_t>& index_set) const
{
    if (index_set.size() < 2)
        return 0.0;

    double maxdist = 0.0;
    for (auto first = index_set.begin(); first != index_set.end(); ++first)
    {
        auto second = first;
        ++second;
        for (; second != index_set.end(); ++second)
        {
            const double distance = dist_mat_.getDistance(*first, *second);
            maxdist = std::max(maxdist, distance);
        }
    }

    return maxdist;
}

template <typename DistMatType>
std::vector<std::vector<size_t>> QuotientAndExpand<DistMatType>::getPVRepresentativeLists(
    const WindowState& win_state) const
{
    std::vector<std::vector<size_t>> pv_rep_lists;
    pv_rep_lists.reserve(win_state.pv_list.size());

    for (const auto& pv : win_state.pv_list)
    {
        pv_rep_lists.emplace_back(pv.flat_index_set.begin(), pv.flat_index_set.end());
        std::sort(pv_rep_lists.back().begin(), pv_rep_lists.back().end());
    }

    return pv_rep_lists;
}

template <typename DistMatType>
std::vector<std::vector<size_t>> QuotientAndExpand<DistMatType>::getPVRepresentativeLists(
    const std::vector<std::unordered_set<size_t>>& pv_index_sets) const
{
    std::vector<std::vector<size_t>> pv_rep_lists;
    pv_rep_lists.reserve(pv_index_sets.size());

    for (const auto& pv_index_set : pv_index_sets)
    {
        pv_rep_lists.emplace_back(pv_index_set.begin(), pv_index_set.end());
        std::sort(pv_rep_lists.back().begin(), pv_rep_lists.back().end());
    }

    return pv_rep_lists;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::runPiecewisePH(const std::vector<double>& eps_breaks, const size_t maxdim, const int threadnumber,
                                                    const double pv_cap_scale, const double pv_min_separation, const bool verbose)
{
    WindowState win_state(dist_mat_.getVertexNumber());
    const double min_separation = pv_min_separation * eps_breaks.back();

    // from eps_0 = 0 to eps_max
    std::vector<double> full_eps_list;
    full_eps_list.reserve(eps_breaks.size() + 1);
    full_eps_list.push_back(0.0);
    full_eps_list.insert(full_eps_list.end(), eps_breaks.begin(), eps_breaks.end());

    for (size_t i = 1; i < full_eps_list.size(); ++i)
    {
        const double eps_lo = full_eps_list[i - 1];
        const double eps_hi = full_eps_list[i];

        const bool collect_pv = (i + 1 < full_eps_list.size());

        auto untrimmed_pv_label_sets = runWindow(win_state, maxdim, eps_lo, eps_hi, threadnumber, collect_pv, verbose);

        if (!collect_pv) break;

        auto new_pv_list = trimPVCandidates(win_state, untrimmed_pv_label_sets, eps_hi, pv_cap_scale, min_separation);
        const size_t new_pv_count = new_pv_list.size();

        rebuildWindowState(win_state, std::move(new_pv_list));

        if (verbose)
        {
            std::cout << "after eps range "<<eps_lo<< "  " <<eps_hi<< "  new pv number = "<<new_pv_count
                      <<"  total pv number = "<<win_state.pv_list.size()<<std::endl;

            std::cout << "pv flat index set statistics:" << std::endl;
            if (win_state.pv_list.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (size_t pv_idx = 0; pv_idx < win_state.pv_list.size(); ++pv_idx)
                {
                    const auto& pv = win_state.pv_list[pv_idx];
                    std::cout << "  [" << pv_idx << "] size = " << pv.flat_index_set.size()
                              << "  max diameter = " << pv.diameter << '\n';
                }
                std::cout << std::endl;
            }
        }
    }
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runWindow(const WindowState& win_state, const size_t maxdim, const double eps_lo, const double eps_hi,
                                                                                  const int threadnumber, const bool collect_pv, const bool verbose)
{
    if (!verbose)
    {
        std::cout << "window lower bound, window upper bound\n";
        std::cout << eps_lo << ", " << eps_hi << '\n';
        std::cout << "dimension, birth weight, death weight\n";
    }

    const size_t original_vt_num = win_state.original_vertex_number;
    const size_t pv_num = win_state.pv_list.size();
    const size_t total_label_count = original_vt_num + pv_num;

    SimplexUtility::updateBinomialTable(binomial_table_, original_vt_num, pv_num, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);

    // Build one deterministic ordering per window for enumeration and packed FFI witnesses.
    const auto pv_rep_lists = getPVRepresentativeLists(win_state);
    const auto& active_labels = win_state.active_label_list;

    auto quotient_edge_data = buildQuotientEdges(active_labels, pv_rep_lists, eps_hi, threadnumber);
    auto pv_label_distance_hash = std::move(quotient_edge_data.pv_label_distance_hash);
    auto sorted_quotient_simplex = std::move(quotient_edge_data.sorted_edges);

    auto active_simplex_hash = getQuotientActiveEdgeIndexHashTable(sorted_quotient_simplex, pv_num);

    // Matching sweeps the complete label range. This mask cheaply rejects labels
    // that are inactive in the current quotient complex before distance/hash work.
    std::vector<uint8_t> active_label_mask(total_label_count, 0);
    for (const size_t label : active_labels)
        active_label_mask[label] = 1;

    const bool ffi_stats_requested = verbose && (pv_num > 0) && (maxdim >= 3);
    const bool collect_ffi_stats = ffi_stats_requested && (maxdim < MAX_FFI_PACKED_LABELS_);

    if (ffi_stats_requested && !collect_ffi_stats)
    {
        std::cout << "  [ffi] diagnostics skipped: packed witness realizations support max dimension "
                  << (MAX_FFI_PACKED_LABELS_ - 1) << '\n';
    }

    struct FfiRealizationState
    {
        robin_hood::unordered_map<int64_t, uint64_t> facet_realization_hash;
        robin_hood::unordered_map<int64_t, uint64_t> cofacet_realization_hash;
    };
    std::unique_ptr<FfiRealizationState> ffi_realizations;
    const auto enumerate_initial_cofacets = [&]()
    {
        if (collect_ffi_stats)
        {
            ffi_realizations = std::make_unique<FfiRealizationState>();
            return simplex_enumerator.getGeometricCofacetListWithRealizations(sorted_quotient_simplex, 
                                                                              active_labels, pv_rep_lists, pv_label_distance_hash,
                                                                              1, eps_hi, threadnumber, ffi_realizations->cofacet_realization_hash);
        }

        return simplex_enumerator.getGeometricCofacetList(sorted_quotient_simplex, 
                                                          active_labels, pv_rep_lists, pv_label_distance_hash,
                                                          1, eps_hi, threadnumber);
    };
    auto sorted_quotient_cofacet = enumerate_initial_cofacets();
    
    BipartiteGraph bi_graph(1, 1, ImplicitConstructionTag{});

    std::vector<MaximumMorseMatching::PersistentPairInfo> dim_persistent_pairs;

    std::vector<std::unordered_set<size_t>> untrimmed_pv_label_sets;

    MaximumMorseMatching morse_matching;
    std::unordered_set<size_t> protected_indices;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        const bool is_top_dimension = (dim == maxdim);

        bi_graph.updateDimensionImplicit(sorted_quotient_cofacet.size(), sorted_quotient_simplex.size());

        auto cofacet_index = SimplexUtility::getCofacetIndex(sorted_quotient_cofacet);

        if (collect_ffi_stats && is_top_dimension)
        {
            reportFalseFacetIdentificationStats(win_state, pv_rep_lists,
                                                sorted_quotient_simplex, sorted_quotient_cofacet,
                                                ffi_realizations->facet_realization_hash,
                                                ffi_realizations->cofacet_realization_hash,
                                                dim, eps_lo, eps_hi, threadnumber);
        }

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_quotient_simplex, sorted_quotient_cofacet,
                                         active_simplex_hash, cofacet_index, total_label_count, dim);

        if constexpr (std::is_same_v<DistMatType, NormalDistMat>)
        {
            matching_context.apparent_pair_search.mode = ApparentPairSearchConfig::Mode::QuotientVR;
            matching_context.apparent_pair_search.dist_mat = &dist_mat_;
            matching_context.apparent_pair_search.pv_label_distance_hash = &pv_label_distance_hash;
            matching_context.apparent_pair_search.active_label_mask = &active_label_mask;
            matching_context.apparent_pair_search.original_vertex_count = original_vt_num;
        }

        MaximumMorseMatching::MatchSupportInfo dim_match_support_info;
        if (collect_pv)
        {
            dim_match_support_info = morse_matching.implicitMatchAndCollectSupportInfo(matching_context,
                                                                                       dim_persistent_pairs,
                                                                                       is_top_dimension);
            collectProtectedIndices(matching_context, dim_match_support_info, protected_indices);
        }
        else
        {
            morse_matching.implicitMatch(matching_context, dim_persistent_pairs);
        }

        if (verbose)
        {
            std::cout << "in eps range "<<eps_lo<< "  " <<eps_hi<< "    dimension = " <<dim
                      << "  cofacet num = " << sorted_quotient_cofacet.size()
                      << "  facet num = " << sorted_quotient_simplex.size() <<'\n'
                      << "   persistent pairs:" << std::endl;
        }

        printPersistentPairs(dim_persistent_pairs, binomial_table_, eps_lo, dim, maxdim, total_label_count, verbose);
        dim_persistent_pairs.clear();

        if (is_top_dimension && collect_pv)
        {
            untrimmed_pv_label_sets = getNonMergingPVSupport(matching_context,
                                                             dim_match_support_info.raw_pv_support_cofacet_indices,
                                                             protected_indices,
                                                             original_vt_num,
                                                             verbose);
        }

        if (!is_top_dimension)
        {
            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_simplex_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_quotient_cofacet);

            // Both structures are dead past this point. Release their backing
            // storage before the next enumeration and merge.
            cofacet_index.release();
            decltype(sorted_quotient_simplex){}.swap(sorted_quotient_simplex);

            if (collect_ffi_stats)
            {
                std::swap(ffi_realizations->facet_realization_hash,
                          ffi_realizations->cofacet_realization_hash);
                sorted_quotient_simplex = simplex_enumerator.getGeometricCofacetListWithRealizations(sorted_quotient_cofacet, 
                                                                                                     active_labels, pv_rep_lists, pv_label_distance_hash,
                                                                                                     dim, eps_hi, threadnumber, ffi_realizations->cofacet_realization_hash);
            }
            else
            {
                sorted_quotient_simplex = simplex_enumerator.getGeometricCofacetList(sorted_quotient_cofacet, 
                                                                                     active_labels, pv_rep_lists, pv_label_distance_hash,
                                                                                     dim, eps_hi, threadnumber);
            }
            std::swap(sorted_quotient_simplex, sorted_quotient_cofacet);
        }
    }

    if (!verbose)
        std::cout << '\n';

    return untrimmed_pv_label_sets;
}

template <typename DistMatType>
std::vector<typename QuotientAndExpand<DistMatType>::SelectedPV> 
QuotientAndExpand<DistMatType>::trimPVCandidates(const WindowState& win_state, const std::vector<std::unordered_set<size_t>>& raw_label_sets,
                                                 const double eps_hi, const double pv_cap_scale, const double min_separation)
{
    std::vector<SelectedPV> selected_new_pv_list;

    const auto minCrossDistance = [this](const std::unordered_set<size_t>& set_a, const std::unordered_set<size_t>& set_b)
    {
        double min_distance = std::numeric_limits<double>::infinity();
        for (const auto vertex_a : set_a)
        {
            for (const auto vertex_b : set_b)
                min_distance = std::min(min_distance, dist_mat_.getDistance(vertex_a, vertex_b));
        }
        return min_distance;
    };

    const auto isSeparated = [&](const std::unordered_set<size_t>& flat_index_set)
    {
        if (min_separation <= 0.0)
            return true;

        for (const auto& carried_pv : win_state.pv_list)
        {
            if (minCrossDistance(flat_index_set, carried_pv.flat_index_set) < min_separation)
                return false;
        }

        for (const auto& selected_pv : selected_new_pv_list)
        {
            if (minCrossDistance(flat_index_set, selected_pv.flat_index_set) < min_separation)
                return false;
        }

        return true;
    };

    std::unordered_set<size_t> claimed_labels;

    for (auto rit = raw_label_sets.rbegin(); rit != raw_label_sets.rend(); ++rit)
    {
        auto& label_set = *rit;
        bool overlap = false;

        for (const auto label: label_set)
        {
            if (claimed_labels.count(label))
            {
                overlap = true;
                break;
            }
        }

        if (!overlap)
        {
            auto flat_index_set = flattenLabelSet(win_state, label_set);

            /********************************** no old pv absorption at this time **********************************/

            if (flat_index_set.size() >= MAX_SIZE_) continue;

            const double diameter = getMaxPairwiseDistance(flat_index_set);
            if (diameter > (pv_cap_scale * eps_hi)) continue;

            if (!isSeparated(flat_index_set)) continue;

            claimed_labels.insert(label_set.begin(), label_set.end());

            selected_new_pv_list.push_back(SelectedPV{std::move(flat_index_set), diameter});
        }
    }

    return selected_new_pv_list;
}

template <typename DistMatType>
std::unordered_set<size_t> QuotientAndExpand<DistMatType>::flattenLabelSet(const WindowState& win_state, const std::unordered_set<size_t>& raw_label_set)
{
    const auto origin_vt_num = win_state.original_vertex_number;

    std::unordered_set<size_t> flat_index_set;

    for (const auto label : raw_label_set)
    {
        if (label < origin_vt_num)    //regular vertex index
        {
            flat_index_set.insert(label);
        }
        else
        {
            const auto& pv_indices = win_state.pv_list[label - origin_vt_num].flat_index_set;
            flat_index_set.insert(pv_indices.begin(), pv_indices.end());
        }
    }

    return flat_index_set;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::rebuildWindowState(WindowState& win_state, std::vector<SelectedPV>&& new_pv_list)
{

    /********************************** no old pv absorption **********************************/

    if (new_pv_list.empty())
        return;

    const size_t origin_vt_num = win_state.original_vertex_number;
    const size_t old_pv_num = win_state.pv_list.size();

    std::vector<bool> is_newly_covered(origin_vt_num, false);
    for (const auto& new_pv : new_pv_list)
    {
        for (const auto idx : new_pv.flat_index_set)
            is_newly_covered[idx] = true;
    }

    std::vector<size_t> active_label_list;
    active_label_list.reserve(win_state.active_label_list.size() + new_pv_list.size());

    for (const auto label : win_state.active_label_list)
    {
        if (label < origin_vt_num && !is_newly_covered[label])
            active_label_list.push_back(label);
    }

    for (const auto label : win_state.active_label_list)
    {
        if (label >= origin_vt_num)
            active_label_list.push_back(label);
    }

    size_t next_pv_label = origin_vt_num + old_pv_num;
    win_state.pv_list.reserve(old_pv_num + new_pv_list.size());
    for (auto& new_pv : new_pv_list)
    {
        win_state.pv_list.push_back(std::move(new_pv));
        active_label_list.push_back(next_pv_label++);
    }

    win_state.active_label_list = std::move(active_label_list);

    return;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::collectProtectedIndices(const MatchingContext& matching_context,
                                                             const MaximumMorseMatching::MatchSupportInfo& match_support_info,
                                                             std::unordered_set<size_t>& protected_indices)
{
    const size_t origin_vt_num = dist_mat_.getVertexNumber();
    const size_t total_label_count = matching_context.total_label_count;
    const size_t facet_dim = matching_context.dim - 1;

    for (const auto protected_facet_list_idx : match_support_info.protected_facet_list_indices)
    {
        const auto facet_bindex = matching_context.sorted_facets[protected_facet_list_idx].first;
        const auto facet_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, facet_bindex, total_label_count, facet_dim);

        for (const auto vertex : facet_vertices)
        {
            if (vertex < origin_vt_num)
                protected_indices.insert(vertex);
        }
    }
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getNonMergingPVSupport(const MatchingContext& matching_context,
                                                                                               const std::vector<std::vector<size_t>>& raw_pv_support_cofacet_indices,
                                                                                               const std::unordered_set<size_t>& protected_indices,
                                                                                               const size_t origin_vt_num,
                                                                                               const bool verbose)
{
    std::vector<std::unordered_set<size_t>> pv_support_label_sets;
    pv_support_label_sets.reserve(raw_pv_support_cofacet_indices.size());

    const size_t total_label_count = matching_context.total_label_count;
    const size_t dim = matching_context.dim;

    for (const auto& pv_support : raw_pv_support_cofacet_indices)
    {
        //need set for trimming
        std::unordered_set<size_t> label_set;

        //do not merge the previous PVs
        bool has_pv = false;

        for (auto cofacetidx : pv_support)
        {
            auto bindex = matching_context.sorted_cofacets[cofacetidx].first;
            auto simplex_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, bindex, total_label_count, dim);

            if (simplex_vertices.front() >= origin_vt_num)
            {
                has_pv = true;
                break;
            }

            label_set.insert(simplex_vertices.begin(), simplex_vertices.end());
        }
        if (has_pv) continue;

        bool has_protected_vertex = false;
        for (const auto label : label_set)
        {
            if (protected_indices.count(label))
            {
                has_protected_vertex = true;
                break;
            }
        }

        if (!has_protected_vertex)
            pv_support_label_sets.push_back(std::move(label_set));
    }
    if (verbose)
    {
        std::cout << "raw pv support cofacets size = " << raw_pv_support_cofacet_indices.size() << '\n';
        std::cout << "pv support label sets with no pv/protected contents size = " << pv_support_label_sets.size() << '\n';
    }

    return pv_support_label_sets;
}

template <typename DistMatType>
double QuotientAndExpand<DistMatType>::computeLabelDistance(
    const size_t i,
    const size_t j,
    const std::vector<std::vector<size_t>>& pv_rep_lists)
{
    const size_t originalvtnum = dist_mat_.getVertexNumber();

    if (i < originalvtnum && j < originalvtnum)
        return dist_mat_.getDistance(i, j);

    double mindist = std::numeric_limits<double>::max();

    if (i < originalvtnum)
    {
        for (const auto vj : pv_rep_lists[j - originalvtnum])
            mindist = std::min(mindist, dist_mat_.getDistance(i, vj));
        return mindist;
    }

    if (j < originalvtnum)
    {
        for (const auto vi : pv_rep_lists[i - originalvtnum])
            mindist = std::min(mindist, dist_mat_.getDistance(vi, j));
        return mindist;
    }

    for (const auto vi : pv_rep_lists[i - originalvtnum])
    {
        for (const auto vj : pv_rep_lists[j - originalvtnum])
            mindist = std::min(mindist, dist_mat_.getDistance(vi, vj));
    }

    return mindist;
}

template <typename DistMatType>
typename QuotientAndExpand<DistMatType>::QuotientEdgeData QuotientAndExpand<DistMatType>::buildQuotientEdges(const std::vector<size_t>& active_labels,
                                                                                                             const std::vector<std::vector<size_t>>& pv_rep_lists,
                                                                                                             const double maxeps, int threadnum)
{
    const int worker_count = threadnum > 0 ? threadnum : 1;
    const size_t originalvtnum = dist_mat_.getVertexNumber();
    omp_set_num_threads(worker_count);

    std::vector<robin_hood::unordered_map<uint64_t, double>> thread_pv_label_distance_hashes(static_cast<size_t>(worker_count));
    std::vector<std::vector<std::pair<int64_t, double>>> thread_edge_workspace(static_cast<size_t>(worker_count));

#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (size_t i = 0; i < active_labels.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_pv_label_distance_hash = thread_pv_label_distance_hashes[threadid];
        auto& thread_edges = thread_edge_workspace[threadid];

        for (size_t j = i + 1; j < active_labels.size(); ++j)
        {
            size_t label_i = active_labels[i];
            size_t label_j = active_labels[j];
            if (label_i > label_j)
                std::swap(label_i, label_j);

            const double weight = computeLabelDistance(label_i, label_j, pv_rep_lists);
            // only hash pv incident edges
            if (label_i >= originalvtnum || label_j >= originalvtnum)
            {
                const uint64_t key = (static_cast<uint64_t>(label_i) << 32) | static_cast<uint64_t>(label_j);
                thread_pv_label_distance_hash.emplace(key, weight);
            }

            if (weight > 0 && weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(this->binomial_table_, label_j, label_i);
                thread_edges.emplace_back(bindex, weight);
            }
        }
    }

    size_t pv_label_pair_count = 0;
    for (const auto& thread_pv_label_distance_hash : thread_pv_label_distance_hashes)
        pv_label_pair_count += thread_pv_label_distance_hash.size();

    robin_hood::unordered_map<uint64_t, double> pv_label_distance_hash;
    pv_label_distance_hash.reserve(pv_label_pair_count);
    for (const auto& thread_pv_label_distance_hash : thread_pv_label_distance_hashes)
    {
        pv_label_distance_hash.insert(thread_pv_label_distance_hash.begin(),
                                      thread_pv_label_distance_hash.end());
    }

    return QuotientEdgeData{std::move(pv_label_distance_hash), SimplexUtility::sortAndMergeSimplexChunks(thread_edge_workspace, threadnum)};
}

template <typename DistMatType>
robin_hood::unordered_map<int64_t, size_t> QuotientAndExpand<DistMatType>::getQuotientActiveEdgeIndexHashTable(const std::vector<std::pair<int64_t, double>>& sorted_quotient_edge, const size_t pvnum)
{
    robin_hood::unordered_map<int64_t, size_t> active_edge_hash_table;
    active_edge_hash_table.reserve(sorted_quotient_edge.size());

    // size_t npts = binomial_table_.size() - 1;
    // original vertices + PV labels
    const size_t npts = dist_mat_.getVertexNumber() + pvnum;

    UnionFind union_find(npts);

    for (size_t i = 0; i < sorted_quotient_edge.size(); ++i)
    {
        const int64_t bindex = sorted_quotient_edge[i].first;

        std::vector<size_t> edge_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, 1);

        auto x = edge_vertices[0];
        auto y = edge_vertices[1];

        if (union_find.unionFind(x) != union_find.unionFind(y))
        {
            union_find.unionSets(x, y);
            continue; // Skip adding to the hash table if they are in different sets
        }

        active_edge_hash_table.emplace(bindex, i);
    }

    return active_edge_hash_table;
}


//legacy QE
template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runQuotient(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> untrimed_pv_index_sets = getPVIndexSets(maxdim, initeps, threadnumber);

    std::vector<std::unordered_set<size_t>> pv_index_sets = trimIndexSets(untrimed_pv_index_sets, initeps);

    /********************************debug*******************************/
    for (auto& pv_vts : pv_index_sets)
    {
        std::cout << "size of the vt = " << pv_vts.size()
                  << "  max pairwise distance = " << getMaxPairwiseDistance(pv_vts)
                  << "  contents :  ";
        for (auto& vt : pv_vts)
            std::cout << vt << "  ";
        std::cout << '\n';
    }

    return pv_index_sets;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::runExpand(const std::vector<std::unordered_set<size_t>>& pv_index_sets, const size_t maxdim, const double maxeps, const int threadnumber)
{
    const size_t originalvtnum = dist_mat_.getVertexNumber(); // number of original vertices
    const size_t pvnum = pv_index_sets.size();
    const size_t total_label_count = originalvtnum + pvnum;

    std::cout << "***********PV num = " << pvnum << "*****************" << '\n';

    SimplexUtility::updateBinomialTable(binomial_table_, originalvtnum, pvnum, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);

    auto active_labels = getActiveLabelIndices(pv_index_sets);
    const auto pv_rep_lists = getPVRepresentativeLists(pv_index_sets);

    auto quotient_edge_data = buildQuotientEdges(active_labels, pv_rep_lists, maxeps, threadnumber);
    auto pv_label_distance_hash = std::move(quotient_edge_data.pv_label_distance_hash);
    auto sorted_quotient_simplex = std::move(quotient_edge_data.sorted_edges);

    auto active_facet_hash = getQuotientActiveEdgeIndexHashTable(sorted_quotient_simplex, pvnum);

    auto sorted_quotient_cofacet = simplex_enumerator.getGeometricCofacetList(
        sorted_quotient_simplex, active_labels, pv_rep_lists, pv_label_distance_hash,
        1, maxeps, threadnumber);

    // Implicit interface graph (no explicit adjacency lists)
    BipartiteGraph bi_graph(1, 1, ImplicitConstructionTag{});

    // Implicit matching object
    MaximumMorseMatching morse_matching;

    std::vector<MaximumMorseMatching::PersistentPairInfo> dim_persistent_pairs;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        // Build the compact cofacet index for implicit adjacency queries.
        auto cofacet_index = SimplexUtility::getCofacetIndex(sorted_quotient_cofacet);

        bi_graph.updateDimensionImplicit(sorted_quotient_cofacet.size(), sorted_quotient_simplex.size());

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_quotient_simplex, sorted_quotient_cofacet,
                                         active_facet_hash, cofacet_index, total_label_count, dim);

        std::cout << "in expand phase (implicit). dim = " << dim
                  << "  cofacet num = " << sorted_quotient_cofacet.size()
                  << "  facet num = " << sorted_quotient_simplex.size() << '\n';

        auto crit_simp_num = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);

        std::cout << "dimensional persistent pairs:" << std::endl;
        if (dim_persistent_pairs.empty())
        {
            std::cout << "  (empty)" << std::endl;
        }
        else
        {
            for (const auto& pair_info : dim_persistent_pairs)
            {
                std::cout << "  (" << pair_info.facet_weight << ", " << pair_info.cofacet_weight << ")" << std::endl;
            }
        }

        dim_persistent_pairs.clear();


        if (dim != maxdim)
        {
            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_facet_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_quotient_cofacet);

            cofacet_index.release();

            // enumerate next cofacet list (geometric PV clique filter)
            sorted_quotient_simplex = simplex_enumerator.getGeometricCofacetList(
                sorted_quotient_cofacet, active_labels, pv_rep_lists, pv_label_distance_hash,
                dim, maxeps, threadnumber);

            std::swap(sorted_quotient_simplex, sorted_quotient_cofacet);
        }
    }

    return;
}

template <typename DistMatType>
std::vector<size_t> QuotientAndExpand<DistMatType>::getActiveLabelIndices(const std::vector<std::unordered_set<size_t>>& pv_index_sets)
{
    const size_t originalvtnum = dist_mat_.getVertexNumber();

    std::vector<bool> is_covered_by_pv(originalvtnum, false);

    for (const auto& pv_indices : pv_index_sets)
    {
        for (const auto vt : pv_indices)
            is_covered_by_pv[vt] = true;
    }

    std::vector<size_t> active_labels;
    active_labels.reserve(originalvtnum);    // number of active original labels is <= original vertex count

    for (size_t i = 0; i < originalvtnum; ++i)
    {
        if (!is_covered_by_pv[i])
            active_labels.push_back(i);
    }

    size_t initpvidx = originalvtnum;
    for (size_t i = 0; i < pv_index_sets.size(); ++i)
    {
        active_labels.push_back(initpvidx);
        initpvidx++;
    }

    // active labels are sorted in ascending order
    return active_labels;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::trimIndexSets(std::vector<std::unordered_set<size_t>>& pv_support_vertex_sets, const double initeps)
{
    // std::sort(gradient_path_vertex_sets.begin(), gradient_path_vertex_sets.end(),
    //           [](const std::unordered_set<size_t>& lhs, const std::unordered_set<size_t>& rhs)
    //           { return lhs.size() > rhs.size(); });

    std::vector<std::unordered_set<size_t>> trimmed_vertex_sets;

    std::unordered_set<size_t> claimed_vertices;

    // const size_t MAX_VERTICES_NUM = MAX_SIZE_ - 1;    // max number of PV representatives. the max value of uint8_t used in EdgeRecord

    for (auto it = pv_support_vertex_sets.rbegin(); it != pv_support_vertex_sets.rend(); ++it)
    {
        auto& vertex_set = *it;
        bool overlap = false;

        for (const auto vertex : vertex_set)
        {
            if (claimed_vertices.count(vertex))
            {
                overlap = true; // discard the vertex set (the smaller set) if it overlaps with claimed vertices
                break;
            }
        }

        if (!overlap)
        {
            if (vertex_set.size() >= MAX_SIZE_) continue;
            if (getMaxPairwiseDistance(vertex_set) > initeps) continue;

            claimed_vertices.insert(vertex_set.begin(), vertex_set.end());
            trimmed_vertex_sets.push_back(std::move(vertex_set));
        }
    }

    return trimmed_vertex_sets;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getPVIndexSets(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> raw_pv_index_sets;

    const size_t originalvtnum = dist_mat_.getVertexNumber(); // number of original vertices

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);

    // 1-simplices (edges) at initeps
    auto sorted_simplex = simplex_enumerator.getSortedVREdges(initeps);

    // active facets for dim=2 are edges not in the MST 
    auto active_facet_hash = SimplexUtility::getActiveEdgeIndexHashTable(binomial_table_, sorted_simplex, originalvtnum);

    // 2-simplices
    auto sorted_cofacet = simplex_enumerator.getSortedVRCofacets(sorted_simplex, 1, initeps, threadnumber);

    // implicit interface graph (no adjacency list)
    BipartiteGraph bi_graph(1, 1, ImplicitConstructionTag{});

    // matching object (implicit)
    MaximumMorseMatching morse_matching;

    // workspace for persistence pairs (optional)
    std::vector<MaximumMorseMatching::PersistentPairInfo> dim_persistent_pairs;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        // Build the compact cofacet index for implicit adjacency queries.
        auto cofacet_index = SimplexUtility::getCofacetIndex(sorted_cofacet);

        bi_graph.updateDimensionImplicit(sorted_cofacet.size(), sorted_simplex.size());

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_simplex, sorted_cofacet,
                                         active_facet_hash, cofacet_index, originalvtnum, dim);

        std::cout << "in quotient phase (implicit), dim = " << dim
                  << "  cofacet num = " << sorted_cofacet.size()
                  << "  facet num = " << sorted_simplex.size() << "\n";

        if (dim != maxdim)
        {
            auto critsimpnum = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);

            std::cout << "dimensional persistent pairs:" << std::endl;
            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& pair_info : dim_persistent_pairs)
                {
                    std::cout << "  (" << pair_info.facet_weight << ", " << pair_info.cofacet_weight << ")" << std::endl;
                }
            }

            dim_persistent_pairs.clear();
        }
        else
        {
            auto match_support_info = morse_matching.implicitMatchAndCollectSupportInfo(matching_context, dim_persistent_pairs, true);
            std::unordered_set<size_t> protected_indices;

            for (const auto protected_facet_list_idx : match_support_info.protected_facet_list_indices)
            {
                const auto facet_bindex = matching_context.sorted_facets[protected_facet_list_idx].first;
                const auto facet_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, facet_bindex, originalvtnum, dim - 1);
                protected_indices.insert(facet_vertices.begin(), facet_vertices.end());
            }

            std::cout << "pv support set num = " << match_support_info.raw_pv_support_cofacet_indices.size() << '\n';

            std::cout << "dimensional persistent pairs:" << std::endl;
            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& pair_info : dim_persistent_pairs)
                {
                    std::cout << "  (" << pair_info.facet_weight << ", " << pair_info.cofacet_weight << ")" << std::endl;
                }
            }

            dim_persistent_pairs.clear();

            // for (auto& support : pv_support_cofacets)
            // {
            //     std::cout << "support cofacet indices: ";
            //     for (auto& idx : support)
            //         std::cout << idx << "  ";
            //     std::cout << '\n';
            // }

            raw_pv_index_sets = getNonMergingPVSupport(matching_context,
                                                       match_support_info.raw_pv_support_cofacet_indices,
                                                       protected_indices,
                                                       originalvtnum,
                                                       true);
        }

        if (dim != maxdim)
        {
            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_facet_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_cofacet);

            cofacet_index.release();

            // enumerate next cofacet list
            sorted_simplex = simplex_enumerator.getSortedVRCofacets(sorted_cofacet, dim, initeps, threadnumber);
            std::swap(sorted_simplex, sorted_cofacet);
        }
    }

    return raw_pv_index_sets;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::reportFalseFacetIdentificationStats(
    const WindowState& win_state,
    const std::vector<std::vector<size_t>>& pv_rep_lists,
    const std::vector<std::pair<int64_t, double>>& facet_list,
    const std::vector<std::pair<int64_t, double>>& cofacet_list,
    const robin_hood::unordered_map<int64_t, uint64_t>& facet_pv_realizations,
    const robin_hood::unordered_map<int64_t, uint64_t>& cofacet_pv_realizations,
    const size_t interface_dim,
    const double eps_lo,
    const double eps_hi,
    const int threadnum)
{
    const size_t origin_vt_num = win_state.original_vertex_number;
    const size_t pv_num = win_state.pv_list.size();

    if (pv_num == 0 || interface_dim < 3 || interface_dim >= MAX_FFI_PACKED_LABELS_)
        return;

    const size_t total_label_count = origin_vt_num + pv_num;
    const size_t cofacet_label_count = interface_dim + 1;
    const size_t facet_label_count = interface_dim;

    const auto facet_hash = SimplexUtility::getSimplexIndexHashTable(facet_list);

    struct FfiStats
    {
        size_t cofacets_with_pv = 0;
        size_t incidence_total = 0;
        size_t facet_all_original = 0;
        size_t missing_facet_realization = 0;
        size_t diff_multi_safe = 0;
        size_t diff_multi_flagged = 0;
        size_t gap_ratio_zero = 0;

        std::vector<float> gap_ratio_samples;

        double worst_overshoot = 0.0;
        double worst_overshoot_cofacet_weight = 0.0;
        double worst_overshoot_facet_weight = 0.0;

        void accumulate(FfiStats& other)
        {
            cofacets_with_pv += other.cofacets_with_pv;
            incidence_total += other.incidence_total;
            facet_all_original += other.facet_all_original;
            missing_facet_realization += other.missing_facet_realization;
            diff_multi_safe += other.diff_multi_safe;
            diff_multi_flagged += other.diff_multi_flagged;
            gap_ratio_zero += other.gap_ratio_zero;
            gap_ratio_samples.insert(gap_ratio_samples.end(),
                                     other.gap_ratio_samples.begin(),
                                     other.gap_ratio_samples.end());
            std::vector<float>().swap(other.gap_ratio_samples);

            if (other.worst_overshoot > worst_overshoot)
            {
                worst_overshoot = other.worst_overshoot;
                worst_overshoot_cofacet_weight = other.worst_overshoot_cofacet_weight;
                worst_overshoot_facet_weight = other.worst_overshoot_facet_weight;
            }
        }
    };

    const int worker_count = threadnum > 0 ? threadnum : 1;
    std::vector<FfiStats> thread_stats(static_cast<size_t>(worker_count));

    struct FfiWorkspace
    {
        std::vector<size_t> cofacet_labels;
        std::vector<size_t> cofacet_witnesses;
        std::vector<size_t> facet_witnesses;
        std::vector<size_t> restricted_witnesses;
        std::vector<size_t> differing_positions;
    };
    std::vector<FfiWorkspace> thread_workspaces(static_cast<size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic, 64) num_threads(worker_count)
    for (size_t ci = 0; ci < cofacet_list.size(); ++ci)
    {
        const int threadid = omp_get_thread_num();
        auto& stats = thread_stats[threadid];
        auto& workspace = thread_workspaces[threadid];

        const int64_t cofacet_bindex = cofacet_list[ci].first;
        const double cofacet_weight = cofacet_list[ci].second;

        const auto realization_iter = cofacet_pv_realizations.find(cofacet_bindex);
        if (realization_iter == cofacet_pv_realizations.end())
            continue;

        ++stats.cofacets_with_pv;
        const uint64_t packed_cofacet_realization = realization_iter->second;

        SimplexUtility::getSimplexVerticesInPlace(binomial_table_, workspace.cofacet_labels,
                                                  cofacet_bindex, total_label_count, interface_dim);

        workspace.cofacet_witnesses.resize(cofacet_label_count);
        for (size_t group = 0; group < cofacet_label_count; ++group)
        {
            const size_t label = workspace.cofacet_labels[group];
            const size_t local_index = (packed_cofacet_realization >> (8 * group)) & 0xFFULL;
            workspace.cofacet_witnesses[group] =
                (label < origin_vt_num) ? label : pv_rep_lists[label - origin_vt_num][local_index];
        }

        int64_t above = 0;
        int64_t below = cofacet_bindex;
        size_t k = interface_dim;

        for (size_t dropped_position = 0; dropped_position < cofacet_label_count; ++dropped_position)
        {
            const size_t dropped_label = workspace.cofacet_labels[dropped_position];
            below -= binomial_table_[dropped_label][k + 1];
            const int64_t facet_bindex = above + below;
            above += binomial_table_[dropped_label][k];
            if (dropped_position + 1 < cofacet_label_count)
                --k;

            ++stats.incidence_total;

            size_t facet_pv_label_count = 0;
            for (size_t q = 0; q < facet_label_count; ++q)
            {
                const size_t label = workspace.cofacet_labels[q < dropped_position ? q : q + 1];
                if (label < origin_vt_num)
                    break;
                ++facet_pv_label_count;
            }

            if (facet_pv_label_count == 0)
            {
                ++stats.facet_all_original;
                continue;
            }

            const auto facet_iter = facet_hash.find(facet_bindex);
            const auto facet_realization_iter = facet_pv_realizations.find(facet_bindex);
            if (facet_iter == facet_hash.end() || facet_realization_iter == facet_pv_realizations.end())
            {
                ++stats.missing_facet_realization;
                continue;
            }

            const double facet_weight = facet_list[facet_iter->second].second;
            const uint64_t packed_facet_realization = facet_realization_iter->second;

            workspace.facet_witnesses.resize(facet_label_count);
            workspace.restricted_witnesses.resize(facet_label_count);
            workspace.differing_positions.clear();

            for (size_t q = 0; q < facet_label_count; ++q)
            {
                const size_t cofacet_position = (q < dropped_position) ? q : q + 1;
                const size_t label = workspace.cofacet_labels[cofacet_position];
                const size_t facet_local_index = (packed_facet_realization >> (8 * q)) & 0xFFULL;

                workspace.facet_witnesses[q] =
                    (label < origin_vt_num) ? label : pv_rep_lists[label - origin_vt_num][facet_local_index];
                workspace.restricted_witnesses[q] = workspace.cofacet_witnesses[cofacet_position];

                if (workspace.facet_witnesses[q] != workspace.restricted_witnesses[q])
                    workspace.differing_positions.push_back(q);
            }

            double restricted_weight = 0.0;
            for (size_t a = 0; a + 1 < facet_label_count; ++a)
            {
                for (size_t b = a + 1; b < facet_label_count; ++b)
                {
                    restricted_weight = std::max(
                        restricted_weight,
                        dist_mat_.getDistance(workspace.restricted_witnesses[a],
                                              workspace.restricted_witnesses[b]));
                }
            }

            const double realization_gap = restricted_weight - facet_weight;
            const double realization_gap_ratio = realization_gap / facet_weight;
            stats.gap_ratio_samples.push_back(static_cast<float>(realization_gap_ratio));

            if (realization_gap_ratio == 0.0)
                ++stats.gap_ratio_zero;

            if (workspace.differing_positions.size() > 1)
            {
                // A clique on Y|F union Z is a one-step filler. Failure of this local
                // sufficient test is only a flag; a filler chain may still exist.
                double max_cross_distance = 0.0;
                for (const auto q : workspace.differing_positions)
                {
                    for (size_t r = 0; r < facet_label_count; ++r)
                    {
                        max_cross_distance = std::max(
                            max_cross_distance,
                            dist_mat_.getDistance(workspace.restricted_witnesses[q],
                                                  workspace.facet_witnesses[r]));
                    }
                }

                if (max_cross_distance <= cofacet_weight)
                {
                    ++stats.diff_multi_safe;
                }
                else
                {
                    ++stats.diff_multi_flagged;
                    const double overshoot = max_cross_distance - cofacet_weight;
                    if (overshoot > stats.worst_overshoot)
                    {
                        stats.worst_overshoot = overshoot;
                        stats.worst_overshoot_cofacet_weight = cofacet_weight;
                        stats.worst_overshoot_facet_weight = facet_weight;
                    }
                }
            }
        }
    }

    FfiStats total;
    for (auto& stats : thread_stats)
        total.accumulate(stats);

    std::sort(total.gap_ratio_samples.begin(), total.gap_ratio_samples.end());
    const auto quantile = [&](const double q) -> double
    {
        if (total.gap_ratio_samples.empty())
            return 0.0;
        const size_t position = std::min(
            total.gap_ratio_samples.size() - 1,
            static_cast<size_t>(q * (total.gap_ratio_samples.size() - 1) + 0.5));
        return total.gap_ratio_samples[position];
    };

    const size_t pv_incident = total.incidence_total - total.facet_all_original;
    const size_t measured_pv_incident = total.gap_ratio_samples.size();
    const size_t multi_incidence = total.diff_multi_safe + total.diff_multi_flagged;
    const auto percentage = [&](const size_t count) -> double
    {
        if (measured_pv_incident == 0)
            return 0.0;
        return 100.0 * count / measured_pv_incident;
    };

    if (pv_incident > 0)
    {
        if (measured_pv_incident > 0)
        {
            std::cout << std::endl << "FFI stats \n";
            std::cout << "  [ffi2] realization gap ratio (w(Y|F) - w(F)) / w(F): zero = "
                      << percentage(total.gap_ratio_zero) << "%"
                      << "  p50 = " << quantile(0.50)
                      << "  p90 = " << quantile(0.90)
                      << "  p99 = " << quantile(0.99) << '\n';

            std::cout << "  [ffi2] differing PV coordinates:"
                      << "  multi indcidence = " << percentage(multi_incidence) << "%"
                      << "  clique-certified = " << percentage(total.diff_multi_safe) << "%"
                      << "  void capable = "
                      << percentage(total.diff_multi_flagged) << "%\n";
        }

        if (total.diff_multi_flagged > 0)
        {
            std::cout << "  [ffi2] worst no-clique-certificate incidence: carrier overshoot = " << total.worst_overshoot
                      << "  w(facet) = " << total.worst_overshoot_facet_weight
                      << "  w(cofacet) = " << total.worst_overshoot_cofacet_weight << '\n';
        }

        if (total.missing_facet_realization > 0)
        {
            std::cout << "  [ffi2] WARNING: missing facet realizations = "
                      << total.missing_facet_realization << '\n';
        }

        std::cout << "  [ffi2-line]"
                  << " eps_lo=" << eps_lo
                  << " eps_hi=" << eps_hi
                  << " dim=" << interface_dim
                  << " pv_cofacets =" << total.cofacets_with_pv
                //   << " inc=" << total.incidence_total
                //   << " all_orig=" << total.facet_all_original
                //   << " gap_ratio_zero_pct=" << percentage(total.gap_ratio_zero) << "%"
                //   << " gap_ratio_p50=" << quantile(0.50)
                //   << " gap_ratio_p90=" << quantile(0.90)
                //   << " gap_ratio_p99=" << quantile(0.99)
                  << " diffm_safe_pct=" << percentage(total.diff_multi_safe) << "%"
                  << " diffm_flagged_pct=" << percentage(total.diff_multi_flagged) << "%"
                //   << " measured=" << measured_pv_incident
                //   << " missing=" << total.missing_facet_realization
                  <<"\nEnd of FFI stats"
                  << std::endl;

    }
}

template class QuotientAndExpand<NormalDistMat>;
// template class QuotientAndExpand<SparseDistMat>;
