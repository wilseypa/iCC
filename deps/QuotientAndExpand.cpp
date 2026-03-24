#include <unordered_map>
#include <stdexcept>

#include "omp.h"

#include "SimplexUtility.hpp"
#include "SimplexEnumerator.hpp"
#include "MaximumMorseMatching.hpp"
#include "QuotientAndExpand.hpp"

template class QuotientAndExpand<NormalDistMat>;
// template class QuotientAndExpand<SparseDistMat>;

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::runPiecewisePH(const std::vector<double>& eps_breaks, const size_t maxdim, const int threadnumber)
{
    WindowState win_state(dist_mat_.getVertexNumber());

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

        auto raw_pv_label_sets = runWindow(win_state, maxdim, eps_lo, eps_hi, threadnumber, collect_pv);

        if (!collect_pv) break;

        auto new_pv_list = trimPVCandidates(win_state, raw_pv_label_sets, eps_hi);

        rebuildWindowState(win_state, new_pv_list);

        /********************************debug*******************************/
        const auto getMaxPairwiseDistance = [this](const std::unordered_set<size_t>& index_set)
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
                    maxdist = std::max(maxdist, dist_mat_.getDistance(*first, *second));
                }
            }

            return maxdist;
        };

        std::cout << "after eps range "<<eps_lo<< "  " <<eps_hi<< "  new pv number = "<<new_pv_list.size()
                  <<"  total pv number = "<<win_state.pv_flat_index_set_list.size()<<std::endl;

        std::cout << "pv flat index sets:" << std::endl;
        if (win_state.pv_flat_index_set_list.empty())
        {
            std::cout << "  (empty)" << std::endl;
        }
        else
        {
            size_t pv_idx = 0;
            for (const auto& pv_index_set : win_state.pv_flat_index_set_list)
            {
                std::cout << "  [" << pv_idx++ << "] ";
                bool first = true;
                for (const auto& flat_index : pv_index_set)
                {
                    if (!first)
                    {
                        std::cout << " ";
                    }
                    std::cout << flat_index;
                    first = false;
                }
                
                std::cout <<"    diameter = "<<getMaxPairwiseDistance(pv_index_set)<<'\n';
            }
            std::cout<<std::endl;
        }


    }

    return;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runQuotient(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> untrimed_pv_index_sets = getPVIndexSets(maxdim, initeps, threadnumber);

    std::vector<std::unordered_set<size_t>> pv_index_sets = trimIndexSets(untrimed_pv_index_sets, initeps);

    /********************************debug*******************************/
    const auto getMaxPairwiseDistance = [this](const std::unordered_set<size_t>& index_set)
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
                maxdist = std::max(maxdist, dist_mat_.getDistance(*first, *second));
            }
        }

        return maxdist;
    };

    for (auto& pv_vts : pv_index_sets)
    {
        std::cout << "size of the vt = " << pv_vts.size()
                  << "  max pairwise distance = " << getMaxPairwiseDistance(pv_vts)
                  << "  contents :  ";
        for (auto &vt : pv_vts)
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
    const size_t npts = originalvtnum + pvnum; // original + virtual vertices

    std::cout << "***********virtual vertex num = " << pvnum << "*****************" << '\n';

    SimplexUtility::updateBinomialTable(binomial_table_, originalvtnum, pvnum, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);

    auto active_vertices = getActiveVertexIndices(pv_index_sets);

    auto virtual_distance_hash_table = getVirtualDistanceHashTable(active_vertices, pv_index_sets, threadnumber);

    auto sorted_virtual_simplex = getSortedVirtualEdgeList(active_vertices, virtual_distance_hash_table, maxeps, threadnumber);

    auto active_facet_hash = getVirtualActiveEdgeIndexHashTable(sorted_virtual_simplex, pvnum);

    // auto sorted_virtual_cofacet = getVirtualCofacetList(sorted_virtual_simplex, active_vertices, virtual_distance_hash_table, 1, maxeps, threadnumber);
    auto sorted_virtual_cofacet = simplex_enumerator.getGeometricVirtualCofacetList(sorted_virtual_simplex, active_vertices, pv_index_sets, 1, maxeps, threadnumber);

    // Implicit interface graph (no explicit adjacency lists)
    BipartiteGraph bi_graph(1, 1, ImplicitConstructionTag{});

    // Implicit matching object
    MaximumMorseMatching morse_matching;

    std::vector<std::pair<double, double>> dim_persistent_pairs;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        // Build cofacet hash for implicit adjacency queries
        auto cofacet_hash = SimplexUtility::getSimplexIndexHashTable(sorted_virtual_cofacet);

        bi_graph.updateDimensionImplicit(sorted_virtual_cofacet.size(), sorted_virtual_simplex.size());

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_virtual_simplex, sorted_virtual_cofacet,
                                         active_facet_hash, cofacet_hash, npts, dim);

        std::cout << "in expand phase (implicit). dim = " << dim
                  << "  cofacet num = " << sorted_virtual_cofacet.size()
                  << "  facet num = " << sorted_virtual_simplex.size() << '\n';

        auto crit_simp_num = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);

        std::cout << "dimensional persistent pairs:" << std::endl;
        if (dim_persistent_pairs.empty())
        {
            std::cout << "  (empty)" << std::endl;
        }
        else
        {
            for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
            {
                std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
            }
        }

        dim_persistent_pairs.clear();


        if (dim != maxdim)
        {
            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_facet_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_virtual_cofacet);

            // enumerate next cofacet list (geometric PV clique filter)
            sorted_virtual_simplex = simplex_enumerator.getGeometricVirtualCofacetList(sorted_virtual_cofacet, active_vertices, pv_index_sets, dim, maxeps, threadnumber);

            std::swap(sorted_virtual_simplex, sorted_virtual_cofacet);
        }
    }

    return;
}

/* private helpers */

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runWindow(const WindowState& win_state, const size_t maxdim, const double eps_lo, const double eps_hi,
                                                                                  const int threadnumber, const bool collect_pv)
{
    const size_t original_vt_num = win_state.original_vertex_number;
    const size_t pv_num = win_state.pv_flat_index_set_list.size();
    const size_t npts = original_vt_num + pv_num;

    SimplexUtility::updateBinomialTable(binomial_table_, original_vt_num, pv_num, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);

    auto& pv_index_sets = win_state.pv_flat_index_set_list;
    auto& active_labels = win_state.active_label_list;

    auto virtual_distance_hash = getVirtualDistanceHashTable(active_labels, pv_index_sets, threadnumber);

    auto sorted_virtual_simplex = getSortedVirtualEdgeList(active_labels, virtual_distance_hash, eps_hi, threadnumber);

    auto active_simplex_hash = getVirtualActiveEdgeIndexHashTable(sorted_virtual_simplex, pv_num);

    auto sorted_virtual_cofacet = simplex_enumerator.getGeometricVirtualCofacetList(sorted_virtual_simplex, 
                                                                                    active_labels, pv_index_sets, 1, eps_hi, threadnumber);
    
    BipartiteGraph bi_graph(1, 1, ImplicitConstructionTag{});

    std::vector<std::pair<double, double>> dim_persistent_pairs;

    std::vector<std::unordered_set<size_t>> raw_pv_label_sets;

    MaximumMorseMatching morse_matching;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        bi_graph.updateDimensionImplicit(sorted_virtual_cofacet.size(), sorted_virtual_simplex.size());

        auto cofacet_hash = SimplexUtility::getSimplexIndexHashTable(sorted_virtual_cofacet);

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_virtual_simplex, sorted_virtual_cofacet,
                                         active_simplex_hash, cofacet_hash, npts, dim);


        if (dim != maxdim)
        {
            auto crit_simp_num = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);

            std::cout << "in eps range "<<eps_lo<< "  " <<eps_hi<< "    dimension = " <<dim
                      << "  cofacet num = " << sorted_virtual_cofacet.size()
                      << "  facet num = " << sorted_virtual_simplex.size() <<'\n'
                      << "   persistent pairs:" << std::endl;


            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
                {
                    std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
                }
            }

            dim_persistent_pairs.clear();

            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_simplex_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_virtual_cofacet);

            // enumerate next cofacet list
            sorted_virtual_simplex = simplex_enumerator.getGeometricVirtualCofacetList(sorted_virtual_cofacet, active_labels, pv_index_sets, dim, eps_hi, threadnumber);
            std::swap(sorted_virtual_simplex, sorted_virtual_cofacet);
        }
        else if (collect_pv)
        {
            auto pv_support_cofacets = morse_matching.implicitMatchAndCollectPVSupports(matching_context, dim_persistent_pairs);

            std::cout << "in eps range "<<eps_lo<< "  " <<eps_hi<< "    dimension = " <<dim
                      << "  cofacet num = " << sorted_virtual_cofacet.size()
                      << "  facet num = " << sorted_virtual_simplex.size() <<'\n'
                      << "   persistent pairs:" << std::endl;

            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
                {
                    std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
                }
            }

            raw_pv_label_sets = getPVSupportLabelSets(matching_context, pv_support_cofacets, npts, dim);
        }
        else
        {
            auto crit_simp_num = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);

            std::cout << "in eps range "<<eps_lo<< "  " <<eps_hi<< "    dimension = " <<dim
                      << "  cofacet num = " << sorted_virtual_cofacet.size()
                      << "  facet num = " << sorted_virtual_simplex.size() <<'\n'
                      << "   persistent pairs:" << std::endl;

            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
                {
                    std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
                }
            }

            dim_persistent_pairs.clear();
        }
    }

    return raw_pv_label_sets;
}

template <typename DistMatType>
std::vector<typename QuotientAndExpand<DistMatType>::PVCandidate> 
QuotientAndExpand<DistMatType>::trimPVCandidates(const WindowState& win_state, const std::vector<std::unordered_set<size_t>>& raw_label_sets, const double eps_hi)
{
    const auto getMaxPairwiseDistance = [this](const std::unordered_set<size_t>& index_set)
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
                maxdist = std::max(maxdist, dist_mat_.getDistance(*first, *second));
            }
        }

        return maxdist;
    };

    std::vector<PVCandidate> accepted_pv_list;

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

            if (flat_index_set.size() >= MAX_SIZE_) continue;

            if (getMaxPairwiseDistance(flat_index_set) > eps_hi) continue;

            claimed_labels.insert(label_set.begin(), label_set.end());

            accepted_pv_list.push_back(PVCandidate{label_set, std::move(flat_index_set)});
        }
    }

    return accepted_pv_list;
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
            const auto& pv_indices = win_state.pv_flat_index_set_list[label - origin_vt_num];
            flat_index_set.insert(pv_indices.begin(), pv_indices.end());
        }
    }

    return flat_index_set;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::rebuildWindowState(WindowState& win_state, const std::vector<PVCandidate>& new_pv_list)
{
    const size_t origin_vt_num = win_state.original_vertex_number;
    const size_t old_pv_num = win_state.pv_flat_index_set_list.size();

    auto old_pv_index_sets = std::move(win_state.pv_flat_index_set_list);

    std::vector<bool> absorbed_old_pv(old_pv_num, false);

    for (const auto& new_pv : new_pv_list)
    {
        for (const auto label : new_pv.label_set)
        {
            if (label >= origin_vt_num)
                absorbed_old_pv[label - origin_vt_num] = true;
        }
    }

    std::vector<std::unordered_set<size_t>> updated_pv_index_sets;

    // old unmerged PVs first
    for (size_t i = 0; i < old_pv_num; ++i)
    {
        if (!absorbed_old_pv[i])
            updated_pv_index_sets.push_back(std::move(old_pv_index_sets[i]));
    }

    // then new PVs
    for (const auto& new_pv : new_pv_list)
    {
        updated_pv_index_sets.push_back(new_pv.flat_index_set);
    }

    // rebuild active original labels from final coverage
    std::vector<bool> is_covered(origin_vt_num, false);
    for (const auto& pv_indices : updated_pv_index_sets)
    {
        for (const auto idx : pv_indices)
            is_covered[idx] = true;
    }

    std::vector<size_t> active_label_list;
    active_label_list.reserve(origin_vt_num + updated_pv_index_sets.size());

    for (size_t i = 0; i < origin_vt_num; ++i)
    {
        if (!is_covered[i]) active_label_list.push_back(i);
    }

    for (size_t i = 0; i < updated_pv_index_sets.size(); ++i)
    {
        active_label_list.push_back(origin_vt_num + i);
    }

    win_state.active_label_list = std::move(active_label_list);
    win_state.pv_flat_index_set_list = std::move(updated_pv_index_sets);

    return;
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
    std::vector<std::pair<double, double>> dim_persistent_pairs;

    for (size_t dim = 2; dim <= maxdim; ++dim)
    {
        // build cofacet hash for implicit adjacency queries
        auto cofacet_hash = SimplexUtility::getSimplexIndexHashTable(sorted_cofacet);

        bi_graph.updateDimensionImplicit(sorted_cofacet.size(), sorted_simplex.size());

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_simplex, sorted_cofacet,
                                         active_facet_hash, cofacet_hash, originalvtnum, dim);

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
                for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
                {
                    std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
                }
            }

            dim_persistent_pairs.clear();
        }
        else
        {
            auto pv_support_cofacets = morse_matching.implicitMatchAndCollectPVSupports(matching_context, dim_persistent_pairs);

            std::cout << "pv support set num = " << pv_support_cofacets.size() << '\n';

            std::cout << "dimensional persistent pairs:" << std::endl;
            if (dim_persistent_pairs.empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& [facetweight, cofacetweight] : dim_persistent_pairs)
                {
                    std::cout << "  (" << facetweight << ", " << cofacetweight << ")" << std::endl;
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

            raw_pv_index_sets = getPVSupportLabelSets(matching_context, pv_support_cofacets, originalvtnum, dim);
        }

        if (dim != maxdim)
        {
            // facets for the next dimension are the unmatched cofacets from the current dimension
            active_facet_hash = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_cofacet);

            // enumerate next cofacet list
            sorted_simplex = simplex_enumerator.getSortedVRCofacets(sorted_cofacet, dim, initeps, threadnumber);
            std::swap(sorted_simplex, sorted_cofacet);
        }
    }

    return raw_pv_index_sets;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getPVSupportLabelSets(const MatchingContext &matching_context, const std::vector<std::vector<size_t>>& pv_support_cofacets,
                                                                                              const size_t npts, const size_t dim)
{
    std::vector<std::unordered_set<size_t>> pv_support_vertex_sets;
    pv_support_vertex_sets.reserve(pv_support_cofacets.size());

    // auto npt = matching_context.binomial_table.size() - 1;

    for (const auto& pv_support : pv_support_cofacets)
    {
        //need set for trimming
        std::unordered_set<size_t> vertex_set;

        for (auto cofacetidx : pv_support)
        {
            auto bindex = matching_context.sorted_cofacets[cofacetidx].first;
            auto simplex_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, bindex, npts, dim);
            vertex_set.insert(simplex_vertices.begin(), simplex_vertices.end());
        }
        if (vertex_set.size() < MAX_SIZE_)
            pv_support_vertex_sets.push_back(std::move(vertex_set));
    }

    return pv_support_vertex_sets;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::trimIndexSets(std::vector<std::unordered_set<size_t>>& pv_support_vertex_sets, const double initeps)
{
    // std::sort(gradient_path_vertex_sets.begin(), gradient_path_vertex_sets.end(),
    //           [](const std::unordered_set<size_t> &lhs, const std::unordered_set<size_t> &rhs)
    //           { return lhs.size() > rhs.size(); });

    const auto getMaxPairwiseDistance = [this](const std::unordered_set<size_t>& index_set)
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
                maxdist = std::max(maxdist, dist_mat_.getDistance(*first, *second));
            }
        }

        return maxdist;
    };

    std::vector<std::unordered_set<size_t>> trimmed_vertex_sets;

    std::unordered_set<size_t> claimed_vertices;

    // const size_t MAX_VERTICES_NUM = MAX_SIZE_ - 1;    //max number of virtual vertex. the max value of uint8_t used in EdgeRecord

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
            if (getMaxPairwiseDistance(vertex_set) > initeps)
                continue;

            claimed_vertices.insert(vertex_set.begin(), vertex_set.end());
            trimmed_vertex_sets.push_back(std::move(vertex_set));
        }
    }

    return trimmed_vertex_sets;
}

template <typename DistMatType>
std::vector<size_t> QuotientAndExpand<DistMatType>::getActiveVertexIndices(const std::vector<std::unordered_set<size_t>>& pv_index_sets)
{
    const size_t originalvtnum = dist_mat_.getVertexNumber();

    std::vector<bool> is_virtualized(originalvtnum, false);

    for (const auto &pv_indices : pv_index_sets)
    {
        for (const auto vt : pv_indices)
            is_virtualized[vt] = true;
    }

    std::vector<size_t> active_vertices;
    active_vertices.reserve(originalvtnum);    //num of active vertices <= original

    for (size_t i = 0; i < originalvtnum; ++i)
    {
        if (!is_virtualized[i])
            active_vertices.push_back(i);
    }

    size_t initpvidx = originalvtnum;
    for (size_t i = 0; i < pv_index_sets.size(); ++i)
    {
        // if (virtual_vt_set.empty()) continue;
        active_vertices.push_back(initpvidx);
        initpvidx++;
    }

    // active vertices indices are sorted in ascending order
    return active_vertices;
}

template <typename DistMatType>
double QuotientAndExpand<DistMatType>::computeVirtualDistance(const size_t i, const size_t j, const std::vector<std::unordered_set<size_t>> &pv_index_sets)
{
    const size_t originalvtnum = dist_mat_.getVertexNumber();

    // passed in as i < j

    if (i < originalvtnum && j < originalvtnum)
        return dist_mat_.getDistance(i, j);

    const std::unordered_set<size_t> temp_i_set = (i < originalvtnum) ? std::unordered_set<size_t>{i} : std::unordered_set<size_t>();
    const std::unordered_set<size_t> temp_j_set = (j < originalvtnum) ? std::unordered_set<size_t>{j} : std::unordered_set<size_t>();

    const std::unordered_set<size_t> &vtset_i = (i < originalvtnum) ? temp_i_set : pv_index_sets[i - originalvtnum];
    const std::unordered_set<size_t> &vtset_j = (j < originalvtnum) ? temp_j_set : pv_index_sets[j - originalvtnum];

    double mindist = std::numeric_limits<double>::max();
    for (const auto &vi : vtset_i)
    {
        for (const auto &vj : vtset_j)
        {
            mindist = std::min(mindist, dist_mat_.getDistance(vi, vj));
        }
    }

    return mindist;
}

template <typename DistMatType>
robin_hood::unordered_map<uint64_t, double> QuotientAndExpand<DistMatType>::getVirtualDistanceHashTable(const std::vector<size_t> &active_vertices, const std::vector<std::unordered_set<size_t>> &pv_index_sets, int threadnum)
{
    omp_set_num_threads(threadnum);

    std::vector<robin_hood::unordered_map<uint64_t, double>> thread_hash_tables(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < active_vertices.size(); ++i)
    {
        size_t threadid = omp_get_thread_num();
        auto &thread_hash_table = thread_hash_tables[threadid];

        // index < std::numeric_limits<uint32_t>::max()
        for (size_t j = i + 1; j < active_vertices.size(); ++j)
        {
            uint64_t key = (static_cast<uint64_t>(active_vertices[i]) << 32) | static_cast<uint64_t>(active_vertices[j]);
            double dist = computeVirtualDistance(active_vertices[i], active_vertices[j], pv_index_sets);
            thread_hash_table.emplace(key, dist);
        }
    }

    robin_hood::unordered_map<uint64_t, double> virtual_distance_hash_table;
    for (const auto &thread_hash_table : thread_hash_tables)
    {
        virtual_distance_hash_table.insert(thread_hash_table.begin(), thread_hash_table.end());
    }

    return virtual_distance_hash_table;
}

template <typename DistMatType>
double QuotientAndExpand<DistMatType>::getVirtualDistance(size_t i, size_t j, const robin_hood::unordered_map<uint64_t, double> &virtual_distance_hash_table)
{
    if (i > j)
        std::swap(i, j);

    uint64_t key = (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);

    auto it = virtual_distance_hash_table.find(key);
    if (it != virtual_distance_hash_table.end())
    {
        return it->second;
    }
    else
    {
        return -1.0; // Return -1.0 if the distance is not found in the hash table
    }
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> QuotientAndExpand<DistMatType>::getSortedVirtualEdgeList(const std::vector<size_t> &active_vertices,
                                                                                             const robin_hood::unordered_map<uint64_t, double> &virtual_distance_hash_table, const double maxeps, int threadnum)
{
    omp_set_num_threads(threadnum);

    std::vector<std::vector<std::pair<int64_t, double>>> thread_edge_workspace(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < active_vertices.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto &thread_edges = thread_edge_workspace[threadid];

        for (size_t j = i + 1; j < active_vertices.size(); ++j)
        {
            double weight = getVirtualDistance(active_vertices[i], active_vertices[j], virtual_distance_hash_table);
            if (weight > 0 && weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(this->binomial_table_, active_vertices[j], active_vertices[i]);
                thread_edges.emplace_back(bindex, weight);
            }
        }
    }

    std::vector<std::pair<int64_t, double>> sorted_edges;
    for (const auto &thread_edges : thread_edge_workspace)
    {
        sorted_edges.insert(sorted_edges.end(), thread_edges.begin(), thread_edges.end());
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_edges);

    return sorted_edges;
}

template <typename DistMatType>
robin_hood::unordered_map<int64_t, size_t> QuotientAndExpand<DistMatType>::getVirtualActiveEdgeIndexHashTable(const std::vector<std::pair<int64_t, double>> &sorted_virtual_edge, const size_t pvnum)
{
    robin_hood::unordered_map<int64_t, size_t> active_edge_hash_table;
    active_edge_hash_table.reserve(sorted_virtual_edge.size());

    // size_t npts = binomial_table_.size() - 1;
    // original number of vertices + number of virtual vertices
    const size_t npts = dist_mat_.getVertexNumber() + pvnum;

    UnionFind union_find(npts);

    for (size_t i = 0; i < sorted_virtual_edge.size(); ++i)
    {
        const int64_t bindex = sorted_virtual_edge[i].first;

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


