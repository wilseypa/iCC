#include <algorithm>

#include "omp.h"

#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"

template class SimplexEnumerator<NormalDistMat>;
// template class SimplexEnumerator<SparseDistMat>;

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVREdges(const double maxeps)
{
    //********************use openmp later************************//

    std::vector<std::pair<int64_t, double>> sorted_edge;

    size_t npt = dist_mat_.getVertexNumber();

    for (size_t i = 0; i < npt - 1; i++)
    {
        for (size_t j = i + 1; j < npt; j++)
        {
            double weight = dist_mat_.getDistance(i, j);
            if (weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(binomial_table_, j, i);
                sorted_edge.push_back(std::make_pair(bindex, weight));
            }
        }
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_edge);

    return sorted_edge;
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex_list, const size_t dim, const double maxeps, const int threadnum)
{
    // dim == simplex dimension == cofacet dimension - 1
    size_t npts = binomial_table_.size() - 1;

    std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(threadnum);

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_simplex_list.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];

        const int64_t bindex = sorted_simplex_list[i].first;
        const double weight = sorted_simplex_list[i].second;

        std::vector<size_t> simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        const size_t minvt = simplex_vertices.back(); // sorted in descending order
        for (size_t covt = 0; covt < minvt; ++covt)
        {
            double newweight = 0.0;
            for (const auto& vt : simplex_vertices)
            {
                newweight = std::max(newweight, dist_mat_.getDistance(covt, vt));
            }

            double cofacetweight = std::max(newweight, weight);

            if (newweight < maxeps)
            {
                int64_t shiftedbindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vertices, 1);
                int64_t cofacetbindex = shiftedbindex + covt;
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    return SimplexUtility::sortAndMergeSimplexChunks(thread_workspace, threadnum);
}

#ifdef BUILD_ALPHA_COMPLEX
template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedAlphaCells(const std::vector<std::vector<int64_t>>& binomial_table,
                                                                                            std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t>& vertex_handle_index,
                                                                                            CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>& delaunay_d, const size_t dim, double maxeps)
{
    std::vector<std::pair<int64_t, double>> sortd_d_cell;

    if (dim == 0)
        return sortd_d_cell;

    int maxdim = delaunay_d.maximal_dimension();

    std::vector<size_t> simplex_pt;
    simplex_pt.reserve(dim + 1);

    if (dim == maxdim - 1)
    {
        for (auto facetit = delaunay_d.finite_facets_begin(); facetit != delaunay_d.finite_facets_end(); facetit++)
        {
            auto fullcellit = facetit->full_cell();
            auto covt = facetit->index_of_covertex();

            auto neighborit = fullcellit->neighbor(covt);
            if (neighborit < fullcellit)
                continue; // skip to avoid double counting the same facet

            for (size_t i = 0; i <= maxdim; i++)
            {
                if (i == covt)
                    continue;                         // skip the co_vertex of the facet
                auto vhandle = fullcellit->vertex(i); // need to use the full cell iter. facet iter does not have vertex method
                simplex_pt.push_back(vertex_handle_index[vhandle]);
            }
            std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());
            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, simplex_pt, 0);
            double weight = getAlphaSimplexWeight(simplex_pt);
            if (weight < maxeps)
                sortd_d_cell.emplace_back(bindex, weight);
            simplex_pt.clear();
        }
        SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
        return sortd_d_cell;
    }

    if (dim == maxdim)
    {
        for (auto fullcellit = delaunay_d.finite_full_cells_begin(); fullcellit != delaunay_d.finite_full_cells_end(); fullcellit++)
        {
            for (size_t i = 0; i <= dim; i++)
            {
                auto vhandle = fullcellit->vertex(i);
                simplex_pt.push_back(vertex_handle_index[vhandle]);
            }
            std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());
            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, simplex_pt, 0);
            double weight = getAlphaSimplexWeight(simplex_pt);
            if (weight < maxeps)
                sortd_d_cell.emplace_back(bindex, weight);
            simplex_pt.clear();
        }
        SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
        return sortd_d_cell;
    }

    // for the rest of the dim
    std::vector<size_t> cell_pt;
    cell_pt.reserve(dim + 1);
    std::unordered_set<int64_t> cell_bindex_lookup;
    for (auto fullcellit = delaunay_d.finite_full_cells_begin(); fullcellit != delaunay_d.finite_full_cells_end(); fullcellit++)
    {
        for (size_t i = 0; i <= maxdim; i++)
        {
            auto vhandle = fullcellit->vertex(i);
            simplex_pt.push_back(vertex_handle_index[vhandle]);
        }
        std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());

        // subroutine from lhf getDimEdges(int dim)
        size_t powsize = pow(2, maxdim + 1);
        for (size_t counter = 1; counter < powsize; counter++)
        {
            // count the number of 1 in binary form of counter
            if (__builtin_popcount(counter) != dim + 1)
                continue;

            // collect the corresponding pt(subset) of full cell
            for (size_t i = 0; i < maxdim + 1; i++)
            {
                if (counter & (1 << i))
                    cell_pt.push_back(simplex_pt[i]);
            }

            std::sort(cell_pt.begin(), cell_pt.end(), std::greater<size_t>());

            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, cell_pt, 0);

            // if (bindex == 0)
            // {
            //     std::cout<<"zero bindex edge = "<<cell_pt[0]<<"  "<<cell_pt[1]<<'\n';
            // }

            // already found
            if (cell_bindex_lookup.find(bindex) != cell_bindex_lookup.end())
            {
                cell_pt.clear();
                continue;
            }

            cell_bindex_lookup.insert(bindex);
            double weight = getAlphaSimplexWeight(cell_pt);
            if (weight < maxeps)
                sortd_d_cell.emplace_back(bindex, weight);
            cell_pt.clear(); // clean up the pt array for d cell
        }
        simplex_pt.clear(); // clean up the pt array for full cell
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
    return sortd_d_cell;
}

template <typename DistMatType>
double SimplexEnumerator<DistMatType>::getAlphaSimplexWeight(const std::vector<size_t>& alpha_simplex)
{
    double weight = 0;
    // simplex pt is sorted in descending order
    for (auto rfirst = alpha_simplex.rbegin(); rfirst != alpha_simplex.rend() - 1; rfirst++)
    {
        for (auto rsecond = rfirst + 1; rsecond != alpha_simplex.rend(); rsecond++)
        {
            weight = std::max(weight, dist_mat_.getDistance(*rfirst, *rsecond));
        }
    }
    return weight;
}
#endif

// ======================= Geometric virtual enumeration (pseudo-vertex) =======================

template <typename DistMatType>
void SimplexEnumerator<DistMatType>::prepareFacetWitnessContext(WitnessWorkspace& ws,
                                                                const std::vector<size_t>& facet_labels,
                                                                const std::vector<std::vector<size_t>>& pv_rep_lists,
                                                                const size_t originalvtnum, const double maxeps) const
{
    const size_t K = facet_labels.size();

    if (ws.singleton_slots.size() < K + 1)
        ws.singleton_slots.resize(K + 1);
    if (ws.rep_ptrs.size() < K + 1)
        ws.rep_ptrs.resize(K + 1);

    for (size_t i = 0; i < K; ++i)
    {
        const size_t label = facet_labels[i];

        if (label < originalvtnum)
        {
            ws.singleton_slots[i].assign(1, label);
            ws.rep_ptrs[i] = &ws.singleton_slots[i];
        }
        else
        {
            ws.rep_ptrs[i] = &pv_rep_lists[label - originalvtnum];
        }
    }

    ws.facet_edges.clear();

    for (size_t i = 0; i + 1 < K; ++i)
    {
        const auto& rep_i = *ws.rep_ptrs[i];

        for (size_t j = i + 1; j < K; ++j)
        {
            const auto& rep_j = *ws.rep_ptrs[j];

            for (size_t a = 0; a < rep_i.size(); ++a)
            {
                for (size_t b = 0; b < rep_j.size(); ++b)
                {
                    const double w = dist_mat_.getDistance(rep_i[a], rep_j[b]);
                    if (w < maxeps)
                    {
                        ws.facet_edges.push_back({w, static_cast<uint8_t>(i), static_cast<uint8_t>(j),
                                                  static_cast<uint8_t>(a), static_cast<uint8_t>(b)});
                    }
                }
            }
        }
    }

    std::sort(ws.facet_edges.begin(), ws.facet_edges.end());
}

template <typename DistMatType>
void SimplexEnumerator<DistMatType>::prepareCovtWitnessGroup(WitnessWorkspace& ws, const size_t covt, const size_t facet_label_count,
                                                             const std::vector<std::vector<size_t>>& pv_rep_lists,
                                                             const size_t originalvtnum, const double maxeps) const
{
    const size_t K = facet_label_count;

    if (covt < originalvtnum)
    {
        ws.singleton_slots[K].assign(1, covt);
        ws.rep_ptrs[K] = &ws.singleton_slots[K];
    }
    else
    {
        ws.rep_ptrs[K] = &pv_rep_lists[covt - originalvtnum];
    }

    const auto& rep_c = *ws.rep_ptrs[K];

    ws.covt_edges.clear();

    for (size_t i = 0; i < K; ++i)
    {
        const auto& rep_i = *ws.rep_ptrs[i];

        for (size_t a = 0; a < rep_i.size(); ++a)
        {
            for (size_t b = 0; b < rep_c.size(); ++b)
            {
                const double w = dist_mat_.getDistance(rep_i[a], rep_c[b]);
                if (w < maxeps)
                {
                    ws.covt_edges.push_back({w, static_cast<uint8_t>(i), static_cast<uint8_t>(K),
                                             static_cast<uint8_t>(a), static_cast<uint8_t>(b)});
                }
            }
        }
    }

    std::sort(ws.covt_edges.begin(), ws.covt_edges.end());
}

template <typename DistMatType>
double SimplexEnumerator<DistMatType>::getGeometricPVSimplexWeight(WitnessWorkspace& ws, const size_t target_simplex_label_count,
                                                                   const double lower_bound, const double maxeps) const
{
    const size_t flattened_adjacency_mask_words = target_simplex_label_count * target_simplex_label_count * MAX_PV_CARDINALITY_;
    if (ws.flattened_adjacency_mask.size() < flattened_adjacency_mask_words)
        ws.flattened_adjacency_mask.resize(flattened_adjacency_mask_words);
    std::fill_n(ws.flattened_adjacency_mask.begin(), flattened_adjacency_mask_words, 0ULL);

    if (ws.candidate_local_index_mask.size() < target_simplex_label_count)
        ws.candidate_local_index_mask.resize(target_simplex_label_count);
    if (ws.current_local_indices.size() < target_simplex_label_count)
        ws.current_local_indices.resize(target_simplex_label_count);
    if (ws.recursion_candidate_local_index_stack.size() < target_simplex_label_count)
        ws.recursion_candidate_local_index_stack.resize(target_simplex_label_count);
    for (auto& level : ws.recursion_candidate_local_index_stack)
    {
        if (level.size() < target_simplex_label_count)
            level.resize(target_simplex_label_count);
    }

    uint64_t* flattened_adjacency_mask = ws.flattened_adjacency_mask.data();

    const auto adjacencyMask = [flattened_adjacency_mask, target_simplex_label_count](const size_t g0, const size_t g1, const size_t local0) -> uint64_t &
    {
        return flattened_adjacency_mask[(g0 * target_simplex_label_count + g1) * MAX_PV_CARDINALITY_ + local0];
    };

    const auto fullMask = [&ws](const size_t g) -> uint64_t
    {
        const size_t repsize = ws.rep_ptrs[g]->size();
        return (repsize >= 64) ? ~0ULL : ((1ULL << repsize) - 1ULL);
    };

    uint64_t coverage = 0ULL;
    bool covered = false;

    const auto& fedges = ws.facet_edges;
    const auto& cedges = ws.covt_edges;
    size_t fi = 0;
    size_t ci = 0;

    while (fi < fedges.size() || ci < cedges.size())
    {
        bool take_facet_edge;
        if (fi >= fedges.size())
            take_facet_edge = false;
        else if (ci >= cedges.size())
            take_facet_edge = true;
        else
            take_facet_edge = (fedges[fi].weight <= cedges[ci].weight);

        const EdgeRecord& edge = take_facet_edge ? fedges[fi++] : cedges[ci++];

        adjacencyMask(edge.virtualidx0, edge.virtualidx1, edge.localidx0) |= (1ULL << edge.localidx1);
        adjacencyMask(edge.virtualidx1, edge.virtualidx0, edge.localidx1) |= (1ULL << edge.localidx0);

        if (!covered)
        {
            coverage |= (1ULL << edge.virtualidx0) | (1ULL << edge.virtualidx1);
            if (static_cast<size_t>(__builtin_popcountll(coverage)) != target_simplex_label_count)
                continue;
            covered = true;
        }

        if (edge.weight < lower_bound)
            continue;

        const size_t g0 = edge.virtualidx0;
        const size_t g1 = edge.virtualidx1;
        const size_t l0 = edge.localidx0;
        const size_t l1 = edge.localidx1;

        // Anchored search: before inserting this edge, no witness clique had been found.
        // Therefore any clique that first appears now must include this newest edge.
        // Pin its two local representatives and pre-intersect every other label's
        // candidate mask with both endpoint adjacency rows before recursing.
        bool feasible = true;
        for (size_t g = 0; g < target_simplex_label_count; ++g)
        {
            if (g == g0)
            {
                ws.candidate_local_index_mask[g] = (1ULL << l0);
                ws.current_local_indices[g] = l0;
                continue;
            }
            if (g == g1)
            {
                ws.candidate_local_index_mask[g] = (1ULL << l1);
                ws.current_local_indices[g] = l1;
                continue;
            }

            ws.current_local_indices[g] = UNCHOSEN_;
            ws.candidate_local_index_mask[g] = fullMask(g) & adjacencyMask(g0, g, l0) & adjacencyMask(g1, g, l1);

            if (ws.candidate_local_index_mask[g] == 0ULL)
            {
                feasible = false;
                break;
            }
        }

        if (feasible && findCliqueRecursive(flattened_adjacency_mask, target_simplex_label_count, ws, ws.candidate_local_index_mask, 2))
            return edge.weight;
    }

    return std::numeric_limits<double>::infinity();
}

template <typename DistMatType>
bool SimplexEnumerator<DistMatType>::findCliqueRecursive(const uint64_t* flattened_adjacency_mask, const size_t target_simplex_label_count, WitnessWorkspace& ws,
                                                         const std::vector<uint64_t>& candidate_local_index_mask,
                                                         const size_t current_local_index_count) const
{
    if (current_local_index_count == target_simplex_label_count)
        return true;

    size_t pivot = UNCHOSEN_;
    int min_candidates = std::numeric_limits<int>::max();

    for (size_t g = 0; g < target_simplex_label_count; ++g)
    {
        if (ws.current_local_indices[g] == UNCHOSEN_)
        {
            const int candnum = __builtin_popcountll(candidate_local_index_mask[g]);
            if (candnum < min_candidates)
            {
                min_candidates = candnum;
                pivot = g;
            }
        }
    }

    if (pivot == UNCHOSEN_)
        return false;

    auto& next_candidate_local_index_mask = ws.recursion_candidate_local_index_stack[current_local_index_count];

    uint64_t opts = candidate_local_index_mask[pivot];
    while (opts != 0ULL)
    {
        const size_t local = static_cast<size_t>(__builtin_ctzll(opts));
        opts &= (opts - 1);

        ws.current_local_indices[pivot] = local;

        bool feasible = true;
        for (size_t g = 0; g < target_simplex_label_count; ++g)
        {
            next_candidate_local_index_mask[g] = candidate_local_index_mask[g];

            if (g == pivot || ws.current_local_indices[g] != UNCHOSEN_)
                continue;

            next_candidate_local_index_mask[g] &= flattened_adjacency_mask[(pivot * target_simplex_label_count + g) * MAX_PV_CARDINALITY_ + local];

            if (next_candidate_local_index_mask[g] == 0ULL)
            {
                feasible = false;
                break;
            }
        }

        if (feasible && findCliqueRecursive(flattened_adjacency_mask, target_simplex_label_count, ws, next_candidate_local_index_mask, current_local_index_count + 1))
            return true;
    }

    ws.current_local_indices[pivot] = UNCHOSEN_;
    return false;
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getGeometricCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex_list,
                                                                                                const std::vector<size_t>& active_vertices,
                                                                                                const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices,
                                                                                                const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash,
                                                                                                const size_t dim, const double maxeps, const int threadnum)
{
    std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(threadnum);
    
    const size_t originalvtnum = dist_mat_.getVertexNumber();
    const size_t npts = originalvtnum + virtual_vertex_indices.size(); // original vertex number + virtual vertex number

    std::vector<std::vector<size_t>> pv_rep_lists(virtual_vertex_indices.size());
    for (size_t i = 0; i < virtual_vertex_indices.size(); ++i)
    {
        pv_rep_lists[i].assign(virtual_vertex_indices[i].begin(), virtual_vertex_indices[i].end());
        std::sort(pv_rep_lists[i].begin(), pv_rep_lists[i].end());
    }

    std::vector<WitnessWorkspace> witness_workspaces(threadnum);

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_virtual_simplex_list.size(); ++i)
    {
        const int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];
        auto& ws = witness_workspaces[threadid];

        const int64_t bindex = sorted_virtual_simplex_list[i].first;
        const double weight = sorted_virtual_simplex_list[i].second;

        SimplexUtility::getSimplexVerticesInPlace(binomial_table_, ws.facet_vertices, bindex, npts, dim);
        const auto& simplex_vertices = ws.facet_vertices;

        const size_t minfacetvt = simplex_vertices.back();

        const auto iter = std::lower_bound(active_vertices.begin(), active_vertices.end(), minfacetvt);

        if (iter == active_vertices.end() || *iter != minfacetvt)
            throw std::out_of_range("vertex not found in active vertex list");

        const size_t vtpos = static_cast<size_t>(std::distance(active_vertices.begin(), iter));

        if (vtpos == 0)
            continue;

        const bool hasvirtual = (simplex_vertices.front() >= originalvtnum);
        const size_t facet_label_count = dim + 1;
        bool witness_context_ready = false;

        for (size_t j = 0; j < vtpos; ++j)
        {
            const size_t covt = active_vertices[j];

            double cofacetweight = 0.0;

            if (hasvirtual)
            {
                double lower_bound = weight;
                for (const auto& vt : simplex_vertices)
                {
                    const double d = SimplexUtility::getLabelDistance(virtual_distance_hash, covt, vt);
                    if (d > lower_bound)
                        lower_bound = d;
                }

                if (lower_bound >= maxeps)
                    continue;

                if (!witness_context_ready)
                {
                    prepareFacetWitnessContext(ws, simplex_vertices, pv_rep_lists, originalvtnum, maxeps);
                    witness_context_ready = true;
                }

                prepareCovtWitnessGroup(ws, covt, facet_label_count, pv_rep_lists, originalvtnum, maxeps);
                cofacetweight = getGeometricPVSimplexWeight(ws, facet_label_count + 1, lower_bound, maxeps);

                if (!(cofacetweight < maxeps))
                    continue;
            }
            else
            {
                cofacetweight = weight;
                for (const auto& vt : simplex_vertices)
                    cofacetweight = std::max(cofacetweight, dist_mat_.getDistance(covt, vt));
            }

            if (cofacetweight > 0.0 && cofacetweight < maxeps)
            {
                ws.cofacet_vertices.assign(simplex_vertices.begin(), simplex_vertices.end());
                ws.cofacet_vertices.push_back(covt); // still descending order because covt < minfacetvt

                const int64_t cofacetbindex = SimplexUtility::getBinomialIndex(binomial_table_, ws.cofacet_vertices, 0);
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    return SimplexUtility::sortAndMergeSimplexChunks(thread_workspace, threadnum);
}
