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
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVRCofacets(const std::vector<std::pair<int64_t, double>> &sorted_simplex_list, const size_t dim, const double maxeps, const int threadnum)
{
    // dim == simplex dimension == cofacet dimension - 1
    std::vector<std::pair<int64_t, double>> cofacet_list;

    size_t npts = binomial_table_.size() - 1;

    std::vector<std::vector<std::pair<size_t, double>>> thread_workspace(threadnum);

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_simplex_list.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto &thread_cofacets = thread_workspace[threadid];

        const int64_t bindex = sorted_simplex_list[i].first;
        const double weight = sorted_simplex_list[i].second;

        std::vector<size_t> simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        const size_t minvt = simplex_vertices.back(); // sorted in descending order
        for (size_t covt = 0; covt < minvt; ++covt)
        {
            double newweight = 0.0;
            for (const auto &vt : simplex_vertices)
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

    for (const auto &thread_cofacets : thread_workspace)
    {
        cofacet_list.insert(cofacet_list.end(), thread_cofacets.begin(), thread_cofacets.end());
    }

    SimplexUtility::sortSimplexByWeightThenIndex(cofacet_list);

    return cofacet_list;
}

#ifdef BUILD_ALPHA_COMPLEX
template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedAlphaCells(const std::vector<std::vector<int64_t>> &binomial_table,
                                                                                            std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t> &vertex_handle_index,
                                                                                            CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> &delaunay_d, const size_t dim, double maxeps)
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
double SimplexEnumerator<DistMatType>::getAlphaSimplexWeight(const std::vector<size_t> &alpha_simplex)
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
bool SimplexEnumerator<DistMatType>::findCliqueRecursive(const std::vector<std::vector<std::vector<uint64_t>>>& adj_mask,
                                                         std::vector<uint64_t> &candidate_mask,
                                                         std::vector<size_t> &current_clique_local_indices,
                                                         size_t depth)
{
    const size_t K = candidate_mask.size(); // clique size

    if (depth == K)
        return true;  // Successfully found a clique of the target size

    // Heuristic: branch on the set with the fewest remaining candidates.
    size_t pivot = std::numeric_limits<size_t>::max();
    size_t mincandidates = std::numeric_limits<size_t>::max();

    const size_t UNCHOSEN = std::numeric_limits<size_t>::max();

    for (size_t i = 0; i < K; ++i)
    {
        if (current_clique_local_indices[i] == UNCHOSEN)
        {
            const int candnum = __builtin_popcountll(candidate_mask[i]);
            if (candnum > 0 && static_cast<size_t>(candnum) < mincandidates)
            {
                mincandidates = static_cast<size_t>(candnum);
                pivot = i;
            }
        }
    }

    if (pivot == std::numeric_limits<size_t>::max())
        return false;  // no valid candidates left

    uint64_t opts = candidate_mask[pivot];
    while (opts)
    {
        const int localidx = __builtin_ctzll(opts);  // start from the local index of the vertex represented by the lowest bit. <64
        opts &= (opts - 1);  // pop the lowest bit

        current_clique_local_indices[pivot] = static_cast<size_t>(localidx);

        std::vector<uint64_t> next_candidate_mask = candidate_mask;
        bool feasible = true;

        // early prune
        for (size_t j = 0; j < K; ++j)
        {
            if (j == pivot || current_clique_local_indices[j] != UNCHOSEN)
                continue;

            next_candidate_mask[j] &= adj_mask[pivot][j][static_cast<size_t>(localidx)];
            if (next_candidate_mask[j] == 0ULL)
            {
                feasible = false;
                break;
            }
        }

        if (feasible && findCliqueRecursive(adj_mask, next_candidate_mask, current_clique_local_indices, depth + 1))
            return true;
    }

    // unfeasible with current options, backtrack
    current_clique_local_indices[pivot] = UNCHOSEN;
    return false;
}

template <typename DistMatType>
double SimplexEnumerator<DistMatType>::getGeometricVirtualSimplexWeight(const std::vector<size_t>& simplex_vertices,
                                                                        const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices,
                                                                        size_t dim)
{
    const size_t UNCHOSEN = std::numeric_limits<size_t>::max();

    const size_t originalvtnum = dist_mat_.getVertexNumber();
    const size_t K = dim + 1;

    if (K < 3)
        return 0.0;

    // Expand each simplex vertex into a list of possible original representatives.
    std::vector<std::vector<size_t>> rep_lists(K);
    for (size_t i = 0; i < K; ++i)
    {
        const size_t vtidx = simplex_vertices[i];

        if (vtidx < originalvtnum)
        {
            rep_lists[i] = {vtidx};
        }
        else
        {
            const auto& idx_set = virtual_vertex_indices[vtidx - originalvtnum];
            rep_lists[i].assign(idx_set.begin(), idx_set.end());
        }
    }

    // Collect and sort all possible cross-edges between representative choices.
    std::vector<EdgeRecord> sorted_edges;
    sorted_edges.reserve(K * (K - 1) / 2 * 16); // mild hint; exact depends on PV sizes

    for (uint8_t i = 0; i < K; ++i)
    {
        for (uint8_t j = i + 1; j < K; ++j)
        {
            for (uint8_t locali = 0; locali < rep_lists[i].size(); ++locali)
            {
                for (uint8_t localj = 0; localj < rep_lists[j].size(); ++localj)
                {
                    const double w = dist_mat_.getDistance(rep_lists[i][locali], rep_lists[j][localj]);
                    sorted_edges.push_back({w, i, j, locali, localj});  // aggregate/list initialization
                }
            }
        }
    }

    if (sorted_edges.empty())
        return std::numeric_limits<double>::infinity();

    std::sort(sorted_edges.begin(), sorted_edges.end());

    // adjacency masks: adj_mask[a][b][local_a] is a bitmask of b's vertices adjacent to local_a
    std::vector<std::vector<std::vector<uint64_t>>> adj_mask(K, std::vector<std::vector<uint64_t>>(K));
    for (size_t i = 0; i < K; ++i)
    {
        for (size_t j = 0; j < K; ++j)
        {
            if (i == j)
                continue;
            adj_mask[i][j].assign(rep_lists[i].size(), 0ULL);
        }
    }

    // Delay the expensive clique check until every group has at least one incident edge.
    uint64_t coverage = 0ULL;
    bool covered = false;

    for (const auto &edge : sorted_edges)
    {
        adj_mask[edge.virtualidx0][edge.virtualidx1][edge.localidx0] |= (1ULL << edge.localidx1);
        adj_mask[edge.virtualidx1][edge.virtualidx0][edge.localidx1] |= (1ULL << edge.localidx0);

        if (!covered)
        {
            coverage |= (1ULL << edge.virtualidx0) | (1ULL << edge.virtualidx1);
            if (static_cast<size_t>(__builtin_popcountll(coverage)) != K)
                continue;
            covered = true;
        }

        std::vector<uint64_t> candidate_mask(K);
        for (size_t i = 0; i < K; ++i)
        {
            // NOTE: PV size is capped (<64) upstream. If this ever becomes 64, this shift is UB.
            candidate_mask[i] = (1ULL << rep_lists[i].size()) - 1ULL;
        }

        std::vector<size_t> current_clique_local_indices(K, UNCHOSEN);

        if (findCliqueRecursive(adj_mask, candidate_mask, current_clique_local_indices, 0))
            return edge.weight;
    }

    return std::numeric_limits<double>::infinity();  // no clique found
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getGeometricVirtualCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex_list,
                                                                                                       const std::vector<size_t>& active_vertices,
                                                                                                       const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices,
                                                                                                       const size_t dim, const double maxeps, const int threadnum)
{
    std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(threadnum);
    
    const size_t originalvtnum = dist_mat_.getVertexNumber();
    const size_t npts = originalvtnum + virtual_vertex_indices.size(); // original vertex number + virtual vertex number
    

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_virtual_simplex_list.size(); ++i)
    {
        const int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];

        const int64_t bindex = sorted_virtual_simplex_list[i].first;
        const double weight = sorted_virtual_simplex_list[i].second;

        auto simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        const size_t minfacetvt = simplex_vertices.back();

        const auto iter = std::lower_bound(active_vertices.begin(), active_vertices.end(), minfacetvt);

        if (iter == active_vertices.end() || *iter != minfacetvt)
            throw std::out_of_range("vertex not found in active vertex list");

        const size_t vtpos = static_cast<size_t>(std::distance(active_vertices.begin(), iter));

        const bool hasvirtual = (simplex_vertices.front() >= originalvtnum);

        for (size_t j = 0; j < vtpos; ++j)
        {
            const size_t covt = active_vertices[j];

            std::vector<size_t> cofacet_vertices = simplex_vertices;
            cofacet_vertices.push_back(covt); // still descending order because covt < minfacetvt

            double cofacetweight = 0.0;

            if (hasvirtual)
            {
                cofacetweight = getGeometricVirtualSimplexWeight(cofacet_vertices, virtual_vertex_indices, dim + 1);
            }
            else
            {
                cofacetweight = weight;
                for (const auto& vt : simplex_vertices)
                    cofacetweight = std::max(cofacetweight, dist_mat_.getDistance(covt, vt));
            }

            if (cofacetweight > 0.0 && cofacetweight < maxeps)
            {
                const int64_t cofacetbindex = SimplexUtility::getBinomialIndex(binomial_table_, cofacet_vertices, 0);
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    std::vector<std::pair<int64_t, double>> sorted_cofacets;
    for (const auto& thread_cofacets : thread_workspace)
        sorted_cofacets.insert(sorted_cofacets.end(), thread_cofacets.begin(), thread_cofacets.end());

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_cofacets);

    return sorted_cofacets;
}
