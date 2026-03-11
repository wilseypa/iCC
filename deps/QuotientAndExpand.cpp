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
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runQuotient(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> untrimed_pv_index_sets = getPVIndexSets(maxdim, initeps, threadnumber);

    std::vector<std::unordered_set<size_t>> pv_index_sets = trimIndexSets(untrimed_pv_index_sets);

    for (auto& pv_vts : pv_index_sets)
    {
        std::cout << "size of the vt = " << pv_vts.size() << "  contents :  ";
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

        auto critsimpnum = morse_matching.implicitMatch(matching_context, dim_persistent_pairs);
        dim_persistent_pairs.clear();

        // std::cout << "dim = " << dim << "  critical simplex number: " << critsimpnum << std::endl;

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

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getPVIndexSets(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> virtual_vertex_indices;

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
            // std::cout << "critical simplex number: " << critsimpnum << std::endl;
            dim_persistent_pairs.clear();
        }
        else
        {
            auto pv_support_cofacets = morse_matching.implicitMatchAndCollectPVSupports(matching_context);

            std::cout << "pv support set num = " << pv_support_cofacets.size() << '\n';

            virtual_vertex_indices = getGradientPathVertexSets(matching_context, pv_support_cofacets, dim);
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

    return virtual_vertex_indices;
}

template <typename DistMatType>
std::vector<std::vector<size_t>> QuotientAndExpand<DistMatType>::extractGradientPaths(const MatchingContext& matching_context, const double minfacetweight)
{
    auto& bi_graph = matching_context.graph;
    auto& sorted_facets = matching_context.sorted_facets;
    auto& sorted_cofacets = matching_context.sorted_cofacets;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& binom_table = matching_context.binomial_table;

    const size_t u = bi_graph.unodes;
    const size_t v = bi_graph.vnodes;

    const size_t dim = matching_context.dim;   // cofacet dimension
    const size_t npts = matching_context.npts; // total vertex number

    // cutoff by weight
    auto facetcutoff = std::lower_bound(sorted_facets.begin(), sorted_facets.end(), minfacetweight,
                                        [](const std::pair<int64_t, double> &pair, double weight)
                                        { return pair.second < weight; });

    const size_t maxfacetindex = static_cast<size_t>(std::distance(sorted_facets.begin(), facetcutoff));

    auto cofacetcutoff = std::lower_bound(sorted_cofacets.begin(), sorted_cofacets.end(), minfacetweight,
                                          [](const std::pair<int64_t, double> &pair, double weight)
                                          { return pair.second < weight; });

    const size_t maxcofacetindex = static_cast<size_t>(std::distance(sorted_cofacets.begin(), cofacetcutoff));

    if (maxfacetindex == 0 || maxcofacetindex == 0)
        return {};

    // Union-find over cofacet indices in [0, maxcofacetindex)
    UnionFind union_find(maxcofacetindex);

    // mark cofacets that have been claimed into some gradient-path component
    std::vector<uint8_t> is_in_gradient_path(maxcofacetindex, 0);

    // workspaces
    std::vector<size_t> facet_indices;
    facet_indices.reserve(dim + 1);
    std::vector<size_t> vertex_workspace;
    vertex_workspace.reserve(dim + 1);

    // scan facets (prefix under cutoff); connect matched cofacets that are linked through active facets
    for (size_t facetidx = 0; facetidx < maxfacetindex; ++facetidx)
    {
        const int64_t cofacetidxint = bi_graph.match_list[facetidx + u];

        if (cofacetidxint < 0)
            continue;

        const size_t cofacetidx = static_cast<size_t>(cofacetidxint);

        if (cofacetidx >= maxcofacetindex)
            continue;

        is_in_gradient_path[cofacetidx] = 1;

        // get active facet neighbors of this cofacet (implicit traversal)
        SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices, vertex_workspace,
                                                   sorted_cofacets[cofacetidx].first, npts, dim);

        for (const auto nextfacetidx : facet_indices)
        {
            if (nextfacetidx == facetidx || nextfacetidx >= maxfacetindex)
                continue;

            const int64_t nextcofacetidxint = bi_graph.match_list[u + nextfacetidx];
            if (nextcofacetidxint < 0)
                continue;

            const size_t nextcofacetidx = static_cast<size_t>(nextcofacetidxint);

            if (nextcofacetidx >= maxcofacetindex)
                continue;

            if (is_in_gradient_path[nextcofacetidx])
                continue; // keep components disjoint

            is_in_gradient_path[nextcofacetidx] = 1;
            union_find.unionSets(cofacetidx, nextcofacetidx);
        }
    }

    // collect components
    std::unordered_map<size_t, std::vector<size_t>> component_map;
    component_map.reserve(maxcofacetindex);

    for (size_t cofacetidx = 0; cofacetidx < maxcofacetindex; ++cofacetidx)
    {
        if (is_in_gradient_path[cofacetidx])
        {
            const size_t root = union_find.unionFind(cofacetidx);
            component_map[root].push_back(cofacetidx);
        }
    }

    std::vector<std::vector<size_t>> gradient_paths;
    gradient_paths.reserve(component_map.size());

    for (auto& [root, path] : component_map) // structured binding requires C++17
    {

        // std::cout<<"in extract grad path function, path size = "<<path.size()<<'\n';

        if (path.size() > 1) // minimum number of cofacets in the path
        {
            gradient_paths.push_back(std::move(path));
        }
    }

    return gradient_paths;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getGradientPathVertexSets(const MatchingContext &matching_context, const std::vector<std::vector<size_t>>& pv_support_cofacets, const size_t dim)
{
    std::vector<std::unordered_set<size_t>> gradient_path_vertex_sets;
    gradient_path_vertex_sets.reserve(pv_support_cofacets.size());

    auto npt = matching_context.binomial_table.size() - 1;

    for (const auto& pv_support : pv_support_cofacets)
    {
        //need set for trimming
        std::unordered_set<size_t> vertex_set;

        for (auto cofacetidx : pv_support)
        {
            auto bindex = matching_context.sorted_cofacets[cofacetidx].first;
            auto simplex_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, bindex, npt, dim);
            vertex_set.insert(simplex_vertices.begin(), simplex_vertices.end());
        }
        if (vertex_set.size() < MAX_SIZE_)
            gradient_path_vertex_sets.push_back(std::move(vertex_set));
    }

    return gradient_path_vertex_sets;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::trimIndexSets(std::vector<std::unordered_set<size_t>> &gradient_path_vertex_sets)
{
    std::sort(gradient_path_vertex_sets.begin(), gradient_path_vertex_sets.end(),
              [](const std::unordered_set<size_t> &lhs, const std::unordered_set<size_t> &rhs)
              { return lhs.size() > rhs.size(); });

    std::vector<std::unordered_set<size_t>> trimmed_vertex_sets;

    std::unordered_set<size_t> claimed_vertices;

    // const size_t MAX_VERTICES_NUM = MAX_SIZE_ - 1;    //max number of virtual vertex. the max value of uint8_t used in EdgeRecord

    for (auto& vertex_set : gradient_path_vertex_sets)
    {
        // not necessary
        // if (trimmed_vertex_sets.size() >= MAX_VERTICES_NUM) 
        // {
        //     std::cout << "Warning: Reached maximum number of virtual vertices (" << MAX_VERTICES_NUM << "). Discarding remaining candidates." << std::endl;
        //     break;
        // }

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

// template <typename DistMatType>
// std::vector<std::pair<int64_t, double>> QuotientAndExpand<DistMatType>::getVirtualCofacetList(const std::vector<std::pair<int64_t, double>> &sorted_virtual_simplex_list,
//                                                                                               const std::vector<size_t> &active_vertices, const robin_hood::unordered_map<uint64_t, double> &virtual_distance_hash_table,
//                                                                                               const size_t dim, const double maxeps, int threadnum)
// {
//     std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(threadnum);

//     const size_t npts = binomial_table_.size() - 1;

//     omp_set_num_threads(threadnum);

// #pragma omp parallel for schedule(dynamic)
//     for (size_t i = 0; i < sorted_virtual_simplex_list.size(); ++i)
//     {
//         int threadid = omp_get_thread_num();
//         auto &thread_cofacets = thread_workspace[threadid];

//         const auto &simplex_pair = sorted_virtual_simplex_list[i];
//         const int64_t bindex = simplex_pair.first;
//         const double weight = simplex_pair.second;

//         auto simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

//         auto iter = std::find(active_vertices.begin(), active_vertices.end(), simplex_vertices.back());

//         if (iter == active_vertices.end())
//             throw std::out_of_range("vertex not found in active vertex list");

//         auto vtpos = std::distance(active_vertices.begin(), iter);

//         for (size_t j = 0; j < vtpos; ++j)
//         {
//             const size_t covt = active_vertices[j];

//             double newweight = 0.0;
//             for (const auto &vt : simplex_vertices)
//             {
//                 newweight = std::max(newweight, getVirtualDistance(covt, vt, virtual_distance_hash_table));
//             }

//             double cofacetweight = std::max(newweight, weight);
//             // std::cout<<"potiential cofacet weight = "<<cofacetweight<<'\n';

//             if (cofacetweight < maxeps)
//             {
//                 int64_t shiftedbindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vertices, 1);
//                 int64_t cofacetbindex = shiftedbindex + covt;
//                 thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
//             }
//         }
//     }

//     std::vector<std::pair<int64_t, double>> sorted_cofacets;
//     for (const auto &thread_cofacets : thread_workspace)
//     {
//         sorted_cofacets.insert(sorted_cofacets.end(), thread_cofacets.begin(), thread_cofacets.end());
//     }

//     SimplexUtility::sortSimplexByWeightThenIndex(sorted_cofacets);

//     return sorted_cofacets;
// }

// template <typename DistMatType>
// bool QuotientAndExpand<DistMatType>::findCliqueRecursive(const std::vector<std::vector<std::vector<uint64_t>>> &adj_mask,
//                                                          std::vector<uint64_t> &candidate_mask,
//                                                          std::vector<size_t> &current_clique_local_indices, size_t depth)
// {
//     const size_t K = candidate_mask.size(); // clique size

//     if (depth == K)
//         return true; // Successfully found a clique of the target size

//     // heuristic: pick the set with the fewest remaining candidates to branch on next.
//     size_t pivot = std::numeric_limits<size_t>::max();
//     size_t mincandidates = std::numeric_limits<size_t>::max(); // size cap of virtual vertex = 64

//     const size_t UNCHOSEN = std::numeric_limits<size_t>::max();

//     for (auto i = 0; i < K; ++i)
//     {
//         if (current_clique_local_indices[i] == UNCHOSEN)
//         {
//             int candnum = __builtin_popcountll(candidate_mask[i]); // popcount returns int

//             if (candnum > 0 && candnum < mincandidates)
//             {
//                 mincandidates = candnum;
//                 pivot = i;
//             }
//         }
//     }

//     if (pivot == std::numeric_limits<size_t>::max())
//         return false; // no valid candidates left

//     uint64_t opts = candidate_mask[pivot];
//     while (opts)
//     {
//         int localidx = __builtin_ctzll(opts); // local index of the vertex represented by the lowest bit. <64
//         opts &= (opts - 1);                   // pop the lowest bit

//         current_clique_local_indices[pivot] = static_cast<size_t>(localidx);
//         std::vector<uint64_t> next_candidate_mask = candidate_mask;
//         bool feasible = true;

//         // early prune the candidates
//         for (auto j = 0; j < K; ++j)
//         {
//             if (j == pivot || current_clique_local_indices[j] != UNCHOSEN)
//                 continue;

//             next_candidate_mask[j] &= adj_mask[pivot][j][static_cast<size_t>(localidx)];
//             if (next_candidate_mask[j] == 0ULL)
//             {
//                 feasible = false;
//                 break;
//             }
//         }

//         if (feasible && findCliqueRecursive(adj_mask, next_candidate_mask, current_clique_local_indices, depth + 1))
//             return true;
//     }

//     // unfeasible with current options, backtrack
//     current_clique_local_indices[pivot] = UNCHOSEN;
//     return false;
// }

// template <typename DistMatType>
// double QuotientAndExpand<DistMatType>::getGeometricVirtualSimplexWeight(const std::vector<size_t> &simplex_vertices,
//                                                                         const std::vector<std::unordered_set<size_t>> &virtual_vertex_indices, size_t dim)
// {
//     // const size_t MAXSIZE = 64;

//     const size_t UNCHOSEN = std::numeric_limits<size_t>::max();

//     const size_t originalvtnum = dist_mat_.getVertexNumber();
//     const size_t K = dim + 1;
//     if (K < 3)
//         return 0;

//     // prepare the indices
//     std::vector<std::vector<size_t>> vvt_idx_vec(K);
//     for (auto i = 0; i < K; ++i)
//     {
//         size_t vtidx = simplex_vertices[i];

//         if (vtidx < originalvtnum) // regular vertex
//         {
//             vvt_idx_vec[i] = {vtidx};
//         }
//         else // virtual vertex vvt
//         {
//             const auto &idx_set = virtual_vertex_indices[vtidx - originalvtnum];
//             // if (idx_set.size() > MAXSIZE) throw std::runtime_error("Virtual vertex set size exceeds the limit.");
//             vvt_idx_vec[i].assign(idx_set.begin(), idx_set.end());
//         }
//     }

//     // collect and sort all the possible edges
//     std::vector<EdgeRecord> sorted_edges;
//     for (uint8_t i = 0; i < K; ++i)
//     {
//         for (uint8_t j = i + 1; j < K; ++j)
//         {
//             for (uint8_t locali = 0; locali < vvt_idx_vec[i].size(); ++locali)
//             {
//                 for (uint8_t localj = 0; localj < vvt_idx_vec[j].size(); ++localj)
//                 {
//                     double weight = dist_mat_.getDistance(vvt_idx_vec[i][locali], vvt_idx_vec[j][localj]);
//                     sorted_edges.push_back({weight, i, j, locali, localj}); // aggregate/list initialization
//                 }
//             }
//         }
//     }
//     std::sort(sorted_edges.begin(), sorted_edges.end());

//     if (sorted_edges.empty())
//         return std::numeric_limits<double>::infinity();

//     // incrementally check for a clique
//     // adj mask size == K * K * VT_SIZE_MAX * uint64_t
//     std::vector<std::vector<std::vector<uint64_t>>> adj_mask(K, std::vector<std::vector<uint64_t>>(K));
//     for (auto i = 0; i < K; ++i)
//     {
//         for (auto j = 0; j < K; ++j)
//         {
//             if (i != j)
//                 adj_mask[i][j].resize(vvt_idx_vec[i].size(), 0);
//         }
//     }

//     // coverage mask to delay the recursive search
//     uint64_t coverage = 0;
//     bool covered = false;

//     auto tempsize = sorted_edges.size();

//     for (const auto &edge : sorted_edges)
//     {
//         // update adj mask with new edge
//         adj_mask[edge.virtualidx0][edge.virtualidx1][edge.localidx0] |= (1ULL << edge.localidx1);
//         adj_mask[edge.virtualidx1][edge.virtualidx0][edge.localidx1] |= (1ULL << edge.localidx0);

//         // coverage check
//         if (!covered)
//         {
//             coverage |= (1ULL << edge.virtualidx0) | (1ULL << edge.virtualidx1);
//             if (static_cast<size_t>(__builtin_popcountll(coverage)) == K)
//             {
//                 covered = true;
//             }
//             else
//             {
//                 continue;
//             }
//         }

//         std::vector<uint64_t> candidate_mask(K);
//         for (auto i = 0; i < K; ++i)
//         {
//             candidate_mask[i] = (1ULL << vvt_idx_vec[i].size()) - 1; // init mask, all available
//         }

//         std::vector<size_t> current_clique_local_indices(K, UNCHOSEN);

//         // std::cout<<"at dim = "<<dim<<"  virtual simp weight recursive call will start with sorted_edges size = "<<tempsize<<'\n';

//         if (findCliqueRecursive(adj_mask, candidate_mask, current_clique_local_indices, 0))
//             return edge.weight;
//     }

//     return std::numeric_limits<double>::infinity(); // no clique found
// }

// template <typename DistMatType>
// std::vector<std::pair<int64_t, double>> QuotientAndExpand<DistMatType>::getGeometricVirtualCofacetList(const std::vector<std::pair<int64_t, double>> &sorted_virtual_simplex_list,
//                                                                                                        const std::vector<size_t> &active_vertices,
//                                                                                                        const std::vector<std::unordered_set<size_t>> &virtual_vertex_indices,
//                                                                                                        const size_t dim, const double maxeps, int threadnum)
// {
//     std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(threadnum);

//     const size_t npts = binomial_table_.size() - 1; // original vertex number + virtual vertex number

//     const size_t originalvtnum = dist_mat_.getVertexNumber();

//     omp_set_num_threads(threadnum);

// #pragma omp parallel for schedule(dynamic)
//     for (auto i = 0; i < sorted_virtual_simplex_list.size(); ++i)
//     {
//         int threadid = omp_get_thread_num();
//         auto &thread_cofacets = thread_workspace[threadid];

//         const auto &simplex_pair = sorted_virtual_simplex_list[i];
//         const int64_t bindex = simplex_pair.first;
//         const double weight = simplex_pair.second;
//         auto simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

//         const size_t minfacetvt = simplex_vertices.back();
//         auto iter = std::find(active_vertices.begin(), active_vertices.end(), minfacetvt);
//         if (iter == active_vertices.end())
//             throw std::out_of_range("vertex not found in active vertex list");
//         auto vtpos = std::distance(active_vertices.begin(), iter);

//         bool hasvirtual = (simplex_vertices.front() >= originalvtnum) ? true : false;

//         for (auto j = 0; j < vtpos; ++j)
//         {
//             const size_t covt = active_vertices[j];

//             std::vector<size_t> cofacet_vertices = simplex_vertices;
//             cofacet_vertices.push_back(covt);

//             double cofacetweight = 0.0;

//             if (hasvirtual)
//             {
//                 cofacetweight = getGeometricVirtualSimplexWeight(cofacet_vertices, virtual_vertex_indices, dim + 1);
//             }
//             else
//             {
//                 // all real/regular vertices
//                 // check distances from the new co-vertex to all vertices of the original facet
//                 cofacetweight = weight;
//                 for (const auto &vt : simplex_vertices)
//                 {
//                     cofacetweight = std::max(cofacetweight, dist_mat_.getDistance(covt, vt));
//                 }
//             }

//             if (cofacetweight > 0 && cofacetweight < maxeps)
//             {
//                 int64_t cofacetbindex = SimplexUtility::getBinomialIndex(binomial_table_, cofacet_vertices, 0);
//                 thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
//             }
//         }
//     }

//     std::vector<std::pair<int64_t, double>> sorted_cofacets;
//     for (const auto &thread_cofacets : thread_workspace)
//     {
//         sorted_cofacets.insert(sorted_cofacets.end(), thread_cofacets.begin(), thread_cofacets.end());
//     }

//     SimplexUtility::sortSimplexByWeightThenIndex(sorted_cofacets);

//     return sorted_cofacets;
// }

// template <typename DistMatType>
// void QuotientAndExpand<DistMatType>::buildInterface(BipartiteGraph &bi_graph, const std::vector<std::pair<int64_t, double>> &sorted_cofacet_list,
//                                                     const robin_hood::unordered_map<int64_t, size_t> &active_facet_index_hash_table, const size_t dim)
// {
//     auto u = bi_graph.unodes;
//     auto v = bi_graph.vnodes;

//     // std::cout<<u<<"    "<<v<<'\n';

//     for (size_t i = 0; i < sorted_cofacet_list.size(); i++)
//     {
//         int64_t bindex = sorted_cofacet_list[i].first;
//         std::vector<int64_t> facet_bindex_list = SimplexUtility::getFacetBinomialIndices(binomial_table_, bindex, dim);

//         for (auto fbidx : facet_bindex_list)
//         {
//             auto mapit = active_facet_index_hash_table.find(fbidx);
//             if (mapit != active_facet_index_hash_table.end())
//             {
//                 size_t fidx = (*mapit).second; // facet idx in facet bindex_weight_pair list
//                 bi_graph.adj_list[i].push_back(u + fidx);
//                 bi_graph.adj_list[u + fidx].push_back(i);
//             }
//         }
//         std::sort(bi_graph.adj_list[i].begin(), bi_graph.adj_list[i].end(), std::greater<size_t>());
//     }

//     return;
// }