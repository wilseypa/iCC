#include <unordered_map>

#include "omp.h"

#include "SimplexUtility.hpp"
#include "SimplexEnumerator.hpp"
#include "MaximumMorseMatching.hpp"
#include "QuotientAndExpand.hpp"

template class QuotientAndExpand<NormalDistMat>;
template class QuotientAndExpand<SparseDistMat>;

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::runQuotient(const size_t maxdim, const double initeps, const int threadnumber)
{
    std::vector<std::unordered_set<size_t>> untrimed_virtual_vertex_indices = getVirtualVertexIndices(maxdim, initeps, threadnumber);

    std::vector<std::unordered_set<size_t>> virtual_vertex_indices = trimVertexSets(untrimed_virtual_vertex_indices);

    return virtual_vertex_indices;
}

template <typename DistMatType>
void QuotientAndExpand<DistMatType>::runExpand(const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices, const size_t maxdim, const double maxeps, const int threadnumber)
{
    size_t originalvtnum = dist_mat_.dist_mat.size();    //number of original vertices
    size_t virtualvtnum = virtual_vertex_indices.size();

    SimplexUtility::updateBinomialTable(binomial_table_, originalvtnum, virtualvtnum, maxdim);

    auto active_vertices = getActiveVertexIndices(virtual_vertex_indices);

    auto virtual_distance_hash_table = getVirtualDistanceHashTable(active_vertices, virtual_vertex_indices, threadnumber);

    auto sorted_virtual_simplex = getVirtualSortedEdge(active_vertices, virtual_distance_hash_table, maxeps, threadnumber);

    auto active_index_hash_table = getVirtualActiveEdgeIndexHashTable(sorted_virtual_simplex, virtualvtnum);

    BipartiteGraph bi_graph(1, 1, 1);

    std::vector<std::pair<double, double>> dim_persistent_pair;

    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        auto sorted_virtual_cofacet = getVirtualSortedCofacetList(sorted_virtual_simplex, active_vertices, virtual_distance_hash_table, dim, maxeps, threadnumber);

        bi_graph.updateDimension(sorted_virtual_cofacet.size(), sorted_virtual_simplex.size());

        buildInterface(bi_graph, sorted_virtual_cofacet, active_index_hash_table, dim);

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_virtual_simplex, sorted_virtual_cofacet);

        MaximumMorseMatching morse_matching(threadnumber);

        auto critsimpnum = morse_matching.matchWithPersistence(matching_context, dim_persistent_pair);
        std::cout << "dim = "<<dim<<"  critical simplex number: " << critsimpnum << std::endl;
        dim_persistent_pair.clear();

        if (dim != maxdim)
        {
            active_index_hash_table = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_virtual_cofacet);

            sorted_virtual_simplex = getVirtualSortedCofacetList(sorted_virtual_simplex, active_vertices, virtual_distance_hash_table, dim + 1, maxeps, threadnumber);
            std::swap(sorted_virtual_simplex, sorted_virtual_cofacet);
        }
    }

    return;
}


template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getVirtualVertexIndices(const size_t maxdim, const double initeps, const int threadnumber)
{
    size_t originalvtnum = dist_mat_.dist_mat.size();    //number of original vertices

    SimplexEnumerator<DistMatType> simplex_enumerator(dist_mat_, binomial_table_);
    
    auto sorted_simplex = simplex_enumerator.getSortedVREdges(initeps);

    auto active_index_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binomial_table_, sorted_simplex, originalvtnum);

    auto sorted_cofacet = simplex_enumerator.getSortedVRCofacets(binomial_table_, sorted_simplex, 1, initeps, threadnumber);

    BipartiteGraph bi_graph(1, 1, 1);

    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        bi_graph.updateDimension(sorted_cofacet.size(), sorted_simplex.size());
        
        buildInterface(bi_graph, sorted_cofacet, active_index_hash_table, dim);

        MatchingContext matching_context(bi_graph, binomial_table_, sorted_simplex, sorted_cofacet);

        MaximumMorseMatching morse_matching(threadnumber);

        std::vector<std::pair<double, double>> dim_persistent_pair;

        if (dim != maxdim)
        {
            dim_persistent_pair.clear();
            auto critsimpnum = morse_matching.matchWithPersistence(matching_context, dim_persistent_pair);
            std::cout << "critical simplex number: " << critsimpnum << std::endl;
        }
        else
        {
            dim_persistent_pair.clear();
            auto critidx = morse_matching.matchWithPersistenceAndReturnMinCriticalIndex(matching_context, dim_persistent_pair);

            auto minfacetindex = critidx - bi_graph.unodes;

            double minfacetweight = -1.0;
            
            if (minfacetindex >= 0) minfacetweight = sorted_simplex[minfacetindex].second;

            if (minfacetweight < 0.0) minfacetweight = initeps;

            //get gradient paths
            auto gradient_paths = extractGradientPahts(matching_context, minfacetweight);
            //get gradient path vertex indices
            auto virtual_vertex_indices = getGradientPathVertexSets(matching_context, gradient_paths, dim);
        }

        if (dim != maxdim)
        {
            active_index_hash_table = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_cofacet);

            sorted_simplex = simplex_enumerator.getSortedVRCofacets(binomial_table_, sorted_cofacet, dim, initeps, threadnumber);
            std::swap(sorted_simplex, sorted_cofacet);
        }

    }

    return virtual_vertex_indices;
}


template <typename DistMatType>
std::vector< std::vector<size_t> > QuotientAndExpand<DistMatType>::extractGradientPahts(const MatchingContext& matching_context, const double minfacetweight)
{
    auto bi_graph = matching_context.graph;
    auto sorted_simplex = matching_context.sorted_facets;
    auto sorted_cofacet = matching_context.sorted_cofacets;

    const size_t u = bi_graph.unodes;
    const size_t v = bi_graph.vnodes;

    auto facetcutoff = std::lower_bound(sorted_simplex.begin(), sorted_simplex.end(), minfacetweight,
        [](const std::pair<int64_t, double>& pair, double weight) { return pair.second < weight; });

    size_t maxfacetindex = std::distance(sorted_simplex.begin(), facetcutoff);

    auto cofacetcutoff = std::lower_bound(sorted_cofacet.begin(), sorted_cofacet.end(), minfacetweight,
        [](const std::pair<int64_t, double>& pair, double weight) { return pair.second < weight; });

    size_t maxcofacetindex = std::distance(sorted_cofacet.begin(), cofacetcutoff);

    UnionFind union_find(u);

    std::vector<bool> is_in_gradient_path(u, false);

    for (size_t facetidx = 0; facetidx < maxfacetindex; ++facetidx)
    {
        size_t cofacetidx = bi_graph.match_list[facetidx + u];

        if (cofacetidx < 0 || cofacetidx >= maxcofacetindex) continue;

        is_in_gradient_path[cofacetidx] = true;

        for (auto vidx : bi_graph.adj_list[cofacetidx])
        {
            size_t nextfacetidx = vidx - u;

            if (nextfacetidx == facetidx || nextfacetidx >= maxfacetindex) continue;

            size_t nextcofacetidx = bi_graph.match_list[vidx];

            if (nextcofacetidx < 0 || nextcofacetidx >= maxcofacetindex) continue;

            is_in_gradient_path[nextcofacetidx] = true;
            union_find.unionSets(cofacetidx, nextcofacetidx);
        }
    }

    std::unordered_map<size_t, std::vector<size_t>> component_map;

    for (size_t cofacetidx = 0; cofacetidx < maxcofacetindex; ++cofacetidx)
    {
        if (is_in_gradient_path[cofacetidx])
        {
            size_t root = union_find.unionFind(cofacetidx);
            component_map[root].push_back(cofacetidx);
        }
    }

    std::vector<std::vector<size_t>> gradient_paths;
    gradient_paths.reserve(component_map.size());

    for (const auto& [root, path] : component_map)    //structured binding requires C++17
    {
        if (!path.empty())
        {
            gradient_paths.push_back(std::move(path));
        }
    }

    return gradient_paths;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::getGradientPathVertexSets(const MatchingContext& matching_context, const std::vector< std::vector<size_t> >& gradient_paths, const size_t dim)
{
    std::vector<std::unordered_set<size_t>> gradient_path_vertex_sets;
    gradient_path_vertex_sets.reserve(gradient_paths.size());

    auto npt = matching_context.graph.unodes + matching_context.graph.vnodes;

    for (const auto& grad_path: gradient_paths)
    {
        std::unordered_set<size_t> vertex_set;

        for (auto cofacetidx : grad_path)
        {
            auto bindex = matching_context.sorted_cofacets[cofacetidx].first;
            auto simplex_vertices = SimplexUtility::getSimplexVertices(matching_context.binomial_table, bindex, npt, dim);
            vertex_set.insert(simplex_vertices.begin(), simplex_vertices.end());
        }
        gradient_path_vertex_sets.push_back(std::move(vertex_set));
    }

    return gradient_path_vertex_sets;
}

template <typename DistMatType>
std::vector<std::unordered_set<size_t>> QuotientAndExpand<DistMatType>::trimVertexSets(std::vector<std::unordered_set<size_t>>& gradient_path_vertex_sets)
{
    std::sort(gradient_path_vertex_sets.begin(), gradient_path_vertex_sets.end(),
        [](const std::unordered_set<size_t>& lhs, const std::unordered_set<size_t>& rhs) 
        { return lhs.size() > rhs.size();});

    std::vector<std::unordered_set<size_t>> trimmed_vertex_sets;

    std::unordered_set<size_t> claimed_vertices;

    for (auto& vertex_set : gradient_path_vertex_sets)
    {
        bool overlap = false;

        for (const auto& vertex : vertex_set)
        {
            if (claimed_vertices.count(vertex))    
            {
                overlap = true;    //discard the vertex set (the smaller set) if it overlaps with claimed vertices
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

template<typename DistMatType>
std::vector<size_t> QuotientAndExpand<DistMatType>::getActiveVertexIndices(const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices)
{
    const size_t originalvtnum = dist_mat_.dist_mat.size();

    std::vector<bool> is_virtualized(originalvtnum, false);

    for (const auto& virtual_vt_set: virtual_vertex_indices)
    {
        for (const auto vt : virtual_vt_set) is_virtualized[vt] = true;
    }

    std::vector<size_t> active_vertices;
    active_vertices.reserve(originalvtnum);

    for (size_t i = 0; i < originalvtnum; ++i)
    {
        if (!is_virtualized[i]) active_vertices.push_back(i);
    }

    size_t initvirtualidx = originalvtnum;
    for (size_t i = 0; i < virtual_vertex_indices.size(); ++i)
    {
        //if (virtual_vt_set.empty()) continue;
        active_vertices.push_back(initvirtualidx);
        initvirtualidx++;
    }

    //active vertices indices are sorted in ascending order
    return active_vertices;
}

template<typename DistMatType>
double QuotientAndExpand<DistMatType>::computeVirtualDistance(const size_t i, const size_t j, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices)
{
    const size_t originalvtnum = dist_mat_.dist_mat.size();

    //passed in as i < j

    if (i < originalvtnum && j < originalvtnum) return dist_mat_.dist_mat[i][j];

    const std::unordered_set<size_t> temp_i_set = (i < originalvtnum) ? std::unordered_set<size_t>{i} : std::unordered_set<size_t>();
    const std::unordered_set<size_t> temp_j_set = (j < originalvtnum) ? std::unordered_set<size_t>{j} : std::unordered_set<size_t>();

    const std::unordered_set<size_t>& vtset_i = (i < originalvtnum) ? temp_i_set : virtual_vertex_indices[i - originalvtnum];
    const std::unordered_set<size_t>& vtset_j = (j < originalvtnum) ? temp_j_set : virtual_vertex_indices[j - originalvtnum];

    double mindist = std::numeric_limits<double>::max();
    for (const auto& vi : vtset_i)
    {
        for (const auto& vj : vtset_j)
        {
            //no intersection. do not need to check for i == j (mindist == 0)
            if (vi < vj)
            {
                mindist = std::min(mindist, dist_mat_.dist_mat[vi][vj]);
            }
            else
            {
                mindist = std::min(mindist, dist_mat_.dist_mat[vj][vi]);
            }
        }
    }

    return mindist;
}

template<typename DistMatType>
robin_hood::unordered_map<uint64_t, double> QuotientAndExpand<DistMatType>::getVirtualDistanceHashTable(const std::vector<size_t>& active_vertices, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices, int threadnum)
{
    omp_set_num_threads(threadnum);

    std::vector< robin_hood::unordered_map<uint64_t, double> > thread_hash_tables(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < active_vertices.size(); ++i)
    {
        size_t threadid = omp_get_thread_num();
        auto& thread_hash_table = thread_hash_tables[threadid];

        for (size_t j = i + 1; j < active_vertices.size(); ++j)
        {
            uint64_t key = (static_cast<uint64_t>(active_vertices[i]) << 32) | static_cast<uint64_t>(active_vertices[j]);
            double dist = computeVirtualDistance(active_vertices[i], active_vertices[j], virtual_vertex_indices);
            thread_hash_table.emplace(key, dist);
        }
    }

    robin_hood::unordered_map<uint64_t, double> virtual_distance_hash_table;
    for (const auto& thread_hash_table : thread_hash_tables)
    {
        virtual_distance_hash_table.insert(thread_hash_table.begin(), thread_hash_table.end());
    }

    return virtual_distance_hash_table;
}

template<typename DistMatType>
double QuotientAndExpand<DistMatType>::getVirtualDistance(size_t i, size_t j, const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table)
{
    if (i > j) std::swap(i, j);

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

template<typename DistMatType>
std::vector<std::pair<int64_t, double>> QuotientAndExpand<DistMatType>::getVirtualSortedEdge(const std::vector<size_t>& active_vertices,
                                                                 const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table, const double maxeps, int threadnum)
{
    omp_set_num_threads(threadnum);

    std::vector< std::vector< std::pair<int64_t, double> > > thread_edge_workspace(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < active_vertices.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_edges = thread_edge_workspace[threadid];

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
    for (const auto& thread_edges : thread_edge_workspace)
    {
        sorted_edges.insert(sorted_edge.end(), thread_edges.begin(), thread_edges.end());
    }

    SimplexUnitility::sortSimplexByWeightThenIndex(sorted_edges);

    return sorted_edges;
}

template<typename DistMatType>
robin_hood::unordered_map<int64_t, size_t> QuotientAndExpand<DistMatType>::getVirtualActiveEdgeIndexHashTable(const std::vector<std::pair<int64_t, double>>& sorted_virtual_edge, const size_t virtualvtnum)
{
    robin_hood::unordered_map<int64_t, size_t> active_edge_hash_table;
    active_edge_hash_table.reserve(sorted_virtual_edge.size());

    size_t npts = binomial_table_.size() - 1;    //original number of vertices + number of virtual vertices

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

template<typename DistMatType>
std::vector<std::pair<int64_t, double>> QuotientAndExpand<DistMatType>::getVirtualSortedCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex_list,
                                                                                                    const std::vector<size_t>& active_vertices, const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table,
                                                                                                    const size_t dim, const double maxeps, int threadnum)
{
    std::vector< std::vector< std::pair<size_t, double> > > thread_workspace(sorted_virtual_simplex.size());

    const size_t npts = binomial_table_.size() - 1;
    
    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_virtual_simplex_list.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];

        const auto& simplex_pair = sorted_virtual_simplex_list[i];
        const int64_t bindex = simplex_pair.first;
        const double weight = simplex_pair.second;

        auto simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        auto iter = std::find(active_vertices.begin(), active_vertices.end(), simplex_vertices.back());

        if (iter == active_vertices.end()) std::error("vertex not found in active vertex list");
        
        auto vtpos = std::distance(active_vertices.begin(), iter);

        for (size_t j = 0; j < vtpos; ++j)
        {
            const size_t covt = active_vertices[j];

            double newweight = 0.0;
            for (const auto& vt: simplex_vertices)
            {
                newweight = std::max(newweight, getVirtualDistance(covt, vt, virtual_distance_hash_table));
            }

            double cofacetweight = std::max(newweight, weight);

            if (cofacetweight < maxeps)
            {
                int64_t shiftedbindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vertices, 1);
                int64_t cofacetbindex = shiftedbindex + covt;
                thread_cofacet.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    std::vector<std::pair<int64_t, double>> sorted_cofacets;
    for (const auto& thread_cofacets : thread_workspace)
    {
        sorted_cofacets.insert(sorted_cofacets.end(), thread_cofacets.begin(), thread_cofacets.end());
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_cofacets);

    return sorted_cofacets;
}

template<typename DistMatType>
void QuotientAndExpand<DistMatType>::buildInterface(const std::vector<std::pair<int64_t, double>>& sorted_cofacet_list,
                                                    const robin_hood::unordered_map<int64_t, size_t>& active_facet_index_hash_table, const size_t dim)
{
    auto u = bi_graph.unodes
    auto v = bi_graph.vnodes;

    // std::cout<<u<<"    "<<v<<'\n';

    for (size_t i = 0; i < sorted_cofacet_list.size(); i++)
    {
        int64_t bindex = sorted_cofacet_list[i].first;
        std::vector<int64_t> facet_bindex_list = SimplexUtility::getFacetBinomialIndex(binomial_table_, bindex, dim);

        for (auto fbidx: facet_bindex_list)
        {
            auto mapit = active_facet_index_hash_table.find(fbidx);
            if (mapit != active_facet_index_hash_table.end())
            {
                size_t fidx = (*mapit).second;    //facet idx in facet bindex_weight_pair list
                bi_graph.adj_list[i].push_back(u + fidx);
                bi_graph.adj_list[u + fidx].push_back(i);
            }
        }
        std::sort(bi_graph.adj_list[i].begin(), bi_graph.adj_list[i].end(), std::greater<size_t>());
    }

    return;
}