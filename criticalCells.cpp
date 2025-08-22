#include "readInput.hpp"
#include "hopcroft_karp.hpp"

#include "SimplexUtility.hpp"
#include "SimplexEnumerator.hpp"
#include "MatchingContext.hpp"
#include "ApproximateMorseMatching.hpp"
#include "MaximumMorseMatching.hpp"

#include "criticalCells.hpp"

#include "omp.h"

#include <numeric>
#include <unordered_set>


#ifdef MPI_ENABLED
#include "mpi.h"
#endif

template class CritCells<VR, NormalDistMat>;
template class CritCells<Alpha, NormalDistMat>;

// template <>
// CritCells<VR, SparseDistMat>::CritCells(Eigen::SparseMatrix<double> &distMat)
// {
//     this->distMatrix = distMat;
// }

template <typename ComplexType, typename DistMatType>
CritCells<ComplexType, DistMatType>::CritCells(const std::string &filename):point_cloud_(readInput::readCSV(filename)), DistMatType(point_cloud_)
{}

template <typename ComplexType, typename DistMatType>
CritCells<ComplexType, DistMatType>::CritCells(std::vector<std::vector<double>> &distMat)
{
    this->distMatrix = distMat;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::run_Compute(int maxDim, int batch_size)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    ComplexType temp_complex;
    auto bins = dsimplices_batches(temp_complex, 1, 0);
    for (size_t dim = 2; dim <= maxDim; dim++)
    {
        std::clog << "Starting for dim " << dim << std::endl;
        auto batch_no = 1;
        ComplexType complex;
        while (complex.active)
        {
            std::clog << "Batch no " << batch_no++ << std::endl;
            auto batch_start_time = std::chrono::high_resolution_clock::now();
            auto weighted_simplicies = dsimplices_batches(complex, dim, batch_size);
            auto batch_end_time = std::chrono::high_resolution_clock::now();
            std::clog << "Batch time: " << std::chrono::duration<double>(batch_end_time - batch_start_time).count() << " seconds." << std::endl;

            auto match_start_time = std::chrono::high_resolution_clock::now();

            // Bin the batches
            binByWeights(weighted_simplicies, bins);
            for (auto &it : bins)
                it.second = dimMatching(it.second, dim, dim == maxDim);
            auto match_end_time = std::chrono::high_resolution_clock::now();
            std::clog << "Match time: " << std::chrono::duration<double>(match_end_time - match_start_time).count() << " seconds." << std::endl;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::clog << "Elapsed time: " << std::chrono::duration<double>(end_time - start_time).count() << " seconds." << std::endl;
    std::cout << bins << std::endl;
}

template <typename ComplexType, typename DistMatType>
std::map<double, std::vector<std::vector<int>>> CritCells<ComplexType, DistMatType>::binEdgeSimplexes() // Direct creation of edgebins to a map
{
    std::map<double, std::vector<std::vector<int>>> binned_edges;
    for (int i = 0; i < this->distMatrix.size() - 1; i++)
        for (int j = i + 1; j < this->distMatrix.size(); j++)
            binned_edges[this->distance(i, j)].push_back({i, j});
    return binned_edges;
}


template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins) // Merged higher dim feature to bins
{
    for (auto &[weight, simplexes] : weighted_simplicies)
        std::move(simplexes.begin(), simplexes.end(), std::back_inserter(bins[weight]));
    return;
}

template <typename ComplexType, typename DistMatType>
std::map<double, std::vector<std::vector<int>>> CritCells<ComplexType, DistMatType>::dsimplices_batches(ComplexType &complex, size_t dim, size_t batch_size) // Worker is invokation counter
{
    std::map<double, std::vector<std::vector<int>>> weighted_simplexes;
    if (!complex.is_initialized)
        complex.initialize(this->distMatrix.size(), dim + 1);
    long long counter = 0;
    do
    {
        weighted_simplexes[getSimplexWeight(complex.simplex)].push_back(complex.simplex);
        complex.next_simplex();
    } while (complex.active && (!batch_size || batch_size > ++counter));
    return weighted_simplexes;
}

template <typename ComplexType, typename DistMatType>
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final)
{
    std::vector<std::vector<int>> critCells, simps, cofaces;
    for (auto &simplex : simplexes)
    {
        if (simplex.size() == dim + 1)
            cofaces.emplace_back(std::move(simplex));
        else if (simplex.size() == dim)
            simps.emplace_back(std::move(simplex));
        else if (simplex.size() < dim)
            critCells.emplace_back(std::move(simplex));
    }
    simplexes.clear();
    if (simps.empty() || cofaces.empty())
    {
        std::move(simps.begin(), simps.end(), std::back_inserter(critCells));
        if (!final)
            std::move(cofaces.begin(), cofaces.end(), std::back_inserter(critCells));
        return critCells;
    }
    HKGraph csr_matrix(simps.size(), cofaces.size());
    for (std::size_t i = 0; i < simps.size(); ++i)
    {
        for (std::size_t j = 0; j < cofaces.size(); ++j)
        {
            if (std::includes(cofaces[j].begin(), cofaces[j].end(), simps[i].begin(), simps[i].end()))
                csr_matrix.addEdge(i + 1, j + 1); // Intial index 1
        }
    }
    auto res = csr_matrix.custom_hopcroftKarpAlgorithm();

    std::for_each(res.first.rbegin(), res.first.rend(), [&simps](auto index)
                  { simps.erase(std::next(simps.begin(), index - 1)); });
    std::move(simps.begin(), simps.end(), std::back_inserter(critCells));
    if (!final)
    {
        std::for_each(res.second.rbegin(), res.second.rend(), [&cofaces](auto index)
                      { cofaces.erase(std::next(cofaces.begin(), index - 1)); });
        std::move(cofaces.begin(), cofaces.end(), std::back_inserter(critCells));
    }
    return critCells;
}



template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::buildInterface(BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& cofacet_list, 
                                                        const robin_hood::unordered_map<int64_t, size_t>& active_facet_index, const std::vector<std::vector<int64_t>>& binomial_table, const size_t dim)
{
    // size_t npt = this->distMatrix.size();
    auto u = bi_graph.unodes
    auto v = bi_graph.vnodes;

    // std::cout<<u<<"    "<<v<<'\n';

    for (size_t i = 0; i < cofacet_list.size(); i++)
    {
        int64_t bindex = cofacet_list[i].first;
        std::vector<int64_t> facet_bindex = getFacetBinomialIndex(binomial_table, bindex, dim);

        for (auto fbidx: facet_bindex)
        {
            auto mapit = active_facet_index.find(fbidx);
            if (mapit != active_facet_index.end())
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



template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::runVRMorseTest(size_t maxdim, double maxeps, int threadnumber)
{
    size_t n = this->dist_matrix.size();
    auto binom_table = SimplexUtility::getBinomialTable(n, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(static_cast<const DistMatType&>(*this), binom_table);

    auto sorted_simplex = simplex_enumerator.getSortedVREdges(binom_table, maxeps);

    auto active_index_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binom_table, sorted_simplex, n);

    auto sorted_cofacet = simplex_enumerator.getSortedVRCofacets(binom_table, sorted_simplex, 1, maxeps, threadnumber);

    std::vector< std::vector< std::pair<double, double> > > persistent_pairs;
    persistent_pairs.push_back(std::vector< std::pair<double, double> >{std::make_pair(0.0, -1.0)});
    
    BipartiteGraph bi_graph(1, 1, 1);

    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        bi_graph.updateDimension(sorted_cofacet.size(), sorted_simplex.size());
        buildInterface(bi_graph, sorted_cofacet, active_index_hash_table, binom_table, dim);
        
        MatchingContext matching_context(bi_graph, binom_table, sorted_simplex, sorted_cofacet);
        
        MaximumMorseMatching morse_matching(threadnumber);

        auto critsimpnum = morse_matching.match(matching_context);

        std::cout << "critical simplex number: " << critsimpnum << std::endl;

        //print DAG edges at the top dim
        // if (dim == maxdim)
        // {
        //     std::ofstream out("dag_edges.csv");

        //     std::cout<<"DAG edges: \n";
        //     for (size_t i = bi_graph.u; i < (bi_graph.u + bi_graph.v); i++)
        //     {
        //         auto mate = bi_graph.match_list[i];
        //         if (mate != -1)
        //         {
        //             for (auto j : bi_graph.adj_list[i])
        //             {
        //                 auto jmate = bi_graph.match_list[j];
        //                 if (jmate != -1 && jmate != i)
        //                 {
        //                     std::cout<<"x = "<<j<<"  y = "<<mate<<'\n';
        //                     out<<j<<","<<mate<<'\n';
        //                 }
        //             }
        //         }
        //     }
        // }


        if (dim != maxdim)
        {
            active_index_hash_table = SimplexUtility::getActiveSimplexIndexHashTable(bi_graph.match_list, sorted_cofacet);

            sorted_simplex = simplex_enumerator.getSortedCofacets(binom_table, sorted_cofacet, dim, maxeps, threadnumber);

            std::swap(sorted_simplex, sorted_cofacet);
        }

    }

    return;
}



template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::runAlphaMorseTest(double maxeps, int threadnumber)
{
    // int threadnumber = 1;

    size_t maxdim = point_cloud_[0].size();

    typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> DT;
    typedef DT::Vertex_handle VT;

    DT delaunay_d(maxdim);

    std::vector<DT::Point> dt_points;
    std::unordered_map<VT, size_t> vertex_handle_index;
    vertex_handle_index.reserve(point_cloud_.size());

    for (size_t i = 0; i < point_cloud_.size(); i++)
    {
        dt_points.emplace_back(point_cloud_[i].begin(), point_cloud_[i].end());
        auto vhandle = delaunay_d.insert(dt_points[i]);
        vertex_handle_index[vhandle] = i;
    }

    dt_points.clear();

    //morse part
    size_t n = point_cloud_.size();
    auto binom_table = SimplexUtility::getBinomialTable(n, maxdim);

    SimplexEnumerator<DistMatType> simplex_enumerator(static_cast<const DistMatType&>(*this), binom_table);

    //get edges
    auto sorted_simplex = simplex_enumerator.getSortedAlphaCells(binom_table, vertex_handle_index, delaunay_d, 1, maxeps);
    //remove mst edges
    auto active_facet_hash_table = SimplexUtility::getActiveEdgeIndexHashTable(binom_table, sorted_simplex, n);

    //get triangles
    auto sorted_cofacet = simplex_enumerator.getSortedAlphaCells(binom_table, vertex_handle_index, delaunay_d, 2, maxeps);

    std::vector< std::vector< std::pair<double, double> > > persistent_pairs;
    persistent_pairs.push_back(std::vector< std::pair<double, double> >{std::make_pair(0.0, -1.0)});

    BipartiteGraph bi_graph(1, 1, 1);
    
    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        bi_graph.updateDimension(sorted_cofacet.size(), sorted_simplex.size());
        buildInterface(bi_graph, binom_table, sorted_cofacet, dim, active_facet_hash_table);

        MatchingContext matching_context(bi_graph, binom_table, sorted_simplex, sorted_cofacet);

        MaximumMorseMatching morse_matching(threadnumber);
        auto critsimpnum = morse_matching.match(matching_context);
        std::cout << "critical simplex number: " << critsimpnum << std::endl;

        if (dim != maxdim)
        {
            active_facet_hash_table = getActiveFacetIndexHashTable(bi_graph, sorted_cofacet);
            sorted_simplex = simplex_enumerator.getSortedAlphaCells(binom_table, vertex_handle_index, delaunay_d, dim + 1, maxeps);
            std::swap(sorted_simplex, sorted_cofacet);
        }
    }

    return;
}


