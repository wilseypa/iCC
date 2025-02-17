#include "distMat.cpp"
#include "readInput.hpp"
#include "hopcroft_karp.hpp"
#include "criticalCells.hpp"

#include "omp.h"

#include <numeric>

#ifdef MPI_ENABLED
#include "mpi.h"
#endif

template <>
CritCells<VR, SparseDistMat>::CritCells(Eigen::SparseMatrix<double> &distMat)
{
    this->distMatrix = distMat;
}

template <typename ComplexType, typename DistMatType>
CritCells<ComplexType, DistMatType>::CritCells(const std::string &fileName)
{
    auto inputData = readInput::readCSV(fileName);
    this->distMatrix = distMat(inputData);
}

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
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::getEdges(double maxweight)
{
    std::vector<std::vector<int>> edge_vec;
    for (int i = 0; i < this->distMatrix.size() - 1; i++)
    {
        for (int j = i + 1; j < this->distMatrix.size(); j++)
        {
            if (this->distMatrix[i][j] < maxweight)
                edge_vec.push_back(std::vector<int>{i, j});
        }
    }
    return edge_vec;
}

template <typename ComplexType, typename DistMatType>
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::getEdgesByWeightRange(std::vector<std::vector<int>> &sorted_edges, double minweight, double maxweight)
{
    auto findmin_lambda = [minweight, this](std::vector<int> &edge)
    { return (this->distMatrix[edge[0]][edge[1]] > minweight); };
    auto findmax_lambda = [maxweight, this](std::vector<int> &edge)
    { return (this->distMatrix[edge[0]][edge[1]] > maxweight); };

    std::vector<std::vector<int>>::iterator minit, maxit;
    if (minweight == 0)
    {
        minit = sorted_edges.begin();
    }
    else
    {
        minit = std::find_if(sorted_edges.begin(), sorted_edges.end(), findmin_lambda);
    }
    maxit = std::find_if(sorted_edges.begin(), sorted_edges.end(), findmax_lambda);

    std::vector<std::vector<int>> edge_bin(minit, maxit);
    return edge_bin;
}

template <typename ComplexType, typename DistMatType>
double CritCells<ComplexType, DistMatType>::getSimplexWeight(const std::vector<int> &simplex)
{
    double maxweight = 0;
    for (int i = 0; i < simplex.size() - 1; i++)
        for (int j = i + 1; j < simplex.size(); j++)
            maxweight = std::max(maxweight, this->distMatrix[simplex[i]][simplex[j]]);
    return maxweight;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::sortEdge(std::vector<std::vector<int>> &edge_bin)
{
    std::sort(edge_bin.begin(), edge_bin.end(), [this](const auto &lhs, const auto &rhs)
              { return getSimplexWeight(lhs) < getSimplexWeight(rhs); });
}

template <typename ComplexType, typename DistMatType>
int CritCells<ComplexType, DistMatType>::getMaxFacetIndex(const std::vector<int>& cofacet, const std::vector<std::vector<int>>& facet_bin)
{
    auto find_lamda = [&](const std::vector<int>& facet)
    { return std::includes(cofacet.begin(), cofacet.end(), facet.begin(), facet.end()); };

    auto riter = std::find_if(facet_bin.rbegin(), facet_bin.rend(), find_lamda);
    int index = std::distance(riter, facet_bin.rend());

    return index;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::sortSimplex(std::vector<std::vector<int>>& cofacet_bin, const std::vector<std::vector<int>>& facet_bin)
{
    auto sort_lambda = [this, &cofacet_bin](const auto& lhs, const auto& rhs)
    {
        double lhsweight = getSimplexWeight(lhs);
        double rhsweight = getSimplexWeight(rhs);

        if (lhsweight == rhsweight)    //reverse lexicographic order
        {
            auto lhsit = lhs.rbegin();
            auto rhsit = rhs.rbegin();
            while (lhsit != lhs.rend())
            {
                if (*lhsit != *rhsit) return *lhsit > *rhsit;
                ++lhsit;
                ++rhsit;
            }
            //duplicate simp. should never be here.
            return false;
        } 
        else
        {
            return lhsweight < rhsweight;
        }
    };

    std::sort(cofacet_bin.begin(), cofacet_bin.end(), sort_lambda);
}


template <typename ComplexType, typename DistMatType>
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::getCofacetBin(std::vector<std::vector<int>> &facet_bin, double maxweight, int threadnum)
{
    std::vector<std::vector<int>> cofacet_bin;
    size_t npts = this->distMatrix.size();

    omp_set_num_threads(threadnum);
    std::vector<std::vector<int>> thread_workspace(facet_bin.size(), std::vector<int>());

#pragma omp parallel for
    for (size_t i = 0; i < facet_bin.size(); i++)
    {
        double simpweight = getSimplexWeight(facet_bin[i]);

        int last_idx = *std::min_element(facet_bin[i].begin(), facet_bin[i].end());
        for (int j = 0; j < last_idx; j++)
        {
            auto maxidx = *std::max_element(facet_bin[i].begin(), facet_bin[i].end(), [this, j](auto first, auto second)
                                            { return this->distMatrix[j][first] < this->distMatrix[j][second]; });

            double weight = this->distMatrix[j][maxidx];

            if (weight < maxweight)
                thread_workspace[i].push_back(j);
        }
    }
    for (size_t i = 0; i < facet_bin.size(); i++)
    {
        for (auto &pt : thread_workspace[i])
        {
            // add pt first, keep the vertices in sorted order
            cofacet_bin.push_back(std::vector<int>{pt});
            std::copy(facet_bin[i].begin(), facet_bin[i].end(), std::back_inserter(cofacet_bin.back()));
        }
    }

    sortSimplex(cofacet_bin, facet_bin);

    return cofacet_bin;
}


template <typename ComplexType, typename DistMatType>
int CritCells<ComplexType, DistMatType>::findRoot(std::vector<int> &parent_idx, int x)
{
    while (parent_idx[x] != x)
    {
        parent_idx[x] = parent_idx[parent_idx[x]];
        x = parent_idx[x];
    }
    return x;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::setUnion(std::vector<int> &parent_idx, int x, int y)
{
    int p = findRoot(parent_idx, x);
    int q = findRoot(parent_idx, y);
    parent_idx[p] = parent_idx[q];
}

template <typename ComplexType, typename DistMatType>
std::vector<int> CritCells<ComplexType, DistMatType>::getMSTEdgeIndices(std::vector<std::vector<int>> &sorted_edges)
{
    int n = this->distMatrix.size();

    std::vector<int> parent_idx(n);
    std::iota(parent_idx.begin(), parent_idx.end(), 0);

    std::vector<int> mst_edge_index;

    int x, y;
    int m = sorted_edges.size();
    for (int i = 0; i < m; i++)
    {
        x = sorted_edges[i][0];
        y = sorted_edges[i][1];
        if (findRoot(parent_idx, x) != findRoot(parent_idx, y))
        {
            setUnion(parent_idx, x, y);
            mst_edge_index.push_back(i);
        }
    }

    return mst_edge_index;
}

// template <typename ComplexType, typename DistMatType>
// std::vector<std::vector<double>> CritCells<ComplexType, DistMatType>::run_MorseMatch(int maxdimension, double mineps, double maxeps, int threadnumber)
// {
//     int threadnum = threadnumber;

//     std::vector<std::vector<int>> simplex_bin = getEdges(maxeps);
//     sortEdge(simplex_bin);

//     std::vector<int> mst_edge_index = getMSTEdgeIndices(simplex_bin);

//     std::vector<int> dim_active_index;
//     for (int i = 0; i < simplex_bin.size(); i++)
//     {
//         if (std::find(mst_edge_index.begin(), mst_edge_index.end(), i) == mst_edge_index.end())
//             dim_active_index.push_back(i);
//     }

//     std::vector<std::vector<double>> critical_weight;
//     critical_weight.push_back(std::vector<double>{0});

//     int initleftdeg = 3; // cofacet of edge, tri

//     std::vector<std::vector<int>> cofacet_bin = getCofacetBin(simplex_bin, maxeps, threadnum);


//     Bi_Graph_Match bi_graph(cofacet_bin.size(), simplex_bin.size(), initleftdeg, threadnum);

//     for (int dim = 2; dim <= maxdimension; dim++)
//     {
//         // std::cout<<"dim active idx size = "<<dim_active_index.size()<<'\n';
//         bi_graph.buildInterface(cofacet_bin, 0, cofacet_bin.size(), simplex_bin, simplex_bin.size(), dim_active_index);

//         std::cout<<"cofacet size = "<<cofacet_bin.size()<<"  facet size = "<<simplex_bin.size()<<'\n';
//         for (int i = 50; i < 65; i++)
//         {
//             for (auto j: bi_graph.adj_list[i]) std::cout<<j<<"  ";
//             std::cout<<'\n';
//         }

//         std::cout <<"dim = " <<dim<< "  simp bin active idx size = " << dim_active_index.size() << "  simp bin size = " << simplex_bin.size();
//         std::cout <<"  cofacet size = " << cofacet_bin.size() << '\n';

//         // bi_graph.parallelKarpSipserInit();

//         bi_graph.parallelMaxFacetInit(0, cofacet_bin.size(), 0, simplex_bin.size());

//         // if (dim)
//         // {
//         //     std::vector<int> unmatched_facet;
//         //     std::vector<double> facet_weight;
//         //     for (int i: dim_active_index)
//         //     {
//         //         int vidx = i + cofacet_bin.size();
//         //         if (bi_graph.match_list[vidx] < 0) 
//         //         {
//         //             unmatched_facet.push_back(i);
//         //             facet_weight.push_back(getSimplexWeight(simplex_bin[i]));
//         //         }
//         //     }
//         //     std::cout<<"dim = "<<dim<<" unmatched facet index after init = ";
//         //     for(int j: unmatched_facet) std::cout<<j<<"  ";
//         //     std::cout<<'\n';
//         //     for(double w: facet_weight) std::cout<<w<<"  ";
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 2)
//         // {   
//         //     std::cout<<"before match/rm  ";
//         //     std::vector<int> target_idx{71, 72, 82, 184, 185};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkCofacetByIndex(cofacet_bin, simplex_bin, target_idx);
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 2)
//         // {   
//         //     std::cout<<"before match/rm  ";
//         //     std::vector<int> target_idx{52, 57, 65};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkSimplexByIndex(cofacet_bin, simplex_bin, target_idx);

//         //     for(auto i: target_idx)
//         //     {
//         //         std::cout<<getSimplexWeight(simplex_bin[i])<<"  ";
//         //     }
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 3)
//         // {   
//         //     std::cout<<"before match/rm  ";
//         //     std::vector<int> target_idx{241, 362, 363};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkCofacetByIndex(cofacet_bin, simplex_bin, target_idx);
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 3)
//         // {   
//         //     std::cout<<"before match/rm  ";
//         //     std::vector<int> target_idx{184, 185, 227};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkSimplexByIndex(cofacet_bin, simplex_bin, target_idx);

//         //     for(auto i: target_idx)
//         //     {
//         //         std::cout<<getSimplexWeight(simplex_bin[i])<<"  ";
//         //     }
//         //     std::cout<<'\n';
//         // }

//         // bi_graph.parallelFacetDFSMatch();

//         bi_graph.serialCofacetDFSMatch();
        
//         // for (int i = 0; i < cofacet_bin.size(); i++)
//         // {
//         //     int imate = bi_graph.match_list[i];
//         //     if (imate != - 1 && bi_graph.match_list[imate] != i)
//         //     {
//         //         std::cout<<"graph error at dim = "<<dim<<"  error uidx = "<<i<<'\n';
//         //     }
//         // }

//         // for (int i = cofacet_bin.size(); i < cofacet_bin.size() + simplex_bin.size(); i++)
//         // {
//         //     int imate = bi_graph.match_list[i];
//         //     if (imate != - 1 && bi_graph.match_list[imate] != i)
//         //     {
//         //         std::cout<<"graph error at dim = "<<dim<<"  error vidx = "<<i<<'\n';
//         //     }
//         // }

//         // int reverted = bi_graph.dfsCycleRemoval();

//         // std::cout << "dim = " << dim << "  reverted = " << reverted << '\n';

//         // if (dim == 2)
//         // {   
//         //     std::cout<<"after match/rm  ";
//         //     std::vector<int> target_idx{71, 72, 82, 184, 185};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkCofacetByIndex(cofacet_bin, simplex_bin, target_idx);
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 2)
//         // {   
//         //     std::cout<<"after match/rm  ";
//         //     std::vector<int> target_idx{52, 57, 65};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkSimplexByIndex(cofacet_bin, simplex_bin, target_idx);
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 3)
//         // {   
//         //     std::cout<<"after match/rm  ";
//         //     std::vector<int> target_idx{241, 361, 362, 363, 364, 365, 366, 367};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkCofacetByIndex(cofacet_bin, simplex_bin, target_idx);

//         //     for(auto i: target_idx)
//         //     {
//         //         std::cout<<getSimplexWeight(cofacet_bin[i])<<"  ";
//         //     }
//         //     std::cout<<'\n';
//         // }

//         // if (dim == 3)
//         // {   
//         //     std::cout<<"after match/rm  ";
//         //     std::vector<int> target_idx{184, 185, 227, 265, 266, 267, 268};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkSimplexByIndex(cofacet_bin, simplex_bin, target_idx);
//         //     for(auto i: target_idx)
//         //     {
//         //         std::cout<<getSimplexWeight(simplex_bin[i])<<"  ";
//         //     }
//         //     std::cout<<'\n';
//         // }


//         std::vector<int> dim_critical_index = bi_graph.getCriticalIndex(dim_active_index, simplex_bin.size());
//         std::cout << "dim = " << dim << "  dim - 1 critical idx size = " << dim_critical_index.size() << '\n';


//         dim_active_index = bi_graph.getActiveIndex();

//         std::cout << "cofact active idx size = " << dim_active_index.size() << "  cofacet size = " << cofacet_bin.size() << '\n';
        

//         // std for each
//         std::vector<double> dim_critical_weight;
//         for (auto i : dim_critical_index)
//         {
//             double weight = getSimplexWeight(simplex_bin[i]);
//             dim_critical_weight.push_back(weight);
//         }
//         critical_weight.push_back(std::move(dim_critical_weight));

//         if (dim < maxdimension)
//         {
//             simplex_bin = getCofacetBin(cofacet_bin, maxeps, threadnum);
//             std::swap(simplex_bin, cofacet_bin);
//             bi_graph.updateDimension(cofacet_bin.size(), simplex_bin.size());
//         }

        
//     }

//     return critical_weight;
// }

template <typename ComplexType, typename DistMatType>
std::vector<double> CritCells<ComplexType, DistMatType>::getApproxDeathWeight(std::vector<std::vector<int>>& facet_bin, std::vector<std::set<int>>& backward_facet_index, double maxeps)
{
    std::vector<double> death_weight;

    size_t npts = this->distMatrix.size();

    omp_set_num_threads(1);

    std::vector<std::vector<int>> thread_workspace(facet_bin.size(), std::vector<int>());

    for (auto& facet_index_set: backward_facet_index)
    {
        std::vector<double> thread_workspace(facet_index_set.size(), maxeps);

        //***********************fix later***********************************//
        double minweight = 9999;

#pragma omp parallel for
        for (size_t i = 0; i < facet_index_set.size(); i++)
        {
            auto idx = *std::next(facet_index_set.begin(), i);

            for (size_t j = 0; j < npts; j++)
            {
                if (std::find(facet_bin[idx].begin(), facet_bin[idx].end(), j) != facet_bin[idx].end()) continue;
                
                auto weight_lambda = [this, j](auto first, auto second)
                                    { double firstweight = (first < j) ? this->distMatrix[first][j] : this->distMatrix[j][first];
                                      double secondweight = (second < j) ? this->distMatrix[second][j] : this->distMatrix[j][second];
                                      return firstweight < secondweight; };

                auto maxidx = *std::max_element(facet_bin[idx].begin(), facet_bin[idx].end(), weight_lambda);

                double cofacetweight = (maxidx < j) ? this->distMatrix[maxidx][j] : this->distMatrix[j][maxidx];

                if (cofacetweight > maxeps && cofacetweight < minweight) minweight = cofacetweight;
            }
            thread_workspace[i] = minweight;
        }

        double dweight = *std::min_element(thread_workspace.begin(), thread_workspace.end());
        death_weight.push_back(dweight);
    }
    return death_weight;
}


// template <typename ComplexType, typename DistMatType>
// std::vector< std::vector< std::pair<double, double> > > CritCells<ComplexType, DistMatType>::run_MorseMatchPersistence(int maxdimension, double mineps, double maxeps)
// {
//     int threadnum = 1;

//     std::vector<std::vector<int>> simplex_bin = getEdges(maxeps);
//     sortEdge(simplex_bin);

//     std::vector<int> mst_edge_index = getMSTEdgeIndices(simplex_bin);

//     std::vector<int> dim_active_index;
//     for (int i = 0; i < simplex_bin.size(); i++)
//     {
//         if (std::find(mst_edge_index.begin(), mst_edge_index.end(), i) == mst_edge_index.end())
//             dim_active_index.push_back(i);
//     }

//     // std::vector<std::vector<double>> critical_weight;
//     // critical_weight.push_back(std::vector<double>{0});

//     std::vector< std::vector< std::pair<double, double> > > persistent_pairs;
//     //H0.
//     persistent_pairs.push_back(std::vector< std::pair<double, double> >{std::make_pair(0, -1)});

//     int initleftdeg = 3; // cofacet of edge, tri

//     std::vector<std::vector<int>> cofacet_bin = getCofacetBin(simplex_bin, maxeps, threadnum);

//     Bi_Graph_Match bi_graph(cofacet_bin.size(), simplex_bin.size(), initleftdeg, threadnum);

//     for (int dim = 2; dim <= maxdimension; dim++)
//     {
//         // std::cout<<"dim active idx size = "<<dim_active_index.size()<<'\n';
//         bi_graph.buildInterface(cofacet_bin, 0, cofacet_bin.size(), simplex_bin, simplex_bin.size(), dim_active_index);

//         std::cout <<"dim = " <<dim<< "  simp bin active idx size = " << dim_active_index.size() << "  simp bin size = " << simplex_bin.size();
//         std::cout <<"  cofacet size = " << cofacet_bin.size() << '\n';

//         bi_graph.parallelMaxFacetInit(0, cofacet_bin.size(), 0, simplex_bin.size());

//         // bi_graph.parallelFacetDFSMatch();

//         // if (dim == 2)
//         // {   
//         //     std::cout<<"before match/rm  ";
//         //     std::vector<int> target_idx{37, 52, 77};
//         //     std::cout<<"dim = "<<dim<<'\n';
//         //     bi_graph.checkSimplexByIndex(cofacet_bin, simplex_bin, target_idx);

//         //     for(auto i: target_idx)
//         //     {
//         //         std::cout<<getSimplexWeight(simplex_bin[i])<<"  ";
//         //     }
//         //     std::cout<<'\n';
//         // }

//         bi_graph.serialCofacetDFSMatch();
         
//         // int reverted = bi_graph.dfsCycleRemoval();

//         // std::cout << "dim = " << dim << "  reverted = " << reverted << '\n';

//         std::vector<int> dim_critical_index = bi_graph.getCriticalIndex(dim_active_index, simplex_bin.size());

//         std::vector<std::set<int>> backward_facet = bi_graph.getBackwardSingleFacetIndex(dim_critical_index);

//         std::vector<double> approx_death_weight = getApproxDeathWeight(simplex_bin, backward_facet, maxeps);

//         for (auto w: approx_death_weight) std::cout<<w<<"   ";
//         std::cout<<'\n';

//         std::vector< std::pair<double, double> > dim_persis_pair;
//         for (size_t i = 0; i < dim_critical_index.size(); i++)
//         {
//             int idx = dim_critical_index[i];
//             double weight = getSimplexWeight(simplex_bin[idx]);
//             dim_persis_pair.push_back(std::make_pair(weight, approx_death_weight[i]));
//         }
//         persistent_pairs.push_back(std::move(dim_persis_pair));
        
//         dim_active_index = bi_graph.getActiveIndex();

//         if (dim < maxdimension)
//         {
//             simplex_bin = getCofacetBin(cofacet_bin, maxeps, threadnum);
//             std::swap(simplex_bin, cofacet_bin);
//             bi_graph.updateDimension(cofacet_bin.size(), simplex_bin.size());
//         }
//     } 
    
//     return persistent_pairs;
// }


template <typename ComplexType, typename DistMatType>
double CritCells<ComplexType, DistMatType>::getMinimumEpsilon(std::vector<int>& mst_edge_index, std::vector<std::vector<int>>& sorted_edges, double stepsize)
{
    std::vector<int>::iterator it;
    it = std::max_element(mst_edge_index.begin(), mst_edge_index.end());
    double weight = getSimplexWeight(sorted_edges[*it]);
    return (weight + stepsize);
}

template <typename ComplexType, typename DistMatType>
std::vector<int> CritCells<ComplexType, DistMatType>::getStepwiseIndex(std::vector<std::vector<int>>& simplex_bin, double mineps, double maxeps, double stepsize)
{
    std::vector<int> simplex_index{0};
    
    double eps = mineps;

    while (eps < maxeps)
    {
        auto weight_lambda = [eps, this](std::vector<int>& simplex) { return (getSimplexWeight(simplex) > eps); };
        std::vector<std::vector<int>>::iterator it = std::find_if(simplex_bin.begin(), simplex_bin.end(), weight_lambda);
        int index = std::distance(simplex_bin.begin(), it);
        simplex_index.push_back(index);
        eps += stepsize;
    }
    
    //last simp index
    simplex_index.push_back(simplex_bin.size());
    
    return simplex_index;
}


// template <typename ComplexType, typename DistMatType>
// std::vector< std::vector< std::pair<double, double> > > CritCells<ComplexType, DistMatType>::run_MorseMatchStepwise(int maxdimension, double maxeps, double stepsize)
// {
//     int threadnum = 1;

//     std::vector<std::vector<int>> simplex_bin = getEdges(maxeps);
//     sortEdge(simplex_bin);

//     std::vector<int> mst_edge_index = getMSTEdgeIndices(simplex_bin);

//     std::vector<int> dim_active_index;
//     for (int i = 0; i < simplex_bin.size(); i++)    //remove mst edges from active edges
//     {
//         if (std::find(mst_edge_index.begin(), mst_edge_index.end(), i) == mst_edge_index.end())
//             dim_active_index.push_back(i);
//     }
    
//     //mineps = mst_max + stepsize
//     double mineps = getMinimumEpsilon(mst_edge_index, simplex_bin, stepsize);

//     std::vector<std::vector<int>> cofacet_bin = getCofacetBin(simplex_bin, maxeps, threadnum);

//     int initleftdeg = 3; // cofacet of edge, tri

//     Bi_Graph_Match bi_graph(cofacet_bin.size(), simplex_bin.size(), initleftdeg, threadnum);    
    
//     std::vector< std::vector< std::pair<double, double> > > persistent_pairs;
//     //H0.
//     persistent_pairs.push_back(std::vector< std::pair<double, double> >{std::make_pair(0, maxeps)});

//     for (int dim = 2; dim <= maxdimension; dim++)
//     {
//         //simplex/cofacet index for each step
//         std::vector<int> simplex_index = getStepwiseIndex(simplex_bin, mineps, maxeps, stepsize);
//         std::vector<int> cofacet_index = getStepwiseIndex(cofacet_bin, mineps, maxeps, stepsize);
 
        
//         for(auto s: simplex_index) std::cout<<s<<"  ";
//         std::cout<<'\n';
//         for(auto s: cofacet_index) std::cout<<s<<"  ";
//         std::cout<<"\n\n";
        
//         if (true)
//         {
//             for(auto it = simplex_index.begin(); it != simplex_index.end() - 1; it++) std::cout<<getSimplexWeight(simplex_bin[*it])<<"  ";
//             std::cout<<'\n';
//             for(auto it = cofacet_index.begin(); it != cofacet_index.end() - 1; it++) std::cout<<getSimplexWeight(cofacet_bin[*it])<<"  ";
//             std::cout<<'\n';
//         }

//         std::vector<int> critical_index_before, critical_index_after;
//         std::vector< std::pair<double, double> > dim_persistent_pairs;

//         int stepnum = simplex_index.size() - 1;

//         for (int m = 1; m <= stepnum; m++)
//         {
//             bi_graph.buildInterface(cofacet_bin, cofacet_index[m - 1], cofacet_index[m], simplex_bin, simplex_index[m], dim_active_index);

//             bi_graph.parallelKarpSipserInit();
//             // bi_graph.parallelDFSMatch();
//             int reverted = bi_graph.serialCycleRemoval();
            
//             std::cout<<"dim = "<<dim<<"  step = "<<m<<" reverted = "<<reverted<<'\n';

//             critical_index_after = bi_graph.getCriticalIndex(dim_active_index, simplex_index[m]);  //double check

//             for(auto& i: critical_index_before)
//             {
//                 if (std::find(critical_index_after.begin(), critical_index_after.end(), i) == critical_index_after.end())
//                 {
//                     double birth = getSimplexWeight(simplex_bin[i]);

//                     double eps = mineps + (m - 1) * stepsize;    //eps of current step
//                     double death = (eps < maxeps) ? eps : maxeps;

//                     std::cout<<"dim = "<<dim<<"  birth = "<<birth<<"  death = "<<death<<'\n';

//                     dim_persistent_pairs.push_back(std::make_pair(birth, death));
//                 }
//             }

//             std::swap(critical_index_before, critical_index_after);

//             if (m == stepnum)
//             {
//                 for(auto& i: critical_index_before)
//                 {
//                     if (i >= simplex_index[m - 1])
//                     {
//                         double birth = getSimplexWeight(simplex_bin[i]);

//                         double eps = mineps + m * stepsize;    //eps of current step
//                         double death = (eps < maxeps) ? eps : maxeps;

//                         std::cout<<"dim = "<<dim<<"  birth = "<<birth<<"  death = "<<death<<'\n';

//                         dim_persistent_pairs.push_back(std::make_pair(birth, death));
//                     }
//                 }
//             }
//         }

//         persistent_pairs.push_back(std::move(dim_persistent_pairs));

//         dim_active_index = bi_graph.getActiveIndex();

//         std::cout<<"dim = "<<dim<<"  cofacet size = "<<cofacet_bin.size()<<"  active cofacet = "<<dim_active_index.size()<<'\n';

//         if (dim < maxdimension)
//         {
//             simplex_bin = getCofacetBin(cofacet_bin, maxeps, threadnum);
//             std::swap(simplex_bin, cofacet_bin);
//             bi_graph.updateDimension(cofacet_bin.size(), simplex_bin.size());
//         }
//     }

//     return persistent_pairs;
// }



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


//rewrite
template <typename ComplexType, typename DistMatType>
std::vector<std::vector<int64_t>> CritCells<ComplexType, DistMatType>::getBinomialTable(const size_t vtnum, const size_t maxdim)
{   
    //binom table by rows. same as lhf
    std::vector<std::vector<int64_t>> binom_table(vtnum + 1, std::vector<int64_t>(maxdim + 1 + 1, 0));

    binom_table[0][0] = 1;
    for(size_t i = 1; i <= vtnum; i++)
    {
        binom_table[i][0] = 1;
        for (size_t j = 1; j < (maxdim + 2); j++)
        {
            binom_table[i][j] = binom_table[i - 1][j - 1] + binom_table[i - 1][j];
            if (binom_table[i][j] < 0) throw std::overflow_error("Binomial Overflow Error!");
        }
    }

    return binom_table;
}

template <typename ComplexType, typename DistMatType>
int64_t CritCells<ComplexType, DistMatType>::getBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<size_t>& simplex_vt, const size_t shift)    //shift is used for cofacet enum
{
    int64_t bindex = 0;
    size_t k = simplex_vt.size();

    //simplex_vt is sorted in descending order
    for (auto it = simplex_vt.begin(); it != simplex_vt.end(); it++)
    {
        bindex += binomial_table[*it][k + shift];
        if (bindex < 0) throw std::overflow_error("Binomial overflow in get simplex binomial index.");
        k--;
    }

    return bindex;
}

template <typename ComplexType, typename DistMatType>
int64_t CritCells<ComplexType, DistMatType>::getEdgeBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, size_t j, size_t i)
{
    //assume j > i
    if (j < i) std::swap(i, j);

    return (binomial_table[j][2] + i);
}


template <typename ComplexType, typename DistMatType>
size_t CritCells<ComplexType, DistMatType>::getSimplexMaxVertex(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t bindex, size_t vtnum, const size_t dim)
{   
    //top = vtnum - 1, bottom = dim
    size_t top = vtnum - 1;
    if (binomial_table[top][dim + 1] > bindex)
    {
        size_t span = top - dim;
        while (span > 0)
        {
            size_t half = (span >> 1);
            size_t mid = top - half;
            if (binomial_table[mid][dim + 1] > bindex)
            {
                top = mid - 1;
                span -= (half + 1);
            }
            else span = half;
        }
    }
    return top;
}

template <typename ComplexType, typename DistMatType>
std::vector<size_t> CritCells<ComplexType, DistMatType>::getSimplexVertices(const std::vector<std::vector<int64_t>>& binomial_table, int64_t bindex, size_t vtnum, const size_t dim)
{
    //max vertex index = vtnum - 1
    std::vector<size_t> simplex_vt;
    simplex_vt.reserve(dim + 1);

    for (size_t k = dim + 1; k > 0; k--)  //size_t is unsigned. cannot use it to check against negative number
    {
        vtnum = getSimplexMaxVertex(binomial_table, bindex, vtnum, k - 1);

        bindex -= binomial_table[vtnum][k];

        simplex_vt.push_back(vtnum);
    }
    return simplex_vt;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::sortSimplexByWeightThenIndex(std::vector<std::pair<int64_t, double>>& simplex_list)
{
    auto sort_lambda = [](const auto& lhs, const auto& rhs){ return (lhs.second < rhs.second) || ((lhs.second == rhs.second) && (lhs.first < rhs.first)); };

    std::sort(simplex_list.begin(), simplex_list.end(), sort_lambda);

    return;
}

template <typename ComplexType, typename DistMatType>
std::vector<std::pair<int64_t, double>> CritCells<ComplexType, DistMatType>::getSortedEdge(const std::vector<std::vector<int64_t>>& binomial_table, const double maxeps)
{
    std::vector<std::pair<int64_t, double>> sorted_edge;

    size_t npt = this->distMatrix.size();

    for (size_t i = 0; i < npt - 1; i++)
    {
        for (size_t j = i + 1; j < npt; j++)
        {
            double weight = this->distMatrix[i][j];
            if (weight < maxeps)
            {
                int64_t bindex = getEdgeBinomialIndex(binomial_table, j, i);
                sorted_edge.push_back(std::make_pair(bindex, weight));
            }
        }
    }

    sortSimplexByWeightThenIndex(sorted_edge);

    return sorted_edge;
}

template <typename ComplexType, typename DistMatType>
size_t CritCells<ComplexType, DistMatType>::mstFindRoot(std::vector<size_t>& parent_idx, size_t x)
{
    while (parent_idx[x] != x)
    {
        parent_idx[x] = parent_idx[parent_idx[x]];
        x = parent_idx[x];
    }

    return x;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::mstSetUnion(std::vector<size_t>& parent_idx, size_t x, size_t y)
{
    size_t p = mstFindRoot(parent_idx, x);
    size_t q = mstFindRoot(parent_idx, y);
    parent_idx[p] = parent_idx[q];
    return;
}

template <typename ComplexType, typename DistMatType>
robin_hood::unordered_map<int64_t, size_t> CritCells<ComplexType, DistMatType>::getActiveEdgeIndexHashTable(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<std::pair<int64_t, double>>& sorted_edge)
{
    robin_hood::unordered_map<int64_t, size_t> active_edge_index_hash;
    active_edge_index_hash.reserve(sorted_edge.size());
    
    size_t npt = this->distMatrix.size();

    std::vector<size_t> parent_idx(npt);
    std::iota(parent_idx.begin(), parent_idx.end(), 0);

    size_t x, y;
    int64_t bindex;
    size_t m = sorted_edge.size();
    std::vector<size_t> edge_vt;
    for (auto i = 0; i < m; i++)
    {
        bindex = sorted_edge[i].first;
        edge_vt = getSimplexVertices(binomial_table, bindex, npt, 1);
        x = edge_vt[0];
        y = edge_vt[1];
        if (mstFindRoot(parent_idx, x) != mstFindRoot(parent_idx, y))
        {
            mstSetUnion(parent_idx, x, y);
            continue;
        }
        active_edge_index_hash.insert({bindex, i});    //need to use {} initializer list
    }

    return active_edge_index_hash;
}

template <typename ComplexType, typename DistMatType>
robin_hood::unordered_map<int64_t, size_t> CritCells<ComplexType, DistMatType>::getActiveFacetIndexHashTable(const std::vector<std::pair<int64_t, double>>& facet_list, const std::vector<size_t>& active_facet_index)
{
    robin_hood::unordered_map<int64_t, size_t> active_facet_index_hash;
    active_facet_index_hash.reserve(facet_list.size());

    for (auto i : active_facet_index)
    {
        int64_t bindex = facet_list[i].first;
        active_facet_index_hash.insert({bindex, i});
    }

    return active_facet_index_hash;
}

template <typename ComplexType, typename DistMatType>
std::vector<std::pair<int64_t, double>> CritCells<ComplexType, DistMatType>::getSortedCofacetList(const std::vector<std::vector<int64_t>>& binomial_table, 
                                                                                                  const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum)
{
    //dim == simplex dimension == cofacet dimension - 1
    std::vector<std::pair<int64_t, double>> cofacet_list;

    size_t npt = this->distMatrix.size();

    std::vector< std::vector< std::pair<size_t, double> > > thread_workspace(sorted_simplex.size(), std::vector<std::pair<size_t, double>>());

    omp_set_num_threads(threadnum);

#pragma omp parallel for
    for (size_t i = 0; i < sorted_simplex.size(); i++)
    {
        int64_t bindex = sorted_simplex[i].first;
        double originalweight = sorted_simplex[i].second;
        std::vector<size_t> simplex_vt = getSimplexVertices(binomial_table, bindex, npt, dim);

        //cofacet of {i, j} = {i, j, ...}. i > j
        // auto minidx = *std::min_element(simplex_vt.begin(), simplex_vt.end());
        auto minidx = simplex_vt.back();
        for (size_t j = 0; j < minidx; j++)
        {
            auto idx = *std::max_element(simplex_vt.begin(), simplex_vt.end(), [this, j](auto first, auto second)
                                                                                { return this->distMatrix[j][first] < this->distMatrix[j][second]; });
            double weight = this->distMatrix[j][idx];    //new facet weight
            double cofacetweight = (weight > originalweight) ? weight : originalweight; 

            if (cofacetweight < maxeps) thread_workspace[i].push_back({j, cofacetweight});
        }
    }

    for (size_t i = 0; i < sorted_simplex.size(); i++)
    {
        int64_t bindex = sorted_simplex[i].first;
        std::vector<size_t> simplex_vt = getSimplexVertices(binomial_table, bindex, npt, dim);
        //left shift bindex. left shift one position for each simplex vt
        bindex = getBinomialIndex(binomial_table, simplex_vt, 1);
        for (auto& idx_weight_pair: thread_workspace[i])
        {
            auto j = idx_weight_pair.first;
            auto w = idx_weight_pair.second;
            //push new bindex = shifted bin
            cofacet_list.push_back({bindex + j, w});
        }
    }

    sortSimplexByWeightThenIndex(cofacet_list);

    return cofacet_list;
}

template <typename ComplexType, typename DistMatType>
std::vector<int64_t> CritCells<ComplexType, DistMatType>::getFacetBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t binomindex, const size_t dim)
{
    //{i, j, k} -> {j, k}, {i, k}, {i, j}. i > j > k
    std::vector<int64_t> facet_bindex;
    facet_bindex.reserve(dim + 1);

    size_t npt = this->distMatrix.size();

    std::vector<size_t> simplex_vt = getSimplexVertices(binomial_table, binomindex, npt, dim);

    int64_t above = 0;
    int64_t below = binomindex;
    size_t k = dim;

    for (auto i = 0; i < simplex_vt.size(); i++)
    {
        size_t vt = simplex_vt[i];
        below -= binomial_table[vt][k + 1];

        facet_bindex.push_back(above + below);

        above += binomial_table[vt][k];

        k--;
    }

    return facet_bindex;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::buildInterface(Bi_Graph_Match& bi_graph, const std::vector<std::vector<int64_t>>& binomial_table, 
                                                         const std::vector<std::pair<int64_t, double>>& cofacet_list, const size_t dim, const robin_hood::unordered_map<int64_t, size_t>& active_facet_index)
{
    size_t npt = this->distMatrix.size();
    size_t u = bi_graph.u;
    size_t v = bi_graph.v;

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
void CritCells<ComplexType, DistMatType>::runTest(size_t maxdim, double maxeps)
{
    size_t n = this->distMatrix.size();
    auto binom_table = getBinomialTable(n, maxdim);

    auto sorted_simplex = getSortedEdge(binom_table, maxeps);

    auto active_index_hash_table = getActiveEdgeIndexHashTable(binom_table, sorted_simplex);
    std::vector<size_t> dim_active_index;
    for(auto& pair: active_index_hash_table) dim_active_index.push_back(pair.second);

    auto sorted_cofacet = getSortedCofacetList(binom_table, sorted_simplex, 1, maxeps, 1);

    std::vector<std::vector<double>> critical_weight;
    critical_weight.push_back(std::vector<double>{0});
    
    Bi_Graph_Match bi_graph(1, 1, 1);

    for (size_t dim = 2; dim <= maxdim; dim++)
    {
        bi_graph.updateDimension(sorted_cofacet.size(), sorted_simplex.size());
        buildInterface(bi_graph, binom_table, sorted_cofacet, dim, active_index_hash_table);

        bi_graph.parallelMaxFacetInit(0, sorted_cofacet.size(), 0, sorted_simplex.size(), 1);

        // bi_graph.parallelKarpSipserInit(1);

        bi_graph.parallelFacetDFSMatch(1);
        // bi_graph.serialCofacetDFSMatch();

        // std::cout<<"check graph dim = "<<dim<<"  cofacet size = "<<sorted_cofacet.size()<<"  facet size = "<<sorted_simplex.size()<<'\n';
        // std::vector<size_t> cof_idx = {269, 270, 271, 272, 273, 274};
        // std::vector<size_t> f_idx = {181, 182, 183, 184, 185, 279, 280, 281};

        // if (dim == 3)
        // {
        //     for (auto j: cof_idx)
        //     {
        //         std::cout<<"idx = "<<j<<"  match = "<<bi_graph.match_list[j]<<" "<<bi_graph.match_list[j] - int64_t(sorted_cofacet.size())<<"  weight = "<<sorted_cofacet[j].second<<"    adj = ";
        //         for (auto t : bi_graph.adj_list[j]) std::cout<<t<<"("<<t - sorted_cofacet.size()<<")"<<"  ";
        //         std::cout<<'\n';
        //     }

        //     for (auto j: f_idx)
        //     {
        //         std::cout<<"facet idx = "<<j<<"  weight = "<<sorted_simplex[j].second<<"  adj = ";
        //         for (auto t: bi_graph.adj_list[j + sorted_cofacet.size()])  std::cout<<t<<"  ";
        //         std::cout<<'\n';
        //     }
        // }
            

        std::cout<<bi_graph.dfsCycleRemoval()<<'\n';

        std::vector<size_t> crit_index = bi_graph.getCriticalIndex(dim_active_index, sorted_simplex.size());
        for(auto t: crit_index) std::cout<<"idx and weight = "<<t<<" "<<sorted_simplex[t].second<<"   ";
        std::cout<<'\n';

        dim_active_index = bi_graph.getActiveIndex();
        active_index_hash_table = getActiveFacetIndexHashTable(sorted_cofacet, dim_active_index);

        if (dim != maxdim)
        {
            sorted_simplex = getSortedCofacetList(binom_table, sorted_cofacet, dim, maxeps, 1);
            std::swap(sorted_simplex, sorted_cofacet);
        }

    }

    return;
}