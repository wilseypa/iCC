#include "distMat.cpp"
#include "readInput.hpp"
#include "hopcroft_karp.hpp"
#include "criticalCells.hpp"
#include "biGraphMorseMatch.hpp"
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
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::getSortedEdges()
{
    auto sort_lambda = [this](std::vector<int> &edge_i, std::vector<int> &edge_j)
    { return (this->distMatrix[edge_i[0]][edge_j[1]] < this->distMatrix[edge_j[0]][edge_j[1]]); };
    std::vector<std::vector<int>> edge_vec;
    for (int i = 0; i < this->distMatrix.size() - 1; i++)
    {
        for (int j = i + 1; j < this->distMatrix.size(); j++)
        {
            edge_vec.push_back(std::vector<int>{i, j});
        }
    }

    //before sort
    for (int i = 0; i < 10; i++)
    {   
        std::cout<<edge_vec[i][0]<<" "<<edge_vec[i][1]<<"  ";
        std::cout<<getSimplexWeight(edge_vec[i])<<'\n';
    }
    std::cout<<'\n';

    std::ranges::sort(edge_vec.begin(), edge_vec.end(), sort_lambda);

    //after sort
    for (int i = 0; i < 10; i++)
    {   
        std::cout<<edge_vec[i][0]<<" "<<edge_vec[i][1]<<"  ";
        std::cout<<getSimplexWeight(edge_vec[i])<<'\n';
    }
    std::cout<<'\n';

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
double CritCells<ComplexType, DistMatType>::getSimplexWeight(std::vector<int> &simplex)
{
    double maxweight = 0;
    for (int i = 0; i < simplex.size() - 1; i++)
    {   
        int m = simplex[i];
        for (int j = i + 1; j < simplex.size(); j++)
        {
            int n = simplex[j];
            if (m > n) std::swap(m, n);
            if (this->distMatrix[m][n] > maxweight)
                maxweight = this->distMatrix[m][n];
        }
    }
    return maxweight;
}

template <typename ComplexType, typename DistMatType>
std::vector<std::vector<int>> CritCells<ComplexType, DistMatType>::getLEWeightCofacet(std::vector<std::vector<int>> &simplexes, int threadnum)
{
    std::vector<std::vector<int>> cofacet_bin;
    size_t npts = this->distMatrix.size();

    omp_set_num_threads(1);
    std::vector<std::vector<int>> thread_workspace(simplexes.size(), std::vector<int>());

#pragma omp parallel for
    for (size_t i = 0; i < simplexes.size(); i++)
    {
        double simpweight = getSimplexWeight(simplexes[i]);

        int last_idx = *std::min_element(simplexes[i].begin(), simplexes[i].end());
        for (int j = 0; j < last_idx; j++)
        {
            auto maxidx = *std::max_element(simplexes[i].begin(), simplexes[i].end(), [this, j](auto first, auto second)
                                            { return this->distMatrix[j][first] < this->distMatrix[j][second]; });
            
            double weight = this->distMatrix[j][maxidx];

            if (weight < simpweight) thread_workspace[i].push_back(j);
        }
    }
    for (size_t i = 0; i < simplexes.size(); i++)
    {
        for (auto &pt : thread_workspace[i])
        {
            cofacet_bin.emplace_back(simplexes[i].begin(), simplexes[i].end());
            cofacet_bin.back().push_back(pt);
        }
    }
    return cofacet_bin;
}

template <typename ComplexType, typename DistMatType>
int CritCells<ComplexType, DistMatType>::findRoot(std::vector<int>& parent_idx, int x)
{
    while (parent_idx[x] != x)
    {
        parent_idx[x] = parent_idx[parent_idx[x]];
        x = parent_idx[x];
    }
    return x;
}

template <typename ComplexType, typename DistMatType>
void CritCells<ComplexType, DistMatType>::setUnion(std::vector<int>& parent_idx, int x, int y)
{
    int p = findRoot(parent_idx, x);
    int q = findRoot(parent_idx, y);
    parent_idx[p] = parent_idx[q];
}

template <typename ComplexType, typename DistMatType>
std::vector<int> CritCells<ComplexType, DistMatType>::getMSTEdgeIndices(std::vector<std::vector<int>>& sorted_edges)
{
    int n = this->distMatrix.size();

    std::vector<int> parent_idx(n);
    std::iota(parent_idx.begin(), parent_idx.end(), 0);

    std::vector<int> active_edge_index;

    int x, y;
    int m = sorted_edges.size();
    for (int i = 0; i < m; i++)
    {
        x = sorted_edges[i][0];
        y = sorted_edges[i][1];
        if (findRoot(parent_idx, x) != findRoot(parent_idx, y))
        {
            setUnion(parent_idx, x, y);
            active_edge_index.push_back(i);
        }
    }

    return active_edge_index;
}


template <typename ComplexType, typename DistMatType>
std::vector<std::vector<double>> CritCells<ComplexType, DistMatType>::run_MorseMatch(int maxdimension, double mineps, double maxeps)
{
    int threadnum = 1;

    std::vector<std::vector<int>> sorted_edges = getSortedEdges();

    std::vector<std::vector<int>> simplex_bin = getEdgesByWeightRange(sorted_edges, mineps, maxeps);

    // for (int i = 0; i < 10; i++)
    // {
    //     std::cout<<getSimplexWeight(sorted_edges[i])<<"  ";
    // }
    // std::cout<<'\n';

    std::vector<int> mst_edge_index = getMSTEdgeIndices(simplex_bin);

    std::vector<int> dim_active_index;
    for (int i = 0; i < simplex_bin.size(); i++)
    {
        if (std::find(mst_edge_index.begin(), mst_edge_index.end(), i) == mst_edge_index.end()) dim_active_index.push_back(i);
    }
    
    std::vector<std::vector<double>> critical_weight;
    critical_weight.push_back(std::vector<double>{0});

    int initleftdeg = 3;    //cofacet of edge, tri
    
    std::vector<std::vector<int>> cofacet_bin = getLEWeightCofacet(simplex_bin, threadnum);

    Bi_Graph_Match bi_graph(cofacet_bin.size(), simplex_bin.size(), initleftdeg, threadnum);

    for (int dim = 2; dim <= maxdimension; dim++)
    {
        // std::cout<<"dim active idx size = "<<dim_active_index.size()<<'\n';
        bi_graph.buildInterface(cofacet_bin, simplex_bin, dim_active_index);

        bi_graph.parallelDFSMatch();
        bi_graph.serialCycleRemoval();

        std::vector<int> dim_critical_index = bi_graph.getCriticalIndex(dim_active_index);
        std::cout<<"active idx"<<'\n';
        for(auto& i: dim_active_index) std::cout<<i<<" ";
        std::cout<<'\n';
        std::cout<<'\n';
        std::cout<<"critical idx"<<'\n';
        for(auto& i: dim_critical_index) std::cout<<i<<" ";
        std::cout<<'\n';
        std::cout<<'\n';

        dim_active_index = bi_graph.getActiveIndex();
        
        // std for each
        std::vector<double> dim_critical_weight;
        for (auto i : dim_critical_index)
        {
            double weight = getSimplexWeight(simplex_bin[i]);
            dim_critical_weight.push_back(weight);
        }
        critical_weight.push_back(std::move(dim_critical_weight));

        if (dim < maxdimension)
        {   
            simplex_bin = getLEWeightCofacet(cofacet_bin, threadnum);
            std::swap(simplex_bin, cofacet_bin);
            bi_graph.updateDimension(cofacet_bin.size(), simplex_bin.size());  
        } 
        else
        {   
            for (auto i: dim_active_index)
            {   
                double weight = getSimplexWeight(cofacet_bin[i]);
                dim_critical_weight.push_back(weight);
            }
            critical_weight.push_back(std::move(dim_critical_weight));
        }
        
    }

    return critical_weight;
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
        double max_dist = 0;
        for (int i = 0; i < dim; i++)
            for (int j = i + 1; j <= dim; j++)
                max_dist = std::max(max_dist, this->distance(complex.simplex[i], complex.simplex[j]));
        weighted_simplexes[max_dist].push_back(complex.simplex);
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