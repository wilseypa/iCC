#include "deps/distMat.cpp"
#include "deps/readInput.hpp"
#include "deps/hopcroft_karp.hpp"
#include "criticalCells.hpp"
#include <omp.h>

#define PARALLEL

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
    auto bins = binEdgeSimplexes();
    for (size_t dim = 2; dim <= maxDim; dim++)
    {
        auto batch_start_time = std::chrono::high_resolution_clock::now();
        VR complex(this->distMatrix.size(), dim + 1);
        auto weighted_simplicies = dsimplices_batches(complex, dim, batch_size); // Worker is invokation counter
        auto batch_end_time = std::chrono::high_resolution_clock::now();
        std::cout << "Batch time: " << std::chrono::duration<double>(batch_end_time - batch_start_time).count() << " seconds." << std::endl;
        auto match_start_time = std::chrono::high_resolution_clock::now();
        // Bin the batches
        binByWeights(weighted_simplicies, bins);

        // Dim Matching functionality
        int num_threads = 20;
#ifdef PARALLEL
#pragma omp parallel for
        for (int i = 0; i < num_threads; i++)
        {
            size_t block_size = floor((double)bins.size() / num_threads);
            auto it = std::next(bins.begin(), block_size * i);
            bool isLastThread = (i + 1 == num_threads);
            auto end = isLastThread ? bins.end() : std::next(bins.begin(), block_size * (i + 1));
            for (; it != end; it++)
                it->second = dimMatching(it->second, dim, dim == maxDim);
        }

#else
        for (auto &it : bins)
            it.second = dimMatching(it.second, dim, dim == maxDim);
#endif
        auto match_end_time = std::chrono::high_resolution_clock::now();
        std::cout << "Match time: " << std::chrono::duration<double>(match_end_time - match_start_time).count() << " seconds." << std::endl;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration<double>(end_time - start_time).count() << " seconds." << std::endl;
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
std::map<double, std::vector<std::vector<int>>> CritCells<ComplexType, DistMatType>::dsimplices_batches(ComplexType complex, size_t dim, size_t batch_size) // Worker is invokation counter
{
    std::map<double, std::vector<std::vector<int>>> weighted_simplexes;
    do
    {
        double max_dist = 0;
        for (int i = 0; i < dim; i++)
            for (int j = i + 1; j <= dim; j++)
                max_dist = std::max(max_dist, this->distance(complex.simplex[j], complex.simplex[i]));
        weighted_simplexes[max_dist].push_back({complex.simplex.rbegin(), complex.simplex.rend()});
    } while (complex.next_simplex());
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
    std::for_each(res.second.rbegin(), res.second.rend(), [&cofaces](auto index)
                  { cofaces.erase(std::next(cofaces.begin(), index - 1)); });

    std::move(simps.begin(), simps.end(), std::back_inserter(critCells));
    if (!final)
        std::move(cofaces.begin(), cofaces.end(), std::back_inserter(critCells));
    return critCells;
}