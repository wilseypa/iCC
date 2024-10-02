#include "deps/distMat.cpp"
#include "deps/readInput.hpp"
#include "deps/hopcroft_karp.hpp"
#include "criticalCells.hpp"
#include <omp.h>

// #define PARALLEL

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
    auto bins = binEdgeSimplexes();
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
                    dimMatching(it->second, dim, dim == maxDim);
            }

#else
            for (auto &it : bins)
                dimMatching(it.second, dim, dim == maxDim);
#endif
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
        bins[weight].emplace_back(std::move(simplexes));
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
void CritCells<ComplexType, DistMatType>::dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final)
{
    std::sort(simplexes.begin(), simplexes.end(), [](const auto &first, const auto &second)
              { return first.size() < second.size(); });
    auto simps_iter = std::find_if(simplexes.begin(), simplexes.end(), [dim](const auto &vect)
                                   { return vect.size() == dim; });

    auto cofaces_iter = std::find_if(simplexes.begin(), simplexes.end(), [dim](const auto &vect)
                                     { return vect.size() == dim + 1; });

    // Check if there are no simplices or cofaces
    if (simps_iter == cofaces_iter || cofaces_iter == simplexes.end())
    {
        if (final)
            simplexes.erase(cofaces_iter, simplexes.end());
        return; // Exit early if there's nothing to process
    }

    size_t num_simps = std::distance(simps_iter, cofaces_iter);
    size_t num_cofaces = std::distance(cofaces_iter, simplexes.end());

    HKGraph csr_matrix(num_simps, num_cofaces);

    for (size_t i = 0; i < num_simps; ++i)
    {
        for (size_t j = 0; j < num_cofaces; ++j)
        {
            if (std::includes((cofaces_iter + j)->begin(), (cofaces_iter + j)->end(), (simps_iter + i)->begin(), (simps_iter + i)->end()))
            {
                csr_matrix.addEdge(i + 1, j + 1); // Initial index 1
            }
        }
    }

    auto res = csr_matrix.custom_hopcroftKarpAlgorithm();

    if (final)
    {
        simplexes.erase(cofaces_iter, simplexes.end());
    }
    else
    {
        for (auto index : res.second)
        {
            simplexes.erase(cofaces_iter + index - 1);
        }
    }

    for (auto index : res.first)
    {
        simplexes.erase(simps_iter + index - 1);
    }
}
