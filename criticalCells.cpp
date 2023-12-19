#include "deps/distMat.cpp"
#include "deps/readInput.hpp"
#include <chrono>

std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes(std::vector<std::vector<double>> &distMat) // Direct creation of edgebins to a map
{
    std::map<double, std::vector<std::vector<int>>> binned_edges;
    for (int i = 0; i < distMat.size(); i++)
        for (int j = i + 1; j < distMat.size(); j++)
            binned_edges[distMat[i][j]].push_back({i, j});
    return binned_edges;
}

void binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins) // Merged higher dim feature to bins
{
    for (auto &[weight, simplexes] : weighted_simplicies)
        for (auto simplex : simplexes)
            bins[weight].push_back(simplex);
    return;
}

unsigned long long combinations(int n, int r)
{
    if (r > n)
        return 0;
    r = std::min(r, n - r);
    unsigned long long result = 1;
    for (int i = 1; i <= r; ++i)
    {
        result *= (n - r + i);
        result /= i;
    }
    return result;
}

struct dsimplexes // Creates a constructor for combinatorial simplexes n_pts = n, dim = r
{
    // Usage
    //  dsimplexes ds(N, R);
    //  do
    //  {
    //      ds.print_simplex();
    //  } while (ds.next_simplex());

    size_t n_pts;
    size_t dim;
    std::vector<int> simplex;

    dsimplexes(size_t n_pts, size_t dim) : n_pts(n_pts), dim(dim)
    {
        for (size_t i = dim; i > 0; i--)
            simplex.push_back(i - 1);
    };

    dsimplexes(size_t n_pts, size_t dim, size_t start_at) : n_pts(n_pts), dim(dim)
    {
        auto temp = start_at;
    };

    bool next_simplex()
    {
        if (simplex.back() == n_pts - dim)
            return false;
        for (size_t i = 0; i < dim; i++)
        {
            bool flag = true;
            if (++simplex[i] == n_pts - i)
            {
                flag = false;
                auto reset_to = n_pts;
                for (size_t j = i + 1; j < dim && reset_to >= n_pts; j++)
                    reset_to = simplex[j] + j + 1;
                simplex[i] = reset_to - i;
            }
            else if (flag)
                return true;
        }
        return false; // Not really needed just avoiding compiler warnings
    };
};

std::map<double, std::vector<std::vector<int>>> dsimplices_batches(std::vector<std::vector<double>> &distMat, u_int dim, size_t batch_size, u_int worker) // Worker is invokation counter
{
    std::map<double, std::vector<std::vector<int>>> weighted_simplexes;
    // worker param not currently in use.. Will get added based on pre initialization of dsimplexes struct
    // dsimplexes ds(distMat.size(), dim, batch_size*worker);
    dsimplexes ds(distMat.size(), dim);

    do
    {
        double max_dist = 0;
        for (int i = 0; i < dim - 1; i++)
            for (int j = i + 1; j < dim; j++)
                max_dist = std::max(max_dist, distMat[ds.simplex[i]][ds.simplex[j]]);
        weighted_simplexes[max_dist].push_back(ds.simplex);
    } while (ds.next_simplex());
    return weighted_simplexes;
}

std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>> simplexes, size_t dim, bool final = false)
{
    return simplexes;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_filename>" << std::endl;
        return 1;
    }
    auto start_time = std::chrono::high_resolution_clock::now();
    auto inputData = readInput::readCSV(argv[1]);
    size_t maxDim = inputData[0].size();
    auto distMatrix = distMat(inputData);
    auto bins = binEdgeSimplexes(distMatrix);

    for (size_t dim = 2; dim < maxDim; dim++)
    {
        auto batch_start_time = std::chrono::high_resolution_clock::now();
        auto weighted_simplicies = dsimplices_batches(distMatrix, dim + 1, 1000, 0); // Worker is invokation counter
        auto batch_end_time = std::chrono::high_resolution_clock::now();
        std::cout << "Batch time" << (batch_end_time - batch_start_time).count() << std::endl;

        auto match_start_time = std::chrono::high_resolution_clock::now();

        // Bin the batches
        binByWeights(weighted_simplicies, bins);

        // Dim Matching functionality
        for (auto &[weight, simplexes] : bins)
            simplexes = dimMatching(simplexes, dim, dim == maxDim - 1);
        auto match_end_time = std::chrono::high_resolution_clock::now();
        std::cout << "Match time" << (match_end_time - match_start_time).count() << std::endl;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = (end_time - start_time).count();
    std::cout << "Elapsed time: " << duration << " ns" << std::endl;
    return 0;
}