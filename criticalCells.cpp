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

std::map<double, std::vector<std::vector<int>>> binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins) // Merged higher dim feature to bins
{
    for (auto &[weight, simplexes] : weighted_simplicies)
        for (auto simplex : simplexes)
            bins[weight].push_back(simplex);
    return bins;
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
    // worker param not currently in use.. Will get added based on pre initialization of dsimplexes struct
    std::map<double, std::vector<std::vector<int>>> weighted_simplexes;
    dsimplexes ds(distMat.size(), dim);
    do
    {
        double max_dist = 0;
        for (int i = 0; i < dim - 1; i++)
            for (int j = i + 1; j < dim; j++)
            {
                auto dist = distMat[ds.simplex[i]][ds.simplex[j]];
                if (max_dist < dist)
                    max_dist = dist;
            }
        weighted_simplexes[max_dist].push_back(ds.simplex);
    } while (ds.next_simplex());
    return weighted_simplexes;
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
    auto maxDim = inputData[0].size() + 2;
    auto distMatrix = distMat(inputData);
    auto bins = binEdgeSimplexes(distMatrix);

    for (size_t dim = 2; dim < maxDim; dim++)
    {
        auto batch_time = std::chrono::high_resolution_clock::now();
        auto weighted_simplicies = dsimplices_batches(distMatrix, 3, 1000, 0); // Worker is invokation counter

        // for (auto &[dist, vect] : weighted_simplicies)
        //    for (auto i : vect)
        //        std::cout << dist << " " << i[0] << " " << i[1]  << " " << i[2]<< std::endl;
        /*     for (int dim = 2; dim < inputData[0].size(); dim++)
            {
            } */
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Elapsed time: " << duration << " milliseconds" << std::endl;
    return 0;
}