#include "deps/distMat.cpp"
#include "deps/readInput.hpp"
#include <chrono>

template <int dim, int n_pts>
struct dsimplicies
{
    std::array<int, dim> curr_state;
    dsimplicies()
    {
        for (int i = 0; i < dim; i++)
            curr_state[i] = 0;
    }
    void generate_next()
    {
        auto dimension = dim - 1;
        while (curr_state[dimension] == n_pts - 1)
            curr_state[dimension--] = 0;
        curr_state[dimension]++;
    }
};

std::map<double, std::vector<std::vector<int>>> binEdges(std::vector<std::vector<double>> &distMat)
{
    std::map<double, std::vector<std::vector<int>>> binned_edges;
    for (int i = 0; i < distMat.size(); i++)
        for (int j = i + 1; j < distMat.size(); j++)
            binned_edges[distMat[i][j]].push_back({i, j});
    return binned_edges;
}

template <int dim>
std::map<double, std::vector<std::vector<int>>> binByWeights(std::vector<std::pair<double, std::vector<std::vector<int>>>> &weighted_simplicies, std::map<double, std::vector<std::array<int, dim>>> &bins)
{
    for (auto &[weight, simplex] : weighted_simplicies)
        bins[weight].push_back(simplex);
    return bins;
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
    auto distMatrix = distMat(inputData);
    auto bins = binEdges(distMatrix);
    for (auto &[dist, vect] : bins)
        for (auto i : vect)
            std::cout << dist << " " << i[0] << " " << i[1] << std::endl;
    /*     for (int dim = 2; dim < inputData[0].size(); dim++)
        {
        } */
    auto end_time = std::chrono::high_resolution_clock::now();
    return 0;
}