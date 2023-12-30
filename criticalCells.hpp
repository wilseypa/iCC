#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>

std::ostream &operator<<(std::ostream &os, const std::map<double, std::vector<std::vector<int>>> &bins)
{
    for (const auto &[weight, simplexes] : bins)
    {
        if (!simplexes.empty())
        {
            os << weight << " ";
            for (const auto &simplex : simplexes)
            {
                os << "(";
                for (const auto &i : simplex)
                    os << i << " ";
                os << ")";
            }
            os << std::endl;
        }
    }
    return os;
}

struct SparesDistMat
{
    std::vector<std::vector<double>> distMatrix;
    double distance(size_t i, size_t j) { return this->distMatrix[i][j]; };
};

struct NominalDistMat
{
    std::vector<std::vector<double>> distMatrix;
    inline double distance(size_t i, size_t j) { return this->distMatrix[i][j]; };
};

template <typename DistMatType>
class CritCells : public DistMatType
{
public:
    CritCells(const std::string &fileName);               // FileName to read InputData from
    CritCells(std::vector<std::vector<double>> &distMat); // Input normal Distance matrix
    void run_Compute(int maxDim, int batchsize = 50000);

private:
    std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                             // Direct creation of edgebins to a map
    void binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins); // Merged higher dim feature to bins
    std::map<double, std::vector<std::vector<int>>> dsimplices_batches(size_t dim, size_t batch_size);                                              // Worker is invokation counter
    std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final);
};

struct dsimplexes // Creates a constructor for combinatorial simplexes n_pts = n, dim = r
{
public:
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

private:
    size_t n_pts;
    size_t dim;
};

unsigned long long combinations(int n, int r)
{
    if (r > n)
        return 0;
    r = std::min(r, n - r);
    unsigned long long result = 1;
    for (int i = 1; i <= r; ++i)
        result *= (n - r + i) / i;
    return result;
}