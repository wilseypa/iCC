#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <chrono>

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

struct NormallDistMat
{
    std::vector<std::vector<double>> distMatrix;
    inline double distance(size_t i, size_t j) { return this->distMatrix[i][j]; };
};

struct VR // Creates a constructor for combinatorial simplexes n_pts = n, dim = r
{
public:
    std::vector<int> simplex;
    bool active;

    VR(size_t n_pts, size_t dim) : n_pts(n_pts), dim(dim), active(true)
    {
        for (size_t i = dim; i > 0; i--)
            simplex.push_back(i - 1);
    };

    bool next_simplex()
    {
        if (simplex.back() != n_pts - dim)
        {
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
        }
        active = false;
        return false; // Not really needed just avoiding compiler warnings
    };

private:
    size_t n_pts;
    size_t dim;
};

template <typename ComplexType, typename DistMatType>
class CritCells : public DistMatType
{
public:
    CritCells(const std::string &fileName);               // FileName to read InputData from
    CritCells(std::vector<std::vector<double>> &distMat); // Input normal Distance matrix
    void run_Compute(int maxDim, int batchsize = 0);

private:
    std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                             // Direct creation of edgebins to a map
    void binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins); // Merged higher dim feature to bins
    std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType &simplex_const, size_t dim, size_t batch_size);                   // Worker is invokation counter
    std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final);
};