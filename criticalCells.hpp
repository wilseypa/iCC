#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <chrono>
#include <Eigen/Sparse>

#include <set>

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

struct SparseDistMat
{
    Eigen::SparseMatrix<double> distMatrix;
    inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
};
struct NormallDistMat
{
    std::vector<std::vector<double>> distMatrix;
    inline double distance(size_t i, size_t j) const { return this->distMatrix[i][j]; };
};

struct VR // Creates a constructor for combinatorial simplexes n_pts = n, dim = r
{
public:
    std::vector<int> simplex;
    bool is_initialized; // Marks if the object is initialized
    bool active;         // Marks end of simplexes
    VR() : active(true), is_initialized(false) {};

    void initialize(size_t n_pts, size_t dim)
    {
        n_pts_ = n_pts;
        dim_ = dim;
        simplex.resize(dim); // Resize the vector to the appropriate size
        std::iota(simplex.begin(), simplex.end(), 0);
        is_initialized = true;
    };

    void next_simplex()
    {
        if (simplex[0] == n_pts_ - dim_)
        {
            active = false;
            return; // No more simplices can be generated.
        }

        size_t j = dim_ - 1;

        while (simplex[j] - j == n_pts_ - dim_)
            j--;

        simplex[j]++;

        for (size_t i = j + 1; i < dim_; i++)
            simplex[i] = simplex[j] + i - j;
    };

private:
    size_t n_pts_;
    size_t dim_;
};

template <typename ComplexType, typename DistMatType>
class CritCells : public DistMatType
{
public:
    CritCells(const std::string &fileName);               // FileName to read InputData from
    CritCells(std::vector<std::vector<double>> &distMat); // Input normal Distance matrix
    CritCells(Eigen::SparseMatrix<double> &distMat);
    void run_Compute(int maxDim, int batchsize = 0);
    
    /*bins by weight range*/
    //change later
    std::vector<std::vector<int>> getSortedEdges();
    std::vector<std::vector<int>> getEdgesByWeightRange(std::vector<std::vector<int>>& sorted_edges, double minweight, double maxweight);
    double getSimplexWeight(std::vector<int>& simplex);
    std::vector<std::vector<int>> getLEWeightCofacet(std::vector<std::vector<int>>& simplex_bin, int threadnum);
    int findRoot(std::vector<int>& parent_idx, int x);
    void setUnion(std::vector<int>& parent_idx, int x, int y);
    std::vector<int> getMSTEdgeIndices(std::vector<std::vector<int>>& sorted_edges);
    std::vector<std::vector<double>> run_MorseMatch(int maxdimension, double mineps, double maxeps);

private:
    std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                             // Direct creation of edgebins to a map
    void binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins); // Merged higher dim feature to bins
    std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType &simplex_const, size_t dim, size_t batch_size);                  // Worker is invokation counter
    std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final);
};