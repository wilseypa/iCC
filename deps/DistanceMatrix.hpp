#pragma once

#include <vector>

// #include <Eigen/Sparse>


struct NormalDistMat
{
    std::vector<std::vector<double>> dist_mat;

    inline double getDistance(const size_t i, const size_t j) const 
    {
        if (i < j) return this->dist_mat[i][j];
        
        return this->dist_mat[j][i];
    };

    NormalDistMat() = default;

    explicit NormalDistMat(const std::vector<std::vector<double>>& point_cloud, int threadnum = 4);
};

struct SparseDistMat
{};

// struct SparseDistMat
// {
//     Eigen::SparseMatrix<double> distMatrix;
//     inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
// };