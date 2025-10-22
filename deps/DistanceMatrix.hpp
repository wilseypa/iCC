#pragma once

#include <vector>

#include <Eigen/Sparse>

struct NormalDistMat
{
public:
    NormalDistMat(const std::vector<std::vector<double>> &point_cloud);

    inline std::size_t getVertexNumber() const { return dist_mat.size(); }
    inline double getDistance(const size_t i, const size_t j) const
    {
        return (i < j) ? this->dist_mat[i][j] : this->dist_mat[j][i];
    };

private:
    std::vector<std::vector<double>> dist_mat;
};

struct SparseDistMat
{
    Eigen::SparseMatrix<double> distMatrix;
    inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
};