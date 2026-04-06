#pragma once

#include <vector>

#ifdef BUILD_ALPHA_COMPLEX
#include <Eigen/Sparse>
#endif

struct NormalDistMat
{
public:
    NormalDistMat(const std::vector<std::vector<double>> &point_cloud);

    inline std::size_t getVertexNumber() const { return dist_mat_.size(); }
    inline double getDistance(const size_t i, const size_t j) const
    {
        return (i < j) ? this->dist_mat_[i][j] : this->dist_mat_[j][i];
    };

private:
    std::vector<std::vector<double>> dist_mat_;
};

struct SparseDistMat
{
#ifdef BUILD_ALPHA_COMPLEX
    Eigen::SparseMatrix<double> distMatrix;
    inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
#endif
};
