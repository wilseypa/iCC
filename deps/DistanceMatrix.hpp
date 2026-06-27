#pragma once

#include <cstddef>
#include <vector>

#ifdef BUILD_ALPHA_COMPLEX
#include <Eigen/Sparse>
#endif

struct NormalDistMat
{
public:
    NormalDistMat(const std::vector<std::vector<double>>& point_cloud);

    inline std::size_t getVertexNumber() const { return vertex_count_; }
    inline double getDistance(const size_t i, const size_t j) const
    {
        if (i == j)
            return 0.0;

        const size_t hi = (i < j) ? j : i;
        const size_t lo = (i < j) ? i : j;
        return dist_mat_[triangularIndex(hi, lo)];
    };

private:
    static inline size_t triangularIndex(const size_t hi, const size_t lo)
    {
        return (hi * (hi - 1)) / 2 + lo;
    }

    size_t vertex_count_ = 0;
    std::vector<double> dist_mat_;
};

struct SparseDistMat
{
#ifdef BUILD_ALPHA_COMPLEX
    Eigen::SparseMatrix<double> distMatrix;
    inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
#endif
};
