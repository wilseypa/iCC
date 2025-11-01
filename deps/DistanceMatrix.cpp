#include <cmath>
#include <algorithm>
#include <execution>
#include <stdexcept>
#include <omp.h>
#include "DistanceMatrix.hpp"

inline double vectors_distance(const std::vector<double> &a, const std::vector<double> &b)
{
#ifdef _GLIBCXX_DEBUG
    if (a.size() != b.size())
    {
        throw std::invalid_argument("Vectors must be of the same length");
    }
    if (a.empty())
    {
        throw std::invalid_argument("Vectors must not be empty");
    }
#endif
    return sqrt(std::transform_reduce(std::execution::par, a.cbegin(), a.cend(), b.cbegin(), 0.0, std::plus<>(),
                                      [](double e1, double e2)
                                      { return (e1 - e2) * (e1 - e2); }));
}

NormalDistMat::NormalDistMat(const std::vector<std::vector<double>> &point_cloud)
{
    this->dist_mat_.resize(point_cloud.size(), std::vector<double>(point_cloud.size(), 0.0));
#pragma omp parallel for
    for (size_t i = 0; i < point_cloud.size() - 1; i++)
    {
        for (size_t j = i + 1; j < point_cloud.size(); j++)
        {
            this->dist_mat_[i][j] = vectors_distance(point_cloud[i], point_cloud[j]);
            // distMatrix[j][i] = dist_mat[i][j] ; // Symmetric matrix
        }
    }
}