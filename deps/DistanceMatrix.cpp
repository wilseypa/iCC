#include <cmath>
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
    double squared_distance = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        const double diff = a[i] - b[i];
        squared_distance += diff * diff;
    }
    return std::sqrt(squared_distance);
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
