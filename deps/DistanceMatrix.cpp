#include <cmath>
#include <stdexcept>
#include <omp.h>
#include "DistanceMatrix.hpp"

inline double vectors_distance(const std::vector<double>& a, const std::vector<double>& b)
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

NormalDistMat::NormalDistMat(const std::vector<std::vector<double>>& point_cloud)
{

    if (point_cloud.empty()) 
    {
        throw std::invalid_argument("Input point cloud is empty.");
    }

    vertex_count_ = point_cloud.size();
    dist_mat_.assign((vertex_count_ * (vertex_count_ - 1)) / 2, 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < vertex_count_ - 1; i++)
    {
        for (size_t j = i + 1; j < vertex_count_; j++)
        {
            dist_mat_[triangularIndex(j, i)] = vectors_distance(point_cloud[i], point_cloud[j]);
        }
    }
}
