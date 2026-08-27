#include "DistanceMatrix.hpp"

#include <cmath>
#include <stdexcept>

namespace
{
double vectorDistance(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    double squared_distance = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        const double difference = lhs[i] - rhs[i];
        squared_distance += difference * difference;
    }
    return std::sqrt(squared_distance);
}
}

DistanceMatrix::DistanceMatrix(
    const std::vector<std::vector<double>>& point_cloud)
{
    if (point_cloud.empty())
        throw std::invalid_argument("Input point cloud is empty.");

    const std::size_t coordinate_count = point_cloud.front().size();
    if (coordinate_count == 0)
        throw std::invalid_argument("Input points must contain at least one coordinate.");

    for (const auto& point : point_cloud)
    {
        if (point.size() != coordinate_count)
            throw std::invalid_argument("Input points must have equal coordinate counts.");

        for (const double coordinate : point)
        {
            if (!std::isfinite(coordinate))
                throw std::invalid_argument("Input point coordinates must be finite.");
        }
    }

    vertex_count_ = point_cloud.size();
    distances_.assign((vertex_count_ * (vertex_count_ - 1)) / 2, 0.0);
    if (vertex_count_ == 1)
        return;

#pragma omp parallel for
    for (std::size_t i = 0; i < vertex_count_ - 1; ++i)
    {
        for (std::size_t j = i + 1; j < vertex_count_; ++j)
        {
            distances_[triangularIndex(j, i)] =
                vectorDistance(point_cloud[i], point_cloud[j]);
        }
    }

    for (const double distance : distances_)
    {
        if (!std::isfinite(distance))
            throw std::invalid_argument("Input coordinates produce a nonfinite Euclidean distance.");
    }
}
