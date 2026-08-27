#pragma once

#include <cstddef>
#include <vector>

class DistanceMatrix
{
public:
    explicit DistanceMatrix(const std::vector<std::vector<double>>& point_cloud);

    [[nodiscard]]
    std::size_t vertexCount() const noexcept
    {
        return vertex_count_;
    }

    [[nodiscard]]
    double distance(const std::size_t i, const std::size_t j) const noexcept
    {
        if (i == j)
            return 0.0;

        const std::size_t hi = (i < j) ? j : i;
        const std::size_t lo = (i < j) ? i : j;
        return distances_[triangularIndex(hi, lo)];
    }

private:
    static std::size_t triangularIndex(
        const std::size_t hi,
        const std::size_t lo) noexcept
    {
        return (hi * (hi - 1)) / 2 + lo;
    }

    std::size_t vertex_count_ = 0;
    std::vector<double> distances_;
};
