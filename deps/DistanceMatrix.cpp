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


NormalDistMat::NormalDistMat(const std::vector<std::vector<double>>& point_cloud, int threadnum)
{
    const size_t n = point_cloud.size();
    if (n == 0) return;

    dist_mat.resize(n, std::vector<double>(n, 0.0));

    omp_set_num_threads(threadnum);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            dist_mat[i][j] = vectors_distance(point_cloud[i], point_cloud[j]);
            // distMatrix[j][i] = dist_mat[i][j] ; // Symmetric matrix
        }
    }
}