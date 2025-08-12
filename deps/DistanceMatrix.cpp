#include <cmath>
#include <algorithm>
#include <execution>

#include <omp.h>

#include "DistanceMatrix.hpp"

namespace    //static namespace
{
    inline double euclideanDistance(const std::vector<double>& vec_0, const std::vector<double>& vec_1)
    {
        if (vec_0.size() != vec_1.size()) 
        {
            throw std::invalid_argument("Vectors must be of the same length for distance calculation.");
        }

        double sum = 0.0;

        for (size_t i = 0; i < vec_0.size(); i++)
        {
            double diff = vec_0[i] - vec_1[i];
            sum += diff * diff;
        }

        return std::sqrt(sum);
    }
}


NormalDistMat::NormalDistMat(const std::vector<std::vector<double>>& point_cloud, int threadnum = 4)
{
    const size_t n = point_cloud.size();
    if (n == 0) return;

    dist_matrix.resize(n, std::vector<double>(n, 0.0));

    omp_set_num_threads(threadnum);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            double dist = euclideanDistance(point_cloud[i], point_cloud[j]);
            dist_matrix[i][j] = dist;
            // distMatrix[j][i] = dist; // Symmetric matrix
        }
    }
}