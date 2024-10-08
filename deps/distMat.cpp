#include <vector>
#include <algorithm>
#include <cmath>
#include <execution>

#ifdef _GLIBCXX_DEBUG 
#include <stdexcept>
#endif

inline double vectors_distance(const std::vector<double>& a, const std::vector<double>& b) {
#ifdef _GLIBCXX_DEBUG 
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must be of the same length");
  }
  if (a.empty()) {
    throw std::invalid_argument("Vectors must not be empty");
  }
#endif
  return sqrt(std::transform_reduce(std::execution::par, a.cbegin(), a.cend(), b.cbegin(), 0.0, std::plus<>(),
    [](double e1, double e2) { return (e1 - e2) * (e1 - e2); }));
}

std::vector<std::vector<double>> distMat(std::vector<std::vector<double>>& inputData)
{
  std::vector<std::vector<double>> distMatrix(inputData.size(), std::vector<double>(inputData.size(), 0));
  for (size_t i = 0; i < inputData.size(); i++)
    for (size_t j = i + 1; j < inputData.size(); j++)
      distMatrix[i][j] = vectors_distance(inputData[i], inputData[j]);
  return distMatrix;
}