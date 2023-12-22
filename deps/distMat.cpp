#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

double vectors_distance(const std::vector<double> &a, const std::vector<double> &b)
{
  std::vector<double> temp;
  if (b.size() == 0)
    return 0;
  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](double e1, double e2)
                 { return pow((e1 - e2), 2); });
  return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}

std::vector<std::vector<double>> distMat(std::vector<std::vector<double>> &inputData)
{
  std::vector<std::vector<double>> distMatrix;
  if (distMatrix.size() > 0)
    distMatrix.clear();
  distMatrix.resize(inputData.size(), std::vector<double>(inputData.size(), 0));
  for (size_t i = 0; i < inputData.size(); i++)
  {
    for (size_t j = i + 1; j < inputData.size(); j++)
    {
      distMatrix[i][j] = vectors_distance(inputData[i], inputData[j]);
      distMatrix[j][i] = distMatrix[i][j];
    }
  }
  return distMatrix;
}