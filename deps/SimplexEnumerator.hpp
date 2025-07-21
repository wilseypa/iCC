#include "SimplexUtility.hpp"

typedef NormallDistMat;

class SimplexEnumerator
{
public:
    SimplexEnumerator(const NormallDistMat& dist_mat, const std::vector<std::vector<int64_t>>& binomial_table)
        : dist_matrix(dist_mat), binomial_table(binomial_table) {}

private:
    const NormallDistMat& dist_matrix;
    const std::vector<std::vector<int64_t>>& binomial_table;
};