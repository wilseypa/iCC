#include "PipelineCommon.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "SimplexUtility.hpp"

namespace
{
constexpr size_t MAX_SIMPLEX_DIMENSION =
    std::numeric_limits<uint64_t>::digits - 1;
}

void PipelineConfig::validate() const
{
    if (eps_breaks.empty())
        throw std::invalid_argument("Pipeline epsilon breaks must not be empty.");

    for (size_t i = 0; i < eps_breaks.size(); ++i)
    {
        const double epsilon = eps_breaks[i];
        if (!std::isfinite(epsilon) || epsilon <= 0.0)
            throw std::invalid_argument("Pipeline epsilon breaks must be finite and positive.");

        if (i > 0 && epsilon <= eps_breaks[i - 1])
            throw std::invalid_argument("Pipeline epsilon breaks must be strictly increasing.");
    }

    if (maxdim == 0 || maxdim > MAX_SIMPLEX_DIMENSION)
        throw std::invalid_argument("Pipeline maximum dimension must be in the range [1, 63].");

    if (threads <= 0)
        throw std::invalid_argument("Pipeline thread count must be positive.");

    if (!std::isfinite(pv_cap_scale) || pv_cap_scale <= 0.0)
        throw std::invalid_argument("Pipeline PV cap scale must be finite and positive.");

    if (!std::isfinite(pv_min_separation) || pv_min_separation < 0.0)
        throw std::invalid_argument("Pipeline PV minimum separation must be finite and nonnegative.");

    if (!std::isfinite(pv_cap_scale * eps_breaks.back()))
        throw std::invalid_argument("Pipeline PV cap at the final epsilon must be finite.");

    if (!std::isfinite(pv_min_separation * eps_breaks.back()))
        throw std::invalid_argument("Pipeline absolute PV minimum separation must be finite.");
}

double PipelineConfig::finalEpsilon() const
{
    if (eps_breaks.empty())
        throw std::logic_error("Final epsilon is unavailable before pipeline configuration validation.");

    return eps_breaks.back();
}

double PipelineConfig::absoluteMinSeparation() const
{
    return pv_min_separation * finalEpsilon();
}

PipelineRuntime::PipelineRuntime(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : distance_matrix_(std::move(distance_matrix)),
      config_(std::move(config))
{
    config_.validate();
    binomial_table_ = SimplexUtility::getBinomialTable(
        distance_matrix_.vertexCount(),
        config_.maxdim);
}

void PipelineRuntime::ensureBinomialCapacity(const size_t total_label_count)
{
    const size_t original_vertex_count = distance_matrix_.vertexCount();
    if (total_label_count < original_vertex_count)
        throw std::invalid_argument("Total label count cannot be smaller than the original vertex count.");

    if (binomial_table_.size() > total_label_count)
        return;

    // Grow a staged copy so an overflow or allocation failure cannot leave the
    // runtime advertising partially initialized rows.
    BinomialTable expanded_table = binomial_table_;
    SimplexUtility::updateBinomialTable(
        expanded_table,
        original_vertex_count,
        total_label_count - original_vertex_count,
        config_.maxdim);
    binomial_table_.swap(expanded_table);
}
