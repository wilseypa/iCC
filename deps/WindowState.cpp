#include "WindowState.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <omp.h>

#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"

namespace
{
void validateBounds(
    const PipelineConfig& config,
    const WindowBounds& bounds)
{
    if (bounds.index >= config.eps_breaks.size())
        throw std::invalid_argument("Window index is outside the configured epsilon schedule.");

    if (!std::isfinite(bounds.eps_lo) || bounds.eps_lo < 0.0 ||
        !std::isfinite(bounds.eps_hi) || bounds.eps_hi <= 0.0 ||
        bounds.eps_lo >= bounds.eps_hi)
    {
        throw std::invalid_argument("Window bounds must be finite and satisfy 0 <= eps_lo < eps_hi.");
    }

    const double expected_lo =
        (bounds.index == 0) ? 0.0 : config.eps_breaks[bounds.index - 1];
    const double expected_hi = config.eps_breaks[bounds.index];
    const bool expected_final = bounds.index + 1 == config.eps_breaks.size();

    if (bounds.eps_lo != expected_lo || bounds.eps_hi != expected_hi)
        throw std::invalid_argument("Window bounds do not match the configured epsilon schedule.");

    if (bounds.is_final != expected_final)
        throw std::invalid_argument("Window final-state flag does not match the configured epsilon schedule.");
}

void validateActiveLabels(
    const std::vector<size_t>& active_labels,
    const size_t total_label_count)
{
    if (!std::is_sorted(active_labels.begin(), active_labels.end()))
        throw std::logic_error("Active labels must be sorted.");

    for (size_t i = 0; i < active_labels.size(); ++i)
    {
        if (active_labels[i] >= total_label_count)
            throw std::logic_error("Active label is outside the complete label range.");

        if (i > 0 && active_labels[i] == active_labels[i - 1])
            throw std::logic_error("Active labels must not contain duplicates.");
    }
}

std::vector<std::vector<size_t>> buildRepresentativeLists(
    const std::vector<SelectedPV>& pv_list,
    const size_t original_vertex_count)
{
    std::vector<std::vector<size_t>> representatives;
    representatives.reserve(pv_list.size());

    for (const auto& pv : pv_list)
    {
        if (pv.flat_index_set.empty())
            throw std::logic_error("A pseudo-vertex cannot have an empty representative set.");

        if (!std::isfinite(pv.diameter) || pv.diameter < 0.0)
            throw std::logic_error("A pseudo-vertex diameter must be finite and nonnegative.");

        auto& sorted_representatives = representatives.emplace_back(
            pv.flat_index_set.begin(),
            pv.flat_index_set.end());
        std::sort(sorted_representatives.begin(), sorted_representatives.end());

        if (sorted_representatives.back() >= original_vertex_count)
            throw std::logic_error("Pseudo-vertex representatives must be original vertex indices.");
    }

    return representatives;
}

LabelPairKey makeLabelPairKey(size_t lhs, size_t rhs)
{
    if (lhs > rhs)
        std::swap(lhs, rhs);

    constexpr size_t MAX_PACKED_LABEL = std::numeric_limits<uint32_t>::max();
    if (lhs > MAX_PACKED_LABEL || rhs > MAX_PACKED_LABEL)
        throw std::overflow_error("PV label-distance keys support labels up to uint32_t max.");

    return (static_cast<LabelPairKey>(lhs) << 32) |
           static_cast<LabelPairKey>(rhs);
}

double computeLabelDistance(
    const PipelineRuntime& runtime,
    const std::vector<std::vector<size_t>>& pv_representative_lists,
    const size_t lhs,
    const size_t rhs)
{
    const auto& distance_matrix = runtime.distanceMatrix();
    const size_t original_vertex_count = distance_matrix.vertexCount();

    if (lhs < original_vertex_count && rhs < original_vertex_count)
        return distance_matrix.distance(lhs, rhs);

    double minimum_distance = std::numeric_limits<double>::infinity();

    if (lhs < original_vertex_count)
    {
        for (const size_t rhs_representative :
             pv_representative_lists[rhs - original_vertex_count])
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs, rhs_representative));
        }
        return minimum_distance;
    }

    if (rhs < original_vertex_count)
    {
        for (const size_t lhs_representative :
             pv_representative_lists[lhs - original_vertex_count])
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs_representative, rhs));
        }
        return minimum_distance;
    }

    for (const size_t lhs_representative :
         pv_representative_lists[lhs - original_vertex_count])
    {
        for (const size_t rhs_representative :
             pv_representative_lists[rhs - original_vertex_count])
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs_representative, rhs_representative));
        }
    }

    return minimum_distance;
}

struct QuotientEdgeData
{
    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    SimplexList sorted_edges;
};

QuotientEdgeData buildQuotientEdges(
    const PipelineRuntime& runtime,
    const std::vector<size_t>& active_labels,
    const std::vector<std::vector<size_t>>& pv_representative_lists,
    const double eps_hi)
{
    const int worker_count = runtime.config().threads;
    const size_t original_vertex_count = runtime.distanceMatrix().vertexCount();

    std::vector<robin_hood::unordered_map<LabelPairKey, double>>
        thread_pv_label_distance_hashes(static_cast<size_t>(worker_count));
    std::vector<SimplexList> thread_edge_workspaces(static_cast<size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (size_t i = 0; i < active_labels.size(); ++i)
    {
        const int thread_id = omp_get_thread_num();
        auto& thread_distances = thread_pv_label_distance_hashes[thread_id];
        auto& thread_edges = thread_edge_workspaces[thread_id];

        for (size_t j = i + 1; j < active_labels.size(); ++j)
        {
            size_t lhs = active_labels[i];
            size_t rhs = active_labels[j];
            if (lhs > rhs)
                std::swap(lhs, rhs);

            const double weight = computeLabelDistance(
                runtime,
                pv_representative_lists,
                lhs,
                rhs);

            if (lhs >= original_vertex_count || rhs >= original_vertex_count)
                thread_distances.emplace(makeLabelPairKey(lhs, rhs), weight);

            // Preserve quotient-edge semantics: zero-distance pairs are not edges.
            if (weight > 0.0 && weight < eps_hi)
            {
                const SimplexBindex bindex = SimplexUtility::getEdgeBinomialIndex(
                    runtime.binomialTable(),
                    rhs,
                    lhs);
                thread_edges.emplace_back(bindex, weight);
            }
        }
    }

    size_t pv_label_pair_count = 0;
    for (const auto& thread_distances : thread_pv_label_distance_hashes)
        pv_label_pair_count += thread_distances.size();

    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    pv_label_distances.reserve(pv_label_pair_count);
    for (const auto& thread_distances : thread_pv_label_distance_hashes)
    {
        pv_label_distances.insert(
            thread_distances.begin(),
            thread_distances.end());
    }

    auto sorted_edges = SimplexUtility::sortAndMergeSimplexChunks(
        thread_edge_workspaces,
        worker_count);

    return QuotientEdgeData{
        std::move(pv_label_distances),
        std::move(sorted_edges)};
}
}

WindowState::WindowState(const size_t original_vertex_count)
    : original_vertex_count_(original_vertex_count),
      active_label_list_(original_vertex_count)
{
    if (original_vertex_count_ == 0)
        throw std::invalid_argument("Window state requires at least one original vertex.");

    std::iota(active_label_list_.begin(), active_label_list_.end(), size_t{0});
}

SimplexList WindowState::prepareWindow(
    const PipelineRuntime& runtime,
    const WindowBounds& bounds,
    const WindowPreparationMode mode)
{
    if (prepared_)
        throw std::logic_error("The current window must be invalidated before preparing another window.");

    if (runtime.distanceMatrix().vertexCount() != original_vertex_count_)
        throw std::invalid_argument("Window state and pipeline runtime have different original vertex counts.");

    validateBounds(runtime.config(), bounds);
    validateActiveLabels(active_label_list_, totalLabelCount());

    if (runtime.binomialTable().size() <= totalLabelCount())
        throw std::logic_error("Pipeline binomial table capacity is too small for the current label range.");

    auto representative_lists = buildRepresentativeLists(
        pv_list_,
        original_vertex_count_);

    std::vector<uint8_t> active_label_mask(totalLabelCount(), 0);
    for (const size_t label : active_label_list_)
        active_label_mask[label] = 1;

    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    SimplexList sorted_edges;

    switch (mode)
    {
    case WindowPreparationMode::OrdinaryVr:
    {
        if (!pv_list_.empty() || active_label_list_.size() != original_vertex_count_)
            throw std::logic_error("Ordinary VR preparation requires the unchanged original label set.");

        for (size_t label = 0; label < active_label_list_.size(); ++label)
        {
            if (active_label_list_[label] != label)
                throw std::logic_error("Ordinary VR preparation requires every original label to be active.");
        }

        SimplexEnumerator simplex_enumerator(
            runtime.distanceMatrix(),
            runtime.binomialTable());
        // The ordinary enumerator intentionally retains zero-distance edges.
        sorted_edges = simplex_enumerator.getSortedVREdges(bounds.eps_hi);
        break;
    }

    case WindowPreparationMode::QuotientVr:
    {
        auto quotient_edges = buildQuotientEdges(
            runtime,
            active_label_list_,
            representative_lists,
            bounds.eps_hi);
        pv_label_distances = std::move(quotient_edges.pv_label_distances);
        sorted_edges = std::move(quotient_edges.sorted_edges);
        break;
    }

    default:
        throw std::invalid_argument("Unknown window preparation mode.");
    }

    active_label_mask_ = std::move(active_label_mask);
    pv_representative_lists_ = std::move(representative_lists);
    pv_label_distance_hash_ = std::move(pv_label_distances);
    bounds_ = bounds;
    preparation_mode_ = mode;
    prepared_ = true;

    return sorted_edges;
}

void WindowState::invalidateCurrentWindow() noexcept
{
    // These caches can be proportional to the complete active-label graph.
    // Empty swaps release their backing storage at the window boundary.
    decltype(active_label_mask_){}.swap(active_label_mask_);
    decltype(pv_representative_lists_){}.swap(pv_representative_lists_);
    decltype(pv_label_distance_hash_){}.swap(pv_label_distance_hash_);
    bounds_ = WindowBounds{};
    preparation_mode_ = WindowPreparationMode::OrdinaryVr;
    prepared_ = false;
}
