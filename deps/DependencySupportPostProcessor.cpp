#include "DependencySupportPostProcessor.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

void DependencySupportPostProcessor::processPwph(
    const PipelineRuntime& runtime,
    WindowState& window,
    DependencySupportBatch&& batch) const
{
    if (!window.prepared())
        throw std::logic_error("PwPH support processing requires a prepared window.");
    if (window.preparationMode() != WindowPreparationMode::QuotientVr)
        throw std::logic_error("PwPH support processing requires quotient-window geometry.");

    const std::size_t raw_support_count = batch.supports.size();
    const std::size_t original_vertex_count = window.originalVertexCount();
    const std::size_t total_label_count = window.totalLabelCount();

    for (const std::size_t label : batch.protected_labels)
    {
        if (label >= total_label_count)
            throw std::logic_error("Protected dependency-support label is outside the global label range.");
    }

    std::vector<std::unordered_set<std::size_t>> eligible_label_sets;
    eligible_label_sets.reserve(batch.supports.size());

    for (auto& support : batch.supports)
    {
        bool contains_pv = false;
        for (const std::size_t label : support.label_set)
        {
            if (label >= total_label_count)
                throw std::logic_error("Dependency support contains an invalid global label.");
            if (label >= original_vertex_count)
                contains_pv = true;
        }

        // Preserve the frozen-PV policy and its ordering. This test precedes
        // protection, so the additional protected PV labels retained by the
        // generic frame cannot affect current PwPH selection.
        if (contains_pv)
            continue;

        bool contains_protected_label = false;
        for (const std::size_t label : support.label_set)
        {
            if (batch.protected_labels.contains(label))
            {
                contains_protected_label = true;
                break;
            }
        }

        if (!contains_protected_label)
            eligible_label_sets.push_back(std::move(support.label_set));
    }

    if (runtime.config().verbose && runtime.config().maxdim >= 2)
    {
        std::cout << "raw pv support cofacets size = " << raw_support_count << '\n';
        std::cout << "pv support label sets with no pv/protected contents size = "
                  << eligible_label_sets.size() << '\n';
    }

    const double diameter_limit =
        runtime.config().pv_cap_scale * window.bounds().eps_hi;
    const double min_separation = runtime.config().absoluteMinSeparation();
    const auto& distance_matrix = runtime.distanceMatrix();

    std::vector<SelectedPV> selected_pvs;
    std::unordered_set<std::size_t> claimed_labels;

    const auto is_separated = [&](const std::unordered_set<std::size_t>& vertices)
    {
        if (min_separation <= 0.0)
            return true;

        // PwPH never deactivates an existing PV. Iterating active PV labels is
        // equivalent today and is the correct lifetime boundary for a future
        // absorbing policy, where historical inactive records must be ignored.
        for (const std::size_t label : window.active_label_list_)
        {
            if (label < original_vertex_count)
                continue;

            const auto& existing =
                window.pv_list_[label - original_vertex_count].flat_index_set;
            if (minCrossDistance(distance_matrix, vertices, existing) <
                min_separation)
            {
                return false;
            }
        }

        for (const auto& selected : selected_pvs)
        {
            if (minCrossDistance(
                    distance_matrix,
                    vertices,
                    selected.flat_index_set) < min_separation)
            {
                return false;
            }
        }

        return true;
    };

    // Matching discovers supports while traversing facets from high to low
    // rank. The existing PwPH policy then traverses that discovery vector in
    // reverse; this second reversal is intentionally preserved.
    for (auto iter = eligible_label_sets.rbegin();
         iter != eligible_label_sets.rend();
         ++iter)
    {
        const auto& labels = *iter;
        bool overlaps_claimed_label = false;
        for (const std::size_t label : labels)
        {
            if (claimed_labels.contains(label))
            {
                overlaps_claimed_label = true;
                break;
            }
        }
        if (overlaps_claimed_label)
            continue;

        auto flattened = flattenLabelSet(window, labels);
        if (flattened.empty() || flattened.size() >= MAX_PV_CARDINALITY)
            continue;

        const double diameter = maxPairwiseDistance(distance_matrix, flattened);
        if (diameter > diameter_limit)
            continue;
        if (!is_separated(flattened))
            continue;

        claimed_labels.insert(labels.begin(), labels.end());
        selected_pvs.push_back(SelectedPV{std::move(flattened), diameter});
    }

    commitSelectedPvs(
        window,
        std::move(selected_pvs),
        claimed_labels);
}

std::unordered_set<std::size_t>
DependencySupportPostProcessor::flattenLabelSet(
    const WindowState& window,
    const std::unordered_set<std::size_t>& labels)
{
    std::unordered_set<std::size_t> flattened;
    for (const std::size_t label : labels)
    {
        if (label < window.original_vertex_count_)
        {
            flattened.insert(label);
            continue;
        }

        const std::size_t pv_index = label - window.original_vertex_count_;
        if (pv_index >= window.pv_list_.size())
            throw std::logic_error("Cannot flatten an invalid PV label.");

        const auto& representatives = window.pv_list_[pv_index].flat_index_set;
        flattened.insert(representatives.begin(), representatives.end());
    }
    return flattened;
}

double DependencySupportPostProcessor::maxPairwiseDistance(
    const DistanceMatrix& distance_matrix,
    const std::unordered_set<std::size_t>& vertices)
{
    double maximum = 0.0;
    for (auto lhs = vertices.begin(); lhs != vertices.end(); ++lhs)
    {
        auto rhs = lhs;
        ++rhs;
        for (; rhs != vertices.end(); ++rhs)
            maximum = std::max(maximum, distance_matrix.distance(*lhs, *rhs));
    }
    return maximum;
}

double DependencySupportPostProcessor::minCrossDistance(
    const DistanceMatrix& distance_matrix,
    const std::unordered_set<std::size_t>& lhs,
    const std::unordered_set<std::size_t>& rhs)
{
    double minimum = std::numeric_limits<double>::infinity();
    for (const std::size_t lhs_vertex : lhs)
    {
        for (const std::size_t rhs_vertex : rhs)
        {
            minimum = std::min(
                minimum,
                distance_matrix.distance(lhs_vertex, rhs_vertex));
        }
    }
    return minimum;
}

void DependencySupportPostProcessor::commitSelectedPvs(
    WindowState& window,
    std::vector<SelectedPV>&& selected_pvs,
    const std::unordered_set<std::size_t>& absorbed_labels)
{
    std::vector<std::size_t> next_active_labels;
    next_active_labels.reserve(
        window.active_label_list_.size() + selected_pvs.size());

    for (const std::size_t label : window.active_label_list_)
    {
        if (!absorbed_labels.contains(label))
            next_active_labels.push_back(label);
    }

    std::size_t next_pv_label = window.totalLabelCount();
    window.pv_list_.reserve(window.pv_list_.size() + selected_pvs.size());
    for (auto& selected : selected_pvs)
    {
        window.pv_list_.push_back(std::move(selected));
        next_active_labels.push_back(next_pv_label++);
    }

    std::sort(next_active_labels.begin(), next_active_labels.end());
    if (std::adjacent_find(
            next_active_labels.begin(),
            next_active_labels.end()) != next_active_labels.end())
    {
        throw std::logic_error("PV commit produced duplicate active labels.");
    }

    window.active_label_list_ = std::move(next_active_labels);

    // Window-scoped caches are never carried across a completed transition,
    // including transitions that selected no PVs.
    window.invalidateCurrentWindow();
}
