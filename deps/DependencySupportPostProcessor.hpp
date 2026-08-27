#pragma once

#include <cstddef>
#include <unordered_set>
#include <vector>

#include "DependencySupport.hpp"
#include "PipelineCommon.hpp"
#include "WindowState.hpp"

class DependencySupportPostProcessor
{
public:
    // Current frozen-PV policy. A future absorbing policy belongs beside this
    // method and can reuse flattenLabelSet()/commitSelectedPvs() without any
    // matching, enumeration, or WindowState representation change.
    void processPwph(
        const PipelineRuntime& runtime,
        WindowState& window,
        DependencySupportBatch&& batch) const;

private:
    static constexpr std::size_t MAX_PV_CARDINALITY = 64;

    static std::unordered_set<std::size_t> flattenLabelSet(
        const WindowState& window,
        const std::unordered_set<std::size_t>& labels);

    static double maxPairwiseDistance(
        const DistanceMatrix& distance_matrix,
        const std::unordered_set<std::size_t>& vertices);

    static double minCrossDistance(
        const DistanceMatrix& distance_matrix,
        const std::unordered_set<std::size_t>& lhs,
        const std::unordered_set<std::size_t>& rhs);

    static void commitSelectedPvs(
        WindowState& window,
        std::vector<SelectedPV>&& selected_pvs,
        const std::unordered_set<std::size_t>& absorbed_labels);
};
