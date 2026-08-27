#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "PipelineCommon.hpp"
#include "robin_hood.h"

class DependencySupportPostProcessor;

struct SelectedPV
{
    // PVs are always flattened to original-vertex indices.
    std::unordered_set<size_t> flat_index_set;
    double diameter = 0.0;
};

class WindowState
{
public:
    explicit WindowState(size_t original_vertex_count);

    [[nodiscard]]
    size_t originalVertexCount() const noexcept
    {
        return original_vertex_count_;
    }

    [[nodiscard]]
    size_t totalLabelCount() const noexcept
    {
        return original_vertex_count_ + pv_list_.size();
    }

    [[nodiscard]]
    const std::vector<SelectedPV>& pseudoVertices() const noexcept
    {
        return pv_list_;
    }

    [[nodiscard]]
    const std::vector<size_t>& activeLabels() const noexcept
    {
        return active_label_list_;
    }

    [[nodiscard]]
    const std::vector<uint8_t>& activeLabelMask() const noexcept
    {
        return active_label_mask_;
    }

    [[nodiscard]]
    const std::vector<std::vector<size_t>>& pvRepresentativeLists() const noexcept
    {
        return pv_representative_lists_;
    }

    [[nodiscard]]
    const robin_hood::unordered_map<LabelPairKey, double>&
    pvLabelDistanceHash() const noexcept
    {
        return pv_label_distance_hash_;
    }

    [[nodiscard]]
    const WindowBounds& bounds() const noexcept
    {
        return bounds_;
    }

    [[nodiscard]]
    bool prepared() const noexcept
    {
        return prepared_;
    }

    [[nodiscard]]
    WindowPreparationMode mode() const noexcept
    {
        return preparation_mode_;
    }

    // Descriptive alias used by dimension-level orchestration.
    [[nodiscard]]
    WindowPreparationMode preparationMode() const noexcept
    {
        return mode();
    }

    SimplexList prepareWindow(
        const PipelineRuntime& runtime,
        const WindowBounds& bounds,
        WindowPreparationMode mode);

private:
    friend class DependencySupportPostProcessor;

    void invalidateCurrentWindow() noexcept;

    size_t original_vertex_count_ = 0;

    // Append-only for the complete pipeline run.
    std::vector<SelectedPV> pv_list_;

    // Sorted labels active in the current quotient complex.
    std::vector<size_t> active_label_list_;

    // Current-window caches, all defined over the complete historical label range.
    std::vector<uint8_t> active_label_mask_;
    std::vector<std::vector<size_t>> pv_representative_lists_;
    robin_hood::unordered_map<LabelPairKey, double> pv_label_distance_hash_;

    WindowBounds bounds_;
    WindowPreparationMode preparation_mode_ = WindowPreparationMode::OrdinaryVr;
    bool prepared_ = false;
};
