#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "PipelineCommon.hpp"
#include "robin_hood.h"

struct SelectedPV
{
    // Canonical sorted original-vertex indices represented by this PV.
    std::vector<size_t> representatives;
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

    SimplexList prepareWindow(
        const PipelineRuntime& runtime,
        const WindowBounds& bounds);

    void commitSelectedPVs(
        std::vector<SelectedPV>&& selected_pvs,
        const std::unordered_set<size_t>& new_absorbed_labels);

private:
    void invalidateCurrentWindow() noexcept;

    size_t original_vertex_count_ = 0;

    // Append-only for the complete pipeline run.
    std::vector<SelectedPV> pv_list_;

    // Sorted labels active in the current quotient complex.
    std::vector<size_t> active_label_list_;

    // Current-window caches, all defined over the complete historical label range.
    std::vector<uint8_t> active_label_mask_;
    robin_hood::unordered_map<LabelPairKey, double> pv_label_distance_hash_;

    WindowBounds bounds_;
};
