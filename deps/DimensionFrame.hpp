#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DependencySupport.hpp"
#include "PipelineCommon.hpp"
#include "SimplexEnumerator.hpp"
#include "WindowState.hpp"
#include "robin_hood.h"

class ImplicitMorseMatching;

enum class FramePhase : std::uint8_t
{
    ReadyToMatch,
    Matched,
    Finished
};

class DimensionFrame
{
public:
    DimensionFrame(
        const PipelineRuntime& runtime,
        const WindowState& window,
        SimplexList prepared_edges);

    DimensionFrame(const DimensionFrame&) = delete;
    DimensionFrame& operator=(const DimensionFrame&) = delete;
    DimensionFrame(DimensionFrame&&) = delete;
    DimensionFrame& operator=(DimensionFrame&&) = delete;

    void matchPersistence(bool collect_dependency_support);
    void matchApparentPairsOnly();

    [[nodiscard]]
    bool advance();

    [[nodiscard]]
    std::size_t dimension() const noexcept
    {
        return dimension_;
    }

    [[nodiscard]]
    std::size_t facetDimension() const noexcept
    {
        return dimension_ == 1 ? 1 : dimension_ - 1;
    }

    [[nodiscard]]
    bool isTopDimension() const noexcept
    {
        return dimension_ == runtime_.config().maxdim;
    }

    [[nodiscard]]
    FramePhase phase() const noexcept
    {
        return phase_;
    }

    [[nodiscard]]
    std::size_t facetCount() const noexcept
    {
        return facet_list_.size();
    }

    [[nodiscard]]
    std::size_t cofacetCount() const noexcept
    {
        return cofacet_list_.size();
    }

    [[nodiscard]]
    const std::vector<PersistentPairInfo>& persistentPairs() const noexcept
    {
        return persistent_pairs_;
    }

    [[nodiscard]]
    const std::vector<double>& unmatchedFacetWeights() const noexcept
    {
        return unmatched_facet_weights_;
    }

    [[nodiscard]]
    const std::vector<double>& unmatchedTopCofacetWeights() const noexcept
    {
        return unmatched_top_cofacet_weights_;
    }

    DependencySupportBatch takeDependencySupportBatch();

private:
    friend class ImplicitMorseMatching;

    class FullCofacetHash
    {
    public:
        explicit FullCofacetHash(const SimplexList& cofacets);

        [[nodiscard]]
        bool tryFindRank(
            SimplexBindex bindex,
            std::size_t& rank_out) const noexcept;

    private:
        robin_hood::unordered_map<SimplexBindex, std::size_t> hash_;
    };

    struct RawMatchSupportInfo
    {
        std::vector<std::vector<std::size_t>> support_cofacet_indices;
        std::vector<std::size_t> protected_facet_indices;
    };

    struct InterfaceWorkspace
    {
        explicit InterfaceWorkspace(const SimplexList& cofacets, std::size_t facet_count);

        FullCofacetHash cofacet_lookup;
        std::vector<std::int64_t> match_list;

        std::vector<std::size_t> cofacet_indices;
        std::vector<std::size_t> facet_indices;
        std::vector<std::size_t> vertex_workspace;
        std::vector<std::size_t> priority_queue_workspace;

        RawMatchSupportInfo raw_support_info;
    };

    struct FfiRealizationState
    {
        robin_hood::unordered_map<SimplexBindex, std::uint64_t> facet_realizations;
        robin_hood::unordered_map<SimplexBindex, std::uint64_t> cofacet_realizations;
    };

    const PipelineRuntime& runtime_;
    const WindowState& window_;
    const BinomialTable& binomial_table_;

    std::size_t dimension_ = 2;
    FramePhase phase_ = FramePhase::ReadyToMatch;

    SimplexEnumerator simplex_enumerator_;

    SimplexList facet_list_;
    SimplexList cofacet_list_;

    robin_hood::unordered_map<SimplexBindex, std::size_t> active_facet_index_;

    std::vector<std::uint8_t> next_active_cofacet_mask_;
    std::size_t next_active_cofacet_count_ = 0;

    std::vector<PersistentPairInfo> persistent_pairs_;
    std::vector<double> unmatched_facet_weights_;
    std::vector<double> unmatched_top_cofacet_weights_;

    DependencySupportBatch dependency_support_batch_;
    std::unordered_set<std::size_t> protected_labels_;

    std::optional<FfiRealizationState> ffi_realizations_;
    std::optional<InterfaceWorkspace> interface_workspace_;

    void requirePhase(FramePhase expected, const char* operation) const;
    void buildInitialActiveFacetIndex();
    void enumerateInitialCofacets();
    SimplexList enumerateNextCofacets(
        robin_hood::unordered_map<SimplexBindex, std::uint64_t>* realization_output);

    void materializeMatchSupport();
    void finalizeNextActiveCofacets();
    void releaseInterfaceWorkspace();
    void clearCurrentInterfaceResults();

    [[nodiscard]]
    bool tryFindCofacetRank(
        SimplexBindex bindex,
        std::size_t& rank_out) const noexcept;

    void reportFalseFacetIdentificationStats() const;
};
