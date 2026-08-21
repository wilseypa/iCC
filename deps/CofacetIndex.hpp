#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "robin_hood.h"

/**
 * Memory-compact lookup from a cofacet's binomial index to its position in the
 * filtration-ordered cofacet list.
 *
 * The set stores only 32-bit list ranks. Its transparent hash and equality
 * functors obtain the corresponding binomial index from the external cofacet
 * list, so that lookups can still be performed directly with an int64_t
 * binomial index. The referenced list and its element storage must remain
 * unchanged for the lifetime of this index.
 */
class CofacetIndex
{
public:
    using Rank = uint32_t;
    using Cofacet = std::pair<int64_t, double>;
    using CofacetList = std::vector<Cofacet>;

private:
    struct IndirectHash
    {
        using is_transparent = void;

        const Cofacet* cofacets = nullptr;

        size_t operator()(const Rank rank) const noexcept
        {
            return robin_hood::hash<int64_t>{}(cofacets[rank].first);
        }

        size_t operator()(const int64_t bindex) const noexcept
        {
            return robin_hood::hash<int64_t>{}(bindex);
        }
    };

    struct IndirectEqual
    {
        using is_transparent = void;

        const Cofacet* cofacets = nullptr;

        bool operator()(const Rank lhs, const Rank rhs) const noexcept
        {
            return cofacets[lhs].first == cofacets[rhs].first;
        }

        bool operator()(const int64_t bindex, const Rank rank) const noexcept
        {
            return bindex == cofacets[rank].first;
        }
    };

    using RankSet = robin_hood::unordered_flat_set<Rank, IndirectHash, IndirectEqual>;

    struct State
    {
        explicit State(const CofacetList& cofacets)
            : ranks(0, IndirectHash{cofacets.data()}, IndirectEqual{cofacets.data()})
        {
            constexpr uint64_t max_rank_count =
                static_cast<uint64_t>(std::numeric_limits<Rank>::max()) + 1ULL;
            if (static_cast<uint64_t>(cofacets.size()) > max_rank_count)
                throw std::length_error("CofacetIndex rank exceeds uint32_t capacity");

            ranks.reserve(cofacets.size());
            for (size_t rank = 0; rank < cofacets.size(); ++rank)
            {
                const auto result = ranks.emplace(static_cast<Rank>(rank));
                if (!result.second)
                    throw std::logic_error("CofacetIndex requires unique binomial indices");
            }
        }

        RankSet ranks;
    };

public:
    CofacetIndex() = default;

    explicit CofacetIndex(const CofacetList& cofacets)
    {
        if (!cofacets.empty())
            state_ = std::make_unique<State>(cofacets);
    }

    CofacetIndex(CofacetList&&) = delete;
    CofacetIndex(const CofacetList&&) = delete;

    CofacetIndex(const CofacetIndex&) = delete;
    CofacetIndex& operator=(const CofacetIndex&) = delete;
    CofacetIndex(CofacetIndex&&) noexcept = default;
    CofacetIndex& operator=(CofacetIndex&&) noexcept = default;

    [[nodiscard]] std::optional<Rank> findRank(const int64_t bindex) const noexcept
    {
        if (state_ == nullptr || state_->ranks.empty())
            return std::nullopt;

        const auto it = state_->ranks.find(bindex);
        return it == state_->ranks.end() ? std::nullopt : std::optional<Rank>(*it);
    }

    [[nodiscard]] size_t size() const noexcept
    {
        return state_ == nullptr ? 0 : state_->ranks.size();
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return state_ == nullptr || state_->ranks.empty();
    }

    void release() noexcept
    {
        state_.reset();
    }

private:
    std::unique_ptr<State> state_;
};
