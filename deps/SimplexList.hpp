#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "FiltrationValueType.hpp"

struct SimplexList
{
    struct EntryRef
    {
        int64_t& first;
        FiltrationValueType& second;
    };

    struct ConstEntryRef
    {
        const int64_t& first;
        const FiltrationValueType& second;
    };

    std::vector<int64_t> bindices;
    std::vector<FiltrationValueType> weights;

    size_t size() const { return bindices.size(); }
    bool empty() const { return bindices.empty(); }

    void reserve(const size_t n)
    {
        bindices.reserve(n);
        weights.reserve(n);
    }

    void clear()
    {
        bindices.clear();
        weights.clear();
    }

    void emplace_back(const int64_t bindex, const FiltrationValueType weight)
    {
        bindices.push_back(bindex);
        weights.push_back(weight);
    }

    EntryRef operator[](const size_t i) { return EntryRef{bindices[i], weights[i]}; }
    ConstEntryRef operator[](const size_t i) const { return ConstEntryRef{bindices[i], weights[i]}; }

    void swap(SimplexList& other) noexcept
    {
        bindices.swap(other.bindices);
        weights.swap(other.weights);
    }
};

inline void swap(SimplexList& lhs, SimplexList& rhs) noexcept
{
    lhs.swap(rhs);
}
