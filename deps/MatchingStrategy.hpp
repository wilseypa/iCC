#pragma once

#include <vector>

#include "MatchingContext.hpp"

class MatchingStrategy
{
public:
    virtual ~MatchingStrategy() = default;

    virtual size_t match(MatchingContext& matching_context) = 0;
};