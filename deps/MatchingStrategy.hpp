#pragma once

#include <vector>

#include "MatchingContext.hpp"


class MatchingStrategy
{
public:
    virtual ~MatchingStrategy() = default;

    virtual void match(MatchingContext& matching_context) = 0;

};