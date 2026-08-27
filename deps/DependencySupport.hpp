#pragma once

#include <cstddef>
#include <unordered_set>
#include <vector>

struct DependencySupport
{
    // Labels are the current global original-or-PV labels, never list ranks.
    std::unordered_set<std::size_t> label_set;
};

struct DependencySupportBatch
{
    // Matching appends supports while active facets are discovered in descending
    // filtration rank. PwPH consumes this vector from back to front; preserving
    // both directions is part of the existing overlap-selection semantics.
    std::vector<DependencySupport> supports;

    // Every current global label belonging to a protected critical simplex.
    // This deliberately includes PV labels. PwPH rejects PV-containing supports
    // before consulting this set, so retaining those labels is behavior-neutral.
    std::unordered_set<std::size_t> protected_labels;
};
