#pragma once

#include "MatchingContext.hpp"
#include "MatchingStrategy.hpp"

class ApproximateMorseMatching : public MatchingStrategy
{
public:
    ApproximateMorseMatching(int threadnum) : threadnum_(threadnum) {};
    size_t match(MatchingContext& matching_context) override;

private:
    int threadnum_;

    void parallelTwoPhaseInit(BipartiteGraph& graph);

    void parallelFacetDFSMatch(BipartiteGraph& graph);

    int64_t facetDfsAugPath(BipartiteGraph& graph, const size_t startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<size_t>& aug_path_tid);

    void dfsCycleRemoval(BipartiteGraph& graph);
};
