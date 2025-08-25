#pragma once

#include "MatchingStrategy.hpp"
#include "BipartiteGraph.hpp"

class MaximumMorseMatching : public MatchingStrategy
{
public:

    explicit MaximumMorseMatching(int threadnum): threadnum_(threadnum) {};

    size_t match(MatchingContext& matching_context) override;

    size_t matchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    size_t matchWithPersistenceBackup(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    int64_t matchWithPersistenceAndReturnMinCriticalIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

private:

    int threadnum_;    //for initialization phase

    void parallelMaxFacetInit(MatchingContext& matching_context);

    void parallelMinCofacetInit(MatchingContext& matching_context);

    int64_t serialCofacetAugPath(BipartiteGraph& graph, const size_t cofacetindex, std::vector<size_t>& aug_path, std::vector<size_t>& facet_stack);
    
    int64_t serialFacetAugPath(BipartiteGraph& graph, const size_t facetindex, std::vector<size_t>& aug_path, std::vector<size_t>& cofacet_stack);

    size_t serialCofacetMorseMatch(BipartiteGraph& graph);

    size_t serialFacetMorseMatch(BipartiteGraph& graph);

    size_t serialCofacetMorseMatchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    size_t serialFacetMorseMatchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    int64_t serialFacetMorseMatchAndPersistenceReturnMinCriticalIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

};