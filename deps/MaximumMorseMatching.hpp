#pragma once

#include <functional>
#include <queue>

#include "MatchingStrategy.hpp"
#include "BipartiteGraph.hpp"

class MaximumMorseMatching : public MatchingStrategy
{
public:

    MaximumMorseMatching() {};

    size_t implicitMatch(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    int64_t implicitMatchAndGetMinCritialIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);    //for the quotient

    std::vector<std::vector<size_t>> implicitMatchAndCollectPVSupports(MatchingContext& matching_context);

    //legacy explicit graph representation  
    MaximumMorseMatching(int threadnum): threadnum_(threadnum) {};

    size_t match(MatchingContext& matching_context) override;

    int64_t matchAndGetMinCriticalIndex(MatchingContext &matching_context);

    std::vector<std::vector<size_t>> matchAndGetAugPath(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

private:
    //helper strcut to allow moving the underlying container from a priority_queue
    template <class StdComparator>
    struct QueueWorkspaceHelper : std::priority_queue<size_t, std::vector<size_t>, StdComparator> 
    {
        using Base = std::priority_queue<size_t, std::vector<size_t>, StdComparator>;
        using Base::Base;    //inherit base constructors

        //verbose form using std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>>::priority_queue;

        std::vector<size_t>& getContainer() { return this->c; }
    };

    using MinIndexQueue = QueueWorkspaceHelper<std::greater<size_t>>;
    using MaxIndexQueue = QueueWorkspaceHelper<std::less<size_t>>;

    //workspaces for matching
    std::vector<size_t> aug_path_;
    std::vector<size_t> cofacet_stack_;
    std::vector<size_t> cofacet_indices_;
    std::vector<size_t> facet_indices_;
    std::vector<size_t> vertex_workspace_;
    std::vector<size_t> pq_workspace_;

    int64_t implicitFacetAugPath(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& facet_list, 
                                 const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, size_t npts, size_t interfacedimension);

    std::vector<size_t> collectReducedColumnSupport(const MatchingContext& matching_context, size_t terminalcofacet);
    

    //legacy explicit graph representation  
    int threadnum_ = 4;    //for explicit initialization phase

    void parallelMinCofacetInit(MatchingContext &matching_context);

    int64_t serialFacetAugPath(BipartiteGraph &graph, const size_t facetindex, std::vector<size_t> &aug_path, std::vector<size_t> &cofacet_stack);

    size_t serialFacetMatch(MatchingContext &matching_context);

    int64_t serialFacetMatchAndGetMinCriticalIndex(MatchingContext &matching_context);

    std::vector<std::vector<size_t>> serialFacetMatchAndGetAugPath(MatchingContext &matching_context);
};