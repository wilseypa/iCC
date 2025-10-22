#pragma once

#include "MatchingStrategy.hpp"
#include "BipartiteGraph.hpp"

class MaximumMorseMatching : public MatchingStrategy
{
public:
    MaximumMorseMatching(int threadnum) : threadnum_(threadnum) {};

    MaximumMorseMatching() {};

    size_t implicitMatch(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);

    int64_t implicitMatchReturnMinCritialIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair);    //for the quotient

    //legacy explicit graph representation  
    MaximumMorseMatching(int threadnum): threadnum_(threadnum) {};

    size_t matchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    size_t matchWithPersistenceBackup(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    int64_t matchWithPersistenceReturnMinCriticalIndex(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    std::vector<std::vector<size_t>> matchWithPersistenceReturnAugPath(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

private:
    int threadnum_; // for initialization phase

    //helper strcut to allow moving the underlying container from a priority_queue
    struct QueueWorkspaceHelper : std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> 
    {
        using Base = std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>>;
        using Base::Base;    //inherit base constructors

        //verbose form using std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>>::priority_queue;

        std::vector<size_t>& getContainer() { return this->c; }
    };


    //workspaces for matching
    std::vector<size_t> aug_path_;
    std::vector<size_t> cofacet_stack_;
    std::vector<size_t> cofacet_indices_;
    std::vector<size_t> facet_indices_;
    std::vector<size_t> vertex_workspace_;
    std::vector<size_t> pq_workspace_;

    int64_t implicitFacetAugPath(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& facet_list, 
                                 const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, size_t npts, size_t interfacedimension);

    //legacy explicit graph representation  
    int threadnum_ = 4;    //for explicit initialization phase

    void parallelMinCofacetInit(MatchingContext &matching_context);

    int64_t serialCofacetAugPath(BipartiteGraph &graph, const size_t cofacetindex, std::vector<size_t> &aug_path, std::vector<size_t> &facet_stack);

    int64_t serialFacetAugPath(BipartiteGraph &graph, const size_t facetindex, std::vector<size_t> &aug_path, std::vector<size_t> &cofacet_stack);

    int64_t serialFacetAugPathTest(BipartiteGraph& graph, const size_t facetindex, std::vector<size_t>& aug_path, std::vector<size_t>& cofacet_stack);

    size_t serialCofacetMatch(BipartiteGraph& graph);

    size_t serialFacetMatch(BipartiteGraph &graph);

    size_t serialCofacetMatchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    size_t serialFacetMatchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    int64_t serialFacetMatchWithPersistenceReturnMinCriticalIndex(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);

    std::vector<std::vector<size_t>> serialFacetMatchWithPersistenceReturnAugPath(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair);
};