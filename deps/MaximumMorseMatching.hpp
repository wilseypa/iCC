#pragma once

#include <cstdint>
#include <functional>
#include <queue>

#include "MatchingStrategy.hpp"
#include "BipartiteGraph.hpp"

class MaximumMorseMatching : public MatchingStrategy
{
public:
    struct PersistentPairInfo
    {
        double facet_weight = -1.0;
        double cofacet_weight = -1.0;
        int64_t facet_bindex = -1;
        int64_t cofacet_bindex = -1;
    };

    struct MatchSupportInfo
    {
        std::vector<std::vector<size_t>> raw_pv_support_cofacet_indices;
        std::vector<size_t> protected_facet_list_indices;
    };

    MaximumMorseMatching() {};

    size_t implicitMatch(MatchingContext& matching_context,
                         std::vector<std::pair<double, double>>& dim_persistent_pair,
                         std::vector<PersistentPairInfo>* persistent_pair_info = nullptr);


    MatchSupportInfo implicitMatchAndCollectSupportInfo(MatchingContext& matching_context,
                                                        std::vector<std::pair<double, double>>& dim_persistent_pair,
                                                        const bool collect_pv_support,
                                                        std::vector<PersistentPairInfo>* persistent_pair_info = nullptr);

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

    int64_t implicitFacetAugPath(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, 
                                 const std::vector<std::pair<int64_t, double>>& facet_list, 
                                 const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, 
                                 size_t npts, size_t interfacedimension);

    int64_t implicitFacetCompressedAugPath(const std::vector<std::vector<int64_t>>& binomial_table, const BipartiteGraph& bi_graph,
                                           const std::vector<std::pair<int64_t, double>>& facet_list,
                                           const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table,
                                           const size_t facetgraphindex, size_t npts, size_t interfacedimension);

    void enqueueReducedCompressedColumnTail(const std::vector<std::vector<int64_t>>& binomial_table,
                                            const BipartiteGraph& bi_graph,
                                            const std::vector<std::pair<int64_t, double>>& facet_list,
                                            const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table,
                                            const size_t facet_list_index,
                                            const size_t pivot_cofacet,
                                            size_t npts,
                                            size_t interfacedimension,
                                            MinIndexQueue& target_queue);

    std::vector<size_t> collectReducedColumnSupport(const MatchingContext& matching_context, const size_t terminalcofacet, 
                                                    const size_t expectedfacet);

    void enqueueReducedCofacetBoundaryTail(const MatchingContext& matching_context,
                                           const size_t cofacet_list_index,
                                           const size_t pivot_facet_list_index,
                                           MaxIndexQueue& target_queue,
                                           std::vector<size_t>& cofacet_trace);
    

    int64_t implicitFacetAugPathDebug(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& facet_list, 
                                 const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, size_t npts, size_t interfacedimension);

    //legacy explicit graph representation  
    int threadnum_ = 4;    //for explicit initialization phase

    void parallelMinCofacetInit(MatchingContext &matching_context);

    int64_t serialFacetAugPath(BipartiteGraph &graph, const size_t facetindex, std::vector<size_t> &aug_path, std::vector<size_t> &cofacet_stack);

    size_t serialFacetMatch(MatchingContext &matching_context);

    int64_t serialFacetMatchAndGetMinCriticalIndex(MatchingContext &matching_context);

    std::vector<std::vector<size_t>> serialFacetMatchAndGetAugPath(MatchingContext &matching_context);
};
