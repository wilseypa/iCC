#pragma once

#include <vector>
#include <unordered_set>

#include"robin_hood.h"



class Bi_Graph_Match 
{
public:
    size_t u, v;

    std::vector<std::vector<size_t>> adj_list;
    std::vector<int64_t> match_list;
    
    size_t udegree;

    // int maxthreadnum;    //affect the memory consumption during bipartite matching

    //parallel init for bipartite matching
    
    void pairDegreeOne(const size_t uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg);
    void pairUnmatched(const size_t uidx, std::vector<int>& visit_flag);
    void parallelKarpSipserInit(const int threadnumber);

    // void elementaryCollapse(int vidx, int umin, int vmin, std::vector<int>& visit_flag_u, std::vector<int>& visit_flag_v, std::vector<int>& node_deg);

    void parallelMaxFacetInit(const size_t cofacet_index_min, const size_t cofacet_index_max, const size_t facet_index_min, const size_t facet_index_max, const int threadnum);

    void parallelTwoPhaseInit(const int threadnum);

    void parallelMinCofacetInit(const size_t cofacet_index_min, const size_t cofacet_index_max, const size_t facet_index_min, const size_t facet_index_max, const int threadnum);

    int64_t facetDfsAugPath(const size_t startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<size_t>& aug_path_tid);

    int64_t serialCofacetDFSAugPath(const size_t cofacetindex, std::vector<size_t>& aug_path, std::vector<size_t>& facet_stack);
    
    int64_t serialFacetDFSAugPath(const size_t facetindex, std::vector<size_t>& aug_path, std::vector<size_t>& cofacet_stack);

    size_t getChild(std::vector<size_t>& child_workspace, const size_t uidx);

    void parallelFacetDFSMatch(const int threadnum);


    std::vector<std::pair<double, double>> serialCofacetDFSMatch(const std::vector<std::pair<int64_t, double>>& sorted_facet, const std::vector<std::pair<int64_t, double>>& sorted_cofacet);

    std::vector<std::pair<double, double>> serialFacetDFSMatch(const std::vector<std::pair<int64_t, double>>& sorted_facet, const std::vector<std::pair<int64_t, double>>& sorted_cofacet);

    int dfsCycleRemoval();

    // void addEdge(int u, int v);

    Bi_Graph_Match(size_t leftnum, size_t rightnum, size_t leftdimension);

    void updateDimension(size_t newleftnum, size_t newrightnum);

    // void buildInterface(const std::vector<std::vector<int>>& cofacet_bin, int cofacet_index_min, int cofacet_index_max, const std::vector<std::vector<int>>& simplex_bin, int simplex_index_max, const std::vector<int>& active_index);

    std::unordered_set<size_t> getActiveIndexSet(); 
    
    std::vector<size_t> getCriticalIndex(const std::unordered_set<size_t>& dim_active_index_set, const size_t simplex_index_max);

    robin_hood::unordered_map<int64_t, size_t> getActiveIndexHashTable(const std::vector<std::pair<int64_t, double>>& cofacet_list);


    void checkSimplex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_simplex);

    void checkCofacet(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_cofacet);

    void checkSimplexByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_simplex_index);

    void checkCofacetByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_cofacet_index);


    std::vector<std::vector<size_t>> serialCofacetDFSReduction(const int maxdim, int dim);

};