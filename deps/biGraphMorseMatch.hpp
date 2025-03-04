#pragma once

#include <vector>
#include <unordered_set>


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


    // bool isRightSinglePath(int cofacetindex);

    //parallel bipartite matching helper function
    int64_t facetDfsAugPath(const size_t startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<size_t>& aug_path_tid);

    int64_t facetRightDfsAugPath(const size_t startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<size_t>& aug_path_tid);

    void add2SingleOrRemove(const size_t index, std::vector<uint64_t>& removed_flag, const size_t round);
    int64_t serialCofacetDFSAugPath(const size_t cofacetindex, std::vector<uint8_t>& cofacet_dfs_flag, std::vector<size_t>& aug_path, std::vector<uint64_t>& removed_flag, const size_t round);

    // std::set<int> criticalFacetBackwardSearch(int facetindex);
    

    //cycle removal helper function
    // int findRoot();
    // int getParent(std::vector<int>& parent_workspace, int uidx);

    size_t getChild(std::vector<size_t>& child_workspace, const size_t uidx);

    // void getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx);
    // bool isBackwardAcyclic(std::vector<int>& ancestor_simp, std::vector<int>& u_child);
    // int lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_simp, std::vector<int>& child_simp, int rootnum, int uidx);
    // int serialBFS(std::vector<int>& ancestor_simp, int rootnum, int root);

    void parallelFacetDFSMatch(const int threadnum);

    void parallelDirectionalFacetDFSMatch(const int threadnum);

    void serialCofacetDFSMatch();

    // int serialCycleRemoval();

    int dfsCycleRemoval();

    void addEdge(int u, int v);

    Bi_Graph_Match(size_t leftnum, size_t rightnum, size_t leftdimension);

    void updateDimension(size_t newleftnum, size_t newrightnum);

    void buildInterface(const std::vector<std::vector<int>>& cofacet_bin, int cofacet_index_min, int cofacet_index_max, const std::vector<std::vector<int>>& simplex_bin, int simplex_index_max, const std::vector<int>& active_index);

    std::unordered_set<size_t> getActiveIndexSet(); 
    
    std::vector<size_t> getCriticalIndex(const std::unordered_set<size_t>& dim_active_index_set, const size_t simplex_index_max);

    // std::vector<std::set<int>> getBackwardSingleFacetIndex(std::vector<int> critical_facet_index);


    void checkSimplex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_simplex);

    void checkCofacet(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_cofacet);

    void checkSimplexByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_simplex_index);

    void checkCofacetByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_cofacet_index);

};