#pragma once

#include <vector>
#include <deque>

class Bi_Graph_Match {
private:
    int u, v;
    
    int udegree;
    int maxthreadnum;    //affect the memory consumption during bipartite matching

    // std::vector<std::vector<int>> adj_list;
    // std::vector<int> match_list;
    std::vector<int> state_flag;

    //parallel init for bipartite matching
    
    void pairDegreeOne(int uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg);
    void pairUnmatched(int uidx, std::vector<int>& visit_flag);

    void elementaryCollapse(int vidx, int umin, int vmin, std::vector<int>& visit_flag_u, std::vector<int>& visit_flag_v, std::vector<int>& node_deg);

    bool isRightSinglePath(int cofacetindex);

    //parallel bipartite matching helper function
    int facetDfsAugPath(int startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    int facetRightSinglePathDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    int facetRightDFSAugPathWithCheck(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    int facetLeftDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    int facetRightDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    bool add2SingleOrRemove(int index, std::set<int>& single_index, std::set<int>& removed_index);
    int serialCofacetLeftDFSAugPath(int cofacetindex, std::vector<int>& cofacet_dfs_flag, std::vector<int>& aug_path);

    std::set<int> criticalFacetBackwardSearch(int facetindex);
    

    //cycle removal helper function
    int findRoot();
    int getParent(std::vector<int>& parent_workspace, int uidx);

    int getChild(std::vector<int>& child_workspace, int uidx);

    void getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx);
    bool isBackwardAcyclic(std::vector<int>& ancestor_simp, std::vector<int>& u_child);
    int lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_simp, std::vector<int>& child_simp, int rootnum, int uidx);
    int serialBFS(std::vector<int>& ancestor_simp, int rootnum, int root);

public:
    void parallelKarpSipserInit();

    void parallelMaxFacetInit(int cofacet_index_min, int cofacet_index_max, int facet_index_min, int facet_index_max);

    void parallelFacetDFSMatch();

    void parallelDirectionalFacetDFSMatch();

    void serialCofacetDFSMatch();

    int serialCycleRemoval();

    int dfsCycleRemoval();

    void addEdge(int u, int v);

    Bi_Graph_Match(int leftnum, int rightnum, int leftdim, int threadnum);

    void updateDimension(int newleftnum, int newrightnum);

    void buildInterface(const std::vector<std::vector<int>>& cofacet_bin, int cofacet_index_min, int cofacet_index_max, const std::vector<std::vector<int>>& simplex_bin, int simplex_index_max, const std::vector<int>& active_index);

    std::vector<int> getActiveIndex(); 
    
    std::vector<int> getCriticalIndex(const std::vector<int>& dim_active_index, int simplex_index_max);

    std::vector<std::set<int>> getBackwardSingleFacetIndex(std::vector<int> critical_facet_index);

    //move to private later
    std::vector<std::vector<int>> adj_list;
    std::vector<int> match_list;

    void checkSimplex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_simplex);

    void checkCofacet(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_cofacet);

    void checkSimplexByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_simplex_index);

    void checkCofacetByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_cofacet_index);

};