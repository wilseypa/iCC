#pragma once

#include <vector>
#include <deque>

class Bi_Graph_Match {
private:
    int u, v;
    
    int udegree;
    int maxthreadnum;    //affect the memory consumption during bipartite matching

    std::vector<std::vector<int>> adj_list;
    std::vector<int> match_list;
    std::vector<int> root_flag;

    //parallel init for bipartite matching
    void parallelKarpSipserInit();
    void pairDegreeOne(int uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg);
    void pairUnmatched(int uidx, std::vector<int>& visit_flag);

    //parallel bipartite matching helper function
    int dfsAugPath(int startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid);

    //cycle removal helper function
    int findRoot();
    int getParent(std::vector<int>& parent_workspace, int uidx);
    int getChild(std::vector<int>& child_workspace, int uidx);
    void getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx);
    bool isBackwardAcyclic(std::vector<int>& ancestor_simp, std::vector<int>& u_child);
    int lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_simp, std::vector<int>& child_simp, int rootnum, int uidx);
    int serialBFS(std::vector<int>& ancestor_simp, int rootnum, int root);

public:
    void parallelDFSMatch();

    int serialCycleRemoval();

    std::vector<int> getCriticalIndex();

    void addEdge(int u, int v);

    Bi_Graph_Match(int leftnum, int rightnum, int leftdim, int threadnum);

    void updateDimension(int newleftnum, int newrightnum);

    void buildInterface(std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& cofacet_bin, std::vector<int>& active_index);

    std::vector<int> getActiveIndex();
};