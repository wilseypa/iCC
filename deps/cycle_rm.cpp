#include <memory>
#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <stack>
#include <queue>
#include <unordered_set>
#include <execution>
#include <atomic>
#include <thread>
#include <set>

#include <iostream>

#include "bi_graph.hpp"
#include "cycle_rm.hpp"

Bi_Graph_Traversal::Bi_Graph_Traversal(Bi_Graph *bi_graph) : graphptr(bi_graph)
{
    visit_flag.resize(graphptr->u, 0);
}

int Bi_Graph_Traversal::findRootIterative(int maxdegree)
{
    int u = graphptr->u;

    // check the number of matched neighbors of d-1 simp
    std::vector<int> u_range(u);
    std::iota(u_range.begin(), u_range.end(), 0);
    auto match_lambda = [&](const auto &i)
    { return (graphptr->match[i] != -1); };

    int root = -1;
    int deg = 1;
    while (deg <= maxdegree)
    {
        for (auto &i : u_range)
        {
            if (visit_flag[i] > 0)
                continue;
            // par?
            int matchnum = std::count_if(std::begin(graphptr->adj_list[i]), std::end(graphptr->adj_list[i]), match_lambda);
            if (matchnum <= deg)
            {
                root = i;
                return root;
            }
        }
        deg += 1;
    }
    return root;
}

// directly return vec. use rvo
std::vector<int> Bi_Graph_Traversal::getParent(int uidx)
{
    std::vector<int> parent_d_simp;
    int vidx = graphptr->match[uidx];
    // uidx is a source node has no parent (umatched)
    if (vidx == -1)
        return parent_d_simp;
    // i is d simp
    for (auto &i : graphptr->adj_list[vidx])
    {
        int imate = graphptr->match[i];
        if (imate != -1 && imate != vidx)
            parent_d_simp.push_back(i);
    }
    return parent_d_simp;
}

// no need to check visit flag. if found child in searched part. it might be merging or cycle
// in later stages, need to check visit flag when pushing to bfs queue
std::vector<int> Bi_Graph_Traversal::getChild(int uidx)
{
    std::vector<int> child_d_simp;
    // i is d-1 simp
    for (auto &i : graphptr->adj_list[uidx])
    {
        int imate = graphptr->match[i];
        if ((imate != -1) && (imate != uidx))
            child_d_simp.push_back(imate);
    }

    // std::cout<<"uidx = "<<uidx<<"  child simp = ";
    // for (auto& i : child_d_simp) std::cout<<i<<"  ";
    // std::cout<<'\n';

    return child_d_simp;
}

// ancestor in bfs searched part
std::vector<int> Bi_Graph_Traversal::getAncestor(int uidx)
{
    // std::fill(std::execution::par, visit_flag.begin(), visit_flag.end(), 0);

    std::vector<int> ancestor_d_simp;
    std::queue<int> bfs_queue;

    std::vector<int> uidx_parent = std::move(getParent(uidx));
    for (auto &i : uidx_parent)
    {
        if (visit_flag[i] > 0)
            bfs_queue.push(i);
    }

    while (!bfs_queue.empty())
    {
        int front = bfs_queue.front();
        ancestor_d_simp.push_back(front);
        bfs_queue.pop();

        std::vector<int> front_parent = std::move(getParent(front));
        for (auto &i : front_parent)
        {
            if (visit_flag[i] > 0)
                bfs_queue.push(i);
        }
    }

    return ancestor_d_simp;
}

bool Bi_Graph_Traversal::isBackwardAcyclic(std::vector<int> &ancestor_d_simp, std::vector<int> &u_child, int uidx)
{
    // par
    for (auto &i : u_child)
    {
        if (std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) != std::end(ancestor_d_simp))
        {
            return false;
        }
    }
    return true;
}

int Bi_Graph_Traversal::lookAheadDFS(std::queue<int> &graph_bfs_queue, std::vector<int> &ancestor_d_simp, int uidx)
{
    std::vector<int> child_d_simp = std::move(getChild(uidx));
    if (child_d_simp.size() == 0)
        return 0;

    std::stack<int> lookahead_stack;
    for (auto &i : child_d_simp)
    {
        std::vector<int> i_parent = std::move(getParent(i));
        // if uidx is the only parent of i (the whole graph). use look ahead shortcut
        if (i_parent.size() == 1)
        {
            lookahead_stack.push(i);
            visit_flag[i] = 1;
        }
        else if (visit_flag[i] < 1)
        {
            // merging node. push to global bfs queue
            graph_bfs_queue.push(i);
            visit_flag[i] = 1;
        }
    }

    // stat var
    int maxsize = 0;

    // start look ahead op on children of uidx
    ancestor_d_simp.push_back(uidx);
    int reverted = 0;
    int flag; // working var to indicate the end of lookahead of each node
    // the nodes appeared in stack at the same time are independent
    while (!lookahead_stack.empty())
    {
        // for stat
        if (lookahead_stack.size() > maxsize)
            maxsize = lookahead_stack.size();

        flag = 0;

        int top = lookahead_stack.top();
        lookahead_stack.pop();

        std::vector<int> top_child_d_simp = std::move(getChild(top));
        if (top_child_d_simp.size() == 0)
            continue;

        // if top is not acyclic. the lookahead of top is over
        // push its unsearched children to bfs queue
        if (!isBackwardAcyclic(ancestor_d_simp, top_child_d_simp, top))
        {
            int temp = graphptr->match[top];
            graphptr->match[top] = -1;
            graphptr->match[temp] = -1;
            reverted += 1;
            for (auto &i : top_child_d_simp)
            {
                if (visit_flag[i] < 1)
                {
                    graph_bfs_queue.push(i);
                    visit_flag[i] = 1;
                }
            }
        }
        else
        {
            // check if the child of top can be pushed to look ahead stack
            for (auto &i : top_child_d_simp)
            {
                std::vector<int> i_parent = std::move(getParent(i));
                // if top is the only parent of i. keep look ahead op
                if (i_parent.size() == 1)
                {
                    lookahead_stack.push(i);
                    visit_flag[i] = 1;
                    flag = 1;
                }
                else
                {
                    // if more than one parent. check visit flag and push to queue
                    if (visit_flag[i] < 1)
                    {
                        graph_bfs_queue.push(i);
                        visit_flag[i] = 1;
                    }
                }
            }
        }
        // if any of top's child is pushed to look ahead stack, update ancestor
        if (flag)
            ancestor_d_simp.push_back(top);
        // not all nodes in the stack have ances-des relationship with top
        // but this does not bring extra cycle rm. simplfied implementation
    }

    std::cout << "max LA stack size = " << maxsize << '\n';

    return reverted;
}

int Bi_Graph_Traversal::serialBFS(std::queue<int> &graph_bfs_queue, int root)
{
    int reverted = 0;

    graph_bfs_queue.push(root);
    visit_flag[root] = 1;

    while (!graph_bfs_queue.empty())
    {
        int front = graph_bfs_queue.front();
        graph_bfs_queue.pop();

        std::vector<int> ancestor_d_simp = std::move(getAncestor(front));
        std::vector<int> front_child = std::move(getChild(front));

        if (!isBackwardAcyclic(ancestor_d_simp, front_child, front))
        {
            int temp = graphptr->match[front];
            graphptr->match[front] = -1;
            graphptr->match[temp] = -1;
            reverted += 1;
            for (auto &i : front_child)
            {
                if (visit_flag[i] < 1)
                {
                    graph_bfs_queue.push(i);
                    visit_flag[i] = 1;
                }
            }
        }
        else
        {
            // look ahead and push descendant to bfs queue
            reverted += lookAheadDFS(graph_bfs_queue, ancestor_d_simp, front);
        }
    }
    return reverted;
}

int Bi_Graph_Traversal::serialCycleRemoval(int maxdegree)
{
    // assume u is the d simplex of d interface
    int u = graphptr->u;

    int reverted = 0;

    int root = findRootIterative(maxdegree);
    std::queue<int> graph_bfs_queue;

    while (root != -1)
    {

        std::cout << "current root = " << root << '\n';

        reverted += serialBFS(graph_bfs_queue, root);

        root = findRootIterative(maxdegree);
    }
    return reverted;
}

void Bi_Graph_Traversal::threadBackwardDFS(int &iteration, int uidx)
{
    std::vector<int> dfs_stack;
    std::set<int> visited_set;

    dfs_stack.push_back(uidx);
    visited_set.insert(uidx);

    bool flag = true;

    while (!dfs_stack.empty())
    {
        int top = dfs_stack.back();
        dfs_stack.erase(--dfs_stack.end());
        std::vector<int> top_child = getChild(top);

        for (auto &i : top_child)
        {
            // only do dfs in searched part of current iteration("bfs componnet")
            if (root_flag[i] == iteration)
            {
                // check if i is seen before
                // maybe set::contains (log complexity)
                if (visited_set.find(i) == visited_set.end())
                {
                    // found cycle
                    int temp = graphptr->match[uidx];
                    graphptr->match[uidx] = -1;
                    graphptr->match[temp] = -1;
                    return;
                }
                dfs_stack.push_back(i);
                visited_set.insert(i);
                flag &= false;
            }
        }

        if (flag)
            visited_set.erase(top);
    }
    return;
}

// Rathod Alg 1 procedure BFSComponent
int Bi_Graph_Traversal::parallelBFS(int iteration, int root)
{
    int reverted = 0;
    root_flag[root] = iteration;

    std::queue<int> graph_bfs_queue;
    graph_bfs_queue.push(root);

    while (!graph_bfs_queue.empty())
    {
        int front = graph_bfs_queue.front();
        graph_bfs_queue.pop();

        // get "leading up-edges"
        std::vector<int> front_child = getChild(front);
        std::remove_if(front_child.begin(), front_child.end(), [&](const auto &i)
                       { return visit_flag[i] > 0; });

#pragma omp parallel for schedule(dynamic) reduction(+ : reverted)
        for (int i = 0; i < front_child.size(); i++)
        {
            threadBackwardDFS(iteration, front_child[i]);
            reverted++;
        }
        // back to the main thread, mark up and push acyclic leading up-edges to the queue
        std::copy_if(front_child.begin(), front_child.end(), std::back_inserter(graph_bfs_queue), [&](const auto &i)
                     { return graphptr->match[i] != -1; });
        for (auto &i : front_child)
        {
            root_flag[i] = iteration;
        }
    }

    return reverted;
}

int Bi_Graph_Traversal::parallelRathod(int maxdegree, int threadnum)
{
    // assume u is the d simplex of d interface
    int u = graphptr->u;

    int reverted = 0;
    int iteration = 0;

    int root = findRootIterative(maxdegree);
    while (root != -1)
    {
        iteration += 1;
        reverted += parallelBFS(iteration, root);
        root = findRootIterative(maxdegree);
    }
    return reverted;
}