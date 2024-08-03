#include "hopcroft_karp.hpp"

#include <queue>
#include <climits>

HKGraph::HKGraph() = default;

HKGraph::HKGraph(int m, int n)
{
    this->m = m;
    this->n = n;
    adj = std::vector<std::vector<int>>(m + 1);
}

void HKGraph::addEdge(int u, int v)
{
    adj[u].push_back(v);
}

std::vector<std::pair<int, int>> HKGraph::hopcroftKarpAlgorithm()
{
    pair_u = std::vector<int>(m + 1, 0);
    pair_v = std::vector<int>(n + 1, 0);

    dist = std::vector<int>(m + 1);
    std::vector<std::pair<int, int>> matchingEdges;
    while (bfs())
    {
        for (int u = 1; u <= m; u++)
        {
            if (pair_u[u] == 0 && dfs(u))
                matchingEdges.emplace_back(u, pair_u[u]);
        }
    }
    return matchingEdges;
}

std::pair<std::set<int>, std::set<int>> HKGraph::custom_hopcroftKarpAlgorithm()
{
    pair_u = std::vector<int>(m + 1, 0);
    pair_v = std::vector<int>(n + 1, 0);

    dist = std::vector<int>(m + 1);
    std::pair<std::set<int>, std::set<int>> matchingEdges;
    while (bfs())
    {
        for (int u = 1; u <= m; u++)
        {
            if (pair_u[u] == 0 && dfs(u))
            {
                matchingEdges.first.insert(u);
                matchingEdges.second.insert(pair_u[u]);
            }
        }
    }
    return matchingEdges;
}

bool HKGraph::bfs()
{
    std::queue<int> q;
    for (int u = 1; u <= m; u++)
    {
        if (pair_u[u] == 0)
        {
            dist[u] = 0;
            q.push(u);
        }
        else
        {
            dist[u] = INT_MAX;
        }
    }
    dist[0] = INT_MAX;
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        if (dist[u] < dist[0])
        {
            for (int v : adj[u])
            {
                if (dist[pair_v[v]] == INT_MAX)
                {
                    dist[pair_v[v]] = dist[u] + 1;
                    q.push(pair_v[v]);
                }
            }
        }
    }
    return (dist[0] != INT_MAX);
}

bool HKGraph::dfs(int u)
{
    if (u != 0)
    {
        for (auto v : adj[u])
        {
            if (dist[pair_v[v]] == dist[u] + 1 && dfs(pair_v[v]))
            {
                pair_v[v] = u;
                pair_u[u] = v;
                return true;
            }
        }
        dist[u] = INT_MAX;
        return false;
    }
    return true;
}
