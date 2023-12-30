

#include "hopcroft_karp.hpp"

std::vector<std::pair<int, int>> HKGraph::hopcroftKarpAlgorithm()
{
    pair_u = std::vector<int>(m + 1, NIL);
    pair_v = std::vector<int>(n + 1, NIL);

    dist = std::vector<int>(m + 1);
    std::vector<std::pair<int, int>> matchingEdges;
    while (bfs())
    {
        for (int u = 1; u <= m; u++)
        {
            if (pair_u[u] == NIL && dfs(u))
                matchingEdges.emplace_back(u, pair_u[u]);
        }
    }
    return matchingEdges;
}

std::pair<std::set<int>, std::set<int>> HKGraph::custom_hopcroftKarpAlgorithm()
{
    pair_u = std::vector<int>(m + 1, NIL);
    pair_v = std::vector<int>(n + 1, NIL);

    dist = std::vector<int>(m + 1);
    std::pair<std::set<int>, std::set<int>> matchingEdges;
    while (bfs())
    {
        for (int u = 1; u <= m; u++)
        {
            if (pair_u[u] == NIL && dfs(u))
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
        if (pair_u[u] == NIL)
        {
            dist[u] = 0;
            q.push(u);
        }
        else
        {
            dist[u] = INF;
        }
    }
    dist[NIL] = INF;
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        if (dist[u] < dist[NIL])
        {
            for (int v : adj[u])
            {
                if (dist[pair_v[v]] == INF)
                {
                    dist[pair_v[v]] = dist[u] + 1;
                    q.push(pair_v[v]);
                }
            }
        }
    }
    return (dist[NIL] != INF);
}

bool HKGraph::dfs(int u)
{
    if (u != NIL)
    {
        for (auto v : adj[u])
        {
            if (dist[pair_v[v]] == dist[u] + 1)
            {
                if (dfs(pair_v[v]) == true)
                {
                    pair_v[v] = u;
                    pair_u[u] = v;
                    return true;
                }
            }
        }
        dist[u] = INF;
        return false;
    }
    return true;
}

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