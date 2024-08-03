#include <iostream>
#include <random>
#include <set>

#include "bi_graph.hpp"

Bi_Graph::Bi_Graph(int leftnum, int rightnum, int maxdeg, bool random) : u(leftnum), v(rightnum), maxdegree(maxdeg) {
    if (maxdegree > v) {
        maxdegree = v;
    }
    adj_list.resize(u + v);
    match.resize(u + v, -1);
    if (random) {
        randomBiGraphGen();
    }
    else {
        std::cout << "only random graph at this time" << std::endl;
    }
}

Bi_Graph::Bi_Graph(int leftnum, int rightnum, int maxdeg = 0) : u(leftnum), v(rightnum), maxdegree(maxdeg)
{
    if (maxdegree == 0 || maxdegree > v) {
        maxdegree = v;
    }
    adj_list.resize(u + v);
    match.resize(u + v, -1);
}

void Bi_Graph::addEdge(int u, int v)
{
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
}

void Bi_Graph::addEdge(int u, std::vector<int>& v_list)
{
    adj_list[u] = v_list;
    for (auto v : v_list)
        adj_list[v].push_back(u);
}

void Bi_Graph::printBiGraph() {
    if (adj_list.size() > MAX_PRINT_SIZE) {
        std::cout << "graph too large" << std::endl;
        return;
    }
    for (int i = 0; i < u; i++) {
        std::cout << "left node " << i << "->";
        for (auto j : adj_list[i]) {
            std::cout << j << " ";
        }
        std::cout << '\n';
    }
    for (int i = u; i < (u + v); i++) {
        std::cout << "right node " << i << "->";
        for (auto j : adj_list[i]) {
            std::cout << j << " ";
        }
        std::cout << '\n';
    }
}

void Bi_Graph::randomBiGraphGen() {
    std::random_device rand_dev;
    std::mt19937 rand_int_gen(rand_dev());
    std::set<int> u_pair;
    for (int i = 0; i < u; i++) {
        std::uniform_int_distribution<> deg_dis(1, maxdegree);
        std::uniform_int_distribution<> pair_dis(0, v - 1);
        u_pair.clear();
        int udeg = deg_dis(rand_int_gen);
        for (int j = 0; j < udeg; j++) {
            int vidx = pair_dis(rand_int_gen) + u;    //global index of nodes in v
            if (u_pair.insert(vidx).second) {
                adj_list[i].push_back(vidx);
                adj_list[vidx].push_back(i);
            }
        }
    }
}

// int main() {
//     Bi_Graph test_rand_graph = Bi_Graph(7, 7, 3, true);
//     test_rand_graph.printBiGraph();
//     return 0;
// }