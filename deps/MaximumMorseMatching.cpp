#include <execution>
#include <iostream>
#include <queue>

#include "omp.h"

#include "MaximumMorseMatching.hpp"

size_t MaximumMorseMatching::match(MatchingContext &matching_context)
{
    // homology
    // parallelMaxFacetInit(matching_context);
    // return serialCofacetMatch(matching_context.graph);

    // cohomology
    parallelMinCofacetInit(matching_context);

    return serialFacetMatch(matching_context.graph);
}

size_t MaximumMorseMatching::matchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    // homology
    // parallelMaxFacetInit(matching_context);
    // return serialCofacetMatchWithPersistence(matching_context, dim_persistent_pair);

    // cohomology
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistence(matching_context, dim_persistent_pair);
}

size_t MaximumMorseMatching::matchWithPersistenceBackup(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    // homology
    parallelMaxFacetInit(matching_context);
    return serialCofacetMatchWithPersistence(matching_context, dim_persistent_pair);

    // cohomology
    // parallelMinCofacetInit(matching_context);
    // return serialFacetMatchWithPersistence(matching_context, dim_persistent_pair);
}

int64_t MaximumMorseMatching::matchWithPersistenceReturnMinCriticalIndex(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    // cohomology
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistenceReturnMinCriticalIndex(matching_context, dim_persistent_pair);
}

std::vector<std::vector<size_t>> MaximumMorseMatching::matchWithPersistenceReturnAugPath(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    // get gradient path for quotient from augmenting path
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistenceReturnAugPath(matching_context, dim_persistent_pair);
}

void MaximumMorseMatching::parallelMaxFacetInit(MatchingContext &matching_context)
{
    auto graph = matching_context.graph;
    // convert index
    auto u = matching_context.graph.unodes;
    auto v = matching_context.graph.vnodes;

    omp_set_num_threads(threadnum_);

// apparent pair
#pragma omp parallel for schedule(dynamic)
    for (size_t i = u; i < u + v; i++)
    {
        for (const auto &uidx : graph.adj_list[i]) // uidx in v adj is in ascending order
        {
            auto vit = graph.adj_list[uidx].begin(); // vidx in u adj is in descending order

            double cofacetweight = matching_context.sorted_cofacets[uidx].second;
            double facetweight = matching_context.sorted_facets[*vit - u].second;

            if (i == *vit) // i == max facet of uidx
            {
                graph.match_list[i] = uidx;
                graph.match_list[uidx] = i;
                break;
            }
        }
    }

    return;
}

void MaximumMorseMatching::parallelMinCofacetInit(MatchingContext &matching_context)
{
    auto graph = matching_context.graph;
    // convert index
    auto u = matching_context.graph.unodes;
    auto v = matching_context.graph.vnodes;

    omp_set_num_threads(threadnum_);

    // apparent pair, need to check the weight/diameter here?
#pragma omp parallel for schedule(dynamic)
    for (size_t i = u; i < u + v; i++)
    {
        if (graph.adj_list[i].empty())
            continue; // skip empty

        auto uidx = *(graph.adj_list[i].begin()); // uidx in v adj is in ascending order

        auto maxfacet = *(graph.adj_list[uidx].begin()); // vidx in u adj is in descending order

        double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        double facetweight = matching_context.sorted_facets[maxfacet - u].second;

        // if (i == maxfacet && cofacetweight == facetweight)    //i == min cofacet of uidx, and same weight/diameter
        if (i == maxfacet)
        {
            graph.match_list[i] = uidx;
            graph.match_list[uidx] = i;
        }
    }

    return;
}

int64_t MaximumMorseMatching::serialCofacetAugPath(BipartiteGraph &graph, const size_t cofacetindex, std::vector<size_t> &aug_path, std::vector<size_t> &facet_stack)
{
    facet_stack.clear();

    auto cmp_lamdab = [](const size_t lhs, const size_t rhs)
    { return lhs < rhs; };

    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp_lamdab)> facet_queue(cmp_lamdab);

    int64_t topindex = -1;

    int64_t unmatchedfacet = -1;

    for (auto &vidx : graph.adj_list[cofacetindex])
    {
        // if (graph.match_list[vidx] != cofacetindex) facet_queue.push(vidx);
        facet_queue.push(vidx);
    }

    // facet queue may have duplicates. skip duplicates
    while (!facet_queue.empty())
    {
        size_t topfacet = facet_queue.top();
        facet_queue.pop();

        // std::cout<<"facet queue top = "<<facet<<" facet match = "<<match_list[facet]<<" start cofacet = "<<cofacetindex<<" q size = "<<facet_queue.size()<<'\n';

        if (facet_queue.empty() || topfacet != facet_queue.top())
        {
            if (graph.match_list[topfacet] < 0)
            {
                // std::cout<<"facet = "<<facet<<"  start cofacet = "<<cofacetindex<<'\n';
                unmatchedfacet = topfacet;
            }
            else
            {
                // std::cout<<"top facet = "<<facet<<"  start cofacet = "<<cofacetindex<<'\n';
                size_t nextcofacet = graph.match_list[topfacet];
                // if (facetmate > cofacetindex) return -1;

                for (auto vidx : graph.adj_list[nextcofacet])
                {
                    // if (graph.match_list[vidx] != nextcofacet) facet_queue.push(vidx);
                    if (vidx != topfacet)
                        facet_queue.push(vidx); // push all facets except the top facet
                }
                facet_stack.push_back(topfacet);
            }
        }
        else
            facet_queue.pop(); // pop/skip the duplicates

        if (unmatchedfacet > 0)
            break;
    }

    // find augmenting path from unmatchedfacet to cofacetindex
    if (unmatchedfacet > 0)
    {
        // std::cout<<"max faccet rel idx = "<<maxfacet - u<<"  ending cofacet idx = "<<cofacetindex<<"  first cofacet of max facet = "<<adj_list[maxfacet][0]<<'\n';

        aug_path[++topindex] = unmatchedfacet;
        size_t topfacet = aug_path[topindex];
        for (auto vit = facet_stack.rbegin(); vit != facet_stack.rend(); ++vit)
        {
            size_t matchedcofacet = graph.match_list[*vit];
            for (auto facet : graph.adj_list[matchedcofacet])
            {
                if (facet == topfacet)
                {
                    aug_path[++topindex] = matchedcofacet;
                    aug_path[++topindex] = *vit;
                    topfacet = *vit;
                    break;
                }
            }
        }

        aug_path[++topindex] = cofacetindex;

        return topindex + 1;
    }

    return -1;
}

int64_t MaximumMorseMatching::serialFacetAugPath(BipartiteGraph &graph, const size_t facetindex, std::vector<size_t> &aug_path, std::vector<size_t> &cofacet_stack)
{
    cofacet_stack.clear();

    auto cmp_lamdab = [](const size_t lhs, const size_t rhs)
    { return lhs > rhs; };

    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp_lamdab)> cofacet_queue(cmp_lamdab);

    int64_t topindex = -1;

    int64_t unmatchedcofacet = -1;

    for (auto &uidx : graph.adj_list[facetindex])
    {
        // if (match_list[uidx] != facetindex) cofacet_queue.push(uidx);
        cofacet_queue.push(uidx);
    }

    while (!cofacet_queue.empty())
    {
        size_t topcofacet = cofacet_queue.top();
        cofacet_queue.pop();

        if (cofacet_queue.empty() || topcofacet != cofacet_queue.top())
        {
            if (graph.match_list[topcofacet] < 0)
            {
                unmatchedcofacet = topcofacet;
            }
            else
            {
                size_t nextfacet = graph.match_list[topcofacet];
                for (auto uidx : graph.adj_list[nextfacet])
                {
                    // if (match_list[ui] != cofacetmate) cofacet_queue.push(ui);
                    if (uidx != topcofacet)
                        cofacet_queue.push(uidx); // push all cofacets except the topcofacet
                }
                cofacet_stack.push_back(topcofacet);
            }
        }
        else
            cofacet_queue.pop(); // pop/skip the duplicates

        if (unmatchedcofacet >= 0)
            break;
    }

    // find augmenting path from mincofacet to facetindex
    if (unmatchedcofacet >= 0)
    {
        aug_path[++topindex] = unmatchedcofacet;
        size_t topcofacet = aug_path[topindex];
        for (auto uit = cofacet_stack.rbegin(); uit != cofacet_stack.rend(); ++uit)
        {
            size_t matchedfacet = graph.match_list[*uit];
            for (auto uidx : graph.adj_list[matchedfacet])
            // if (std::find(adj_list[vidx].begin(), adj_list[vidx].end(), topcofacet) != adj_list[vidx].end())
            {
                if (uidx == topcofacet)
                {
                    aug_path[++topindex] = matchedfacet;
                    aug_path[++topindex] = *uit;
                    topcofacet = *uit;
                    break;
                }
            }
        }
        aug_path[++topindex] = facetindex;
        return topindex + 1;
    }

    return -1;
}

size_t MaximumMorseMatching::serialCofacetMatch(BipartiteGraph &graph)
{
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> facet_stack;
    facet_stack.reserve(graph.vnodes);

    size_t count = 0;

    for (size_t i = 0; i < graph.unodes; i++)
    {
        if (graph.match_list[i] >= 0)
            continue; // skip matched

        int64_t augpathlen = serialCofacetAugPath(graph, i, aug_path, facet_stack);

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }
    }

    for (size_t i = graph.unodes; i < graph.unodes + graph.vnodes; i++)
    {
        if (graph.match_list[i] < 0 && graph.adj_list[i].size() > 0)
            count += 1; // count unmatched facets
    }

    return count;
}

size_t MaximumMorseMatching::serialFacetMatch(BipartiteGraph &graph)
{
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    size_t count = 0;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0)
            continue; // skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        if (augpathlen < 0)
        {
            count += 1;
            continue;
        }

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }
    }

    return count;
}

size_t MaximumMorseMatching::serialCofacetMatchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    auto &graph = matching_context.graph;

    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> facet_stack;
    facet_stack.reserve(graph.vnodes);

    size_t count = 0;

    for (size_t i = 0; i < graph.unodes; i++)
    {
        if (graph.match_list[i] >= 0)
            continue; // skip matched

        int64_t augpathlen = serialCofacetAugPath(graph, i, aug_path, facet_stack);

        double facetweight = -1;
        double cofacetweight = -1;

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            if (j == 0)
            {
                auto facetidx = aug_path[0];
                facetweight = matching_context.sorted_facets[facetidx - graph.unodes].second;
                auto cofacetidx = aug_path[augpathlen - 1];
                cofacetweight = matching_context.sorted_cofacets[cofacetidx].second;
            }
            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }

        // the essential feature (live feature) is not pushed to the persistent pair
        if (facetweight != cofacetweight)
        {
            dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout << "facetweight = " << facetweight << "  cofacetweight = " << cofacetweight << '\n';
        }
    }

    for (size_t i = graph.unodes; i < graph.unodes + graph.vnodes; i++)
    {
        if (graph.match_list[i] < 0 && graph.adj_list[i].size() > 0)
            count += 1; // count unmatched facets
    }

    return count;
}

size_t MaximumMorseMatching::serialFacetMatchWithPersistence(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    auto &graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    size_t count = 0;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0)
            continue; // skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        if (augpathlen < 0)
        {
            count += 1;
            continue;
        }

        double facetweight = -1;
        double cofacetweight = -1;

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            if (j == 0)
            {
                auto cofacetidx = aug_path[0];
                cofacetweight = matching_context.sorted_cofacets[cofacetidx].second;
                auto facetidx = aug_path[augpathlen - 1];
                facetweight = matching_context.sorted_facets[facetidx - graph.unodes].second;
            }

            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }

        if (facetweight != cofacetweight)
        {
            dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout << "facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }
    }

    return count;
}

int64_t MaximumMorseMatching::serialFacetMatchWithPersistenceReturnMinCriticalIndex(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    auto &graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    // std::cout<<"last dim matching function is called. cofacets num = "<<graph.unodes<<'\n';

    int64_t minindex = -1;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0)
            continue; // skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        if (augpathlen < 0)
        {
            minindex = vidx;
            continue;
        }

        double facetweight = -1;
        double cofacetweight = -1;

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            if (j == 0)
            {
                auto cofacetidx = aug_path[0];
                cofacetweight = matching_context.sorted_cofacets[cofacetidx].second;
                auto facetidx = aug_path[augpathlen - 1];
                facetweight = matching_context.sorted_facets[facetidx - graph.unodes].second;
            }

            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }

        if (facetweight != cofacetweight)
        {
            dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout << "maxdim pair:  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }
    }

    return minindex;
}

std::vector<std::vector<size_t>> MaximumMorseMatching::serialFacetMatchWithPersistenceReturnAugPath(MatchingContext &matching_context, std::vector<std::pair<double, double>> &dim_persistent_pair)
{
    auto &graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    // std::cout<<"last dim matching function is called. cofacets num = "<<graph.unodes<<'\n';

    std::vector<std::vector<size_t>> aug_path_cofacet_vec;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0)
            continue; // skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        if (augpathlen < 0)
            continue;

        double facetweight = -1;
        double cofacetweight = -1;

        std::vector<size_t> aug_path_cofacets;

        for (int64_t j = 0; j < augpathlen; j += 2)
        {
            if (j == 0)
            {
                auto cofacetidx = aug_path[0];
                cofacetweight = matching_context.sorted_cofacets[cofacetidx].second;
                auto facetidx = aug_path[augpathlen - 1];
                facetweight = matching_context.sorted_facets[facetidx - graph.unodes].second;
            }

            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];

            aug_path_cofacets.push_back(aug_path[j]);
        }

        if (facetweight != cofacetweight)
        {
            aug_path_cofacet_vec.push_back(std::move(aug_path_cofacets));
            dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout << "maxdim pair:  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }
    }

    return aug_path_cofacet_vec;
}