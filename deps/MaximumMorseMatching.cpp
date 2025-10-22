#include <execution>
#include <iostream>
#include <queue>

#include "omp.h"

#include "MaximumMorseMatching.hpp"
#include "SimplexUtility.hpp"

size_t MaximumMorseMatching::implicitMatch(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    auto& binom_table = matching_context.binomial_table;
    auto& cofacet_list = matching_context.sorted_cofacets;
    auto& facet_list = matching_context.sorted_facets;
    auto& cofacet_hash = matching_context.cofacet_bindex_to_list_index;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& bi_graph = matching_context.graph;
    auto u = bi_graph.unodes;
    auto v = bi_graph.vnodes;
    auto dim = matching_context.dim;    //interface dim
    auto npts = matching_context.npts;    //number of points: original + virtual

    //init workspace
    // cofacet_indices_.clear();
    cofacet_indices_.reserve(dim < 5 ? 32 : 64);

    // facet_indices_.clear();
    facet_indices_.reserve(dim + 1);

    // aug_path_.clear();
    aug_path_.reserve(u + v);

    // cofacet_stack_.clear();
    cofacet_stack_.reserve(u);

    // vertex_workspace_.clear();
    vertex_workspace_.reserve(dim + 1);

    // pq_workspace.clear();
    pq_workspace_.reserve(u);

    size_t count = 0;    //critial/unmatched count

    int k = 0;

    //process facet in reverse order
    for (int64_t i = v - 1; i >= 0; --i)
    {
        //covert list index to graph index
        size_t facetidx = i + u;
        int64_t facetbindex = facet_list[i].first;

        //skip non active facets
        auto fit = facet_hash.find(facetbindex);
        if (fit == facet_hash.end()) continue;

        SimplexUtility::getCofacetListIndicesInPlace(binom_table, cofacet_hash, cofacet_indices_, vertex_workspace_, facetbindex, npts, dim-1);

        // if (k < 10)
        // {
        //     k += 1;
        //     std::cout<<"facet list idx = "<<i<<"   "<<"cofacet list size = "<<cofacet_indices_.size()<<'\n';
        //     size_t mincofacetidx = *std::min_element(cofacet_indices_.begin(), cofacet_indices_.end());
        //     SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices_, vertex_workspace_, cofacet_list[mincofacetidx].first, npts, dim);
        //     std::cout<<"mincofacet's facets idx =  ";
        //     for (auto t : facet_indices_) std::cout<<t + u<<"  ";
        //     std::cout<<'\n';

        // }

        if (cofacet_indices_.empty())
        {
            count += 1;
            continue;
        }

        size_t mincofacetidx = *std::min_element(cofacet_indices_.begin(), cofacet_indices_.end());
        //mincofacetidx == list index == graph index for cofacet

        if (bi_graph.match_list[mincofacetidx] < 0)
        {
            SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices_, vertex_workspace_, cofacet_list[mincofacetidx].first, npts, dim);
            
            if (!facet_indices_.empty())
            {
                size_t maxfacetlistidx = *std::max_element(facet_indices_.begin(), facet_indices_.end());

                //compare the list index
                if (i == maxfacetlistidx)
                {
                    double cofacetweight = cofacet_list[mincofacetidx].second;
                    double facetweight = facet_list[i].second;

                    if (cofacetweight == facetweight)
                    {
                        //apparent pair
                        bi_graph.match_list[facetidx] = mincofacetidx;
                        bi_graph.match_list[mincofacetidx] = facetidx;

                        k+=1;

                        continue;
                    }
                }
            }
        }

        //not an apparent pair, need to find augmenting path
        // std::cout<<"will start aug path"<<'\n';
        int64_t pathlen = implicitFacetAugPath(binom_table, bi_graph, facet_list, cofacet_hash, facetidx, npts, dim);
        // std::cout<<"aug path done, path len = "<<pathlen<<'\n';

        // int64_t pathlen = -1;

        if (pathlen > 0)
        {
            auto uidx = aug_path_[0];
            double cofacetweight = matching_context.sorted_cofacets[uidx].second;
            auto vidx = aug_path_[pathlen - 1];
            double facetweight = matching_context.sorted_facets[vidx - u].second;

            for (int64_t j = 0; j < pathlen; j += 2) 
            {
                bi_graph.match_list[aug_path_[j]] = aug_path_[j + 1];
                bi_graph.match_list[aug_path_[j + 1]] = aug_path_[j];
            }

            if (facetweight != cofacetweight)
            {
                dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
                std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
            }
        }
        else 
        {
            count += 1;
        }
    }

    std::cout<<"implicit apparent pair count = "<<k<<'\n';

    return count;

}

int64_t MaximumMorseMatching::implicitFacetAugPath(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& facet_list,
                                                   const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, size_t npts, size_t interfacedimension)
{
    aug_path_.clear();
    cofacet_stack_.clear();
    pq_workspace_.clear();

    QueueWorkspaceHelper cofacet_queue(std::greater<size_t>{}, std::move(pq_workspace_));

    int64_t unmatchedcofacet = -1;

    //reuse the immediate cofacets computed for apparent pair check
    for (auto cofidx : cofacet_indices_)
    {
        cofacet_queue.push(cofidx);
    }

    while (!cofacet_queue.empty())
    {
        size_t topcofacet = cofacet_queue.top();
        cofacet_queue.pop();

        // std::cout<<"current top cofacet of queue = "<<topcofacet<<'\n';

        //skip duplicates
        bool duplicate = false;
        while (!cofacet_queue.empty() && cofacet_queue.top() == topcofacet)
        {
            cofacet_queue.pop();
            duplicate = true;
        }

        if (duplicate) continue;

        if (bi_graph.match_list[topcofacet] < 0)
        {
            unmatchedcofacet = topcofacet;
            break;
        }

        size_t nextfacet = bi_graph.match_list[topcofacet];

        //get cofacets of the nextfacet
        int64_t facetbindex = facet_list[nextfacet - bi_graph.unodes].first;
        SimplexUtility::getCofacetListIndicesInPlace(binomial_table, cofacet_hash_table, cofacet_indices_, vertex_workspace_, facetbindex, npts, interfacedimension - 1);

        for (auto cofidx : cofacet_indices_)
        {
            if (cofidx != topcofacet) cofacet_queue.push(cofidx);
        }

        cofacet_stack_.push_back(topcofacet);
    }

    //no aug path found
    if (unmatchedcofacet < 0) return -1;

    //reconstruct the actual path
    // int64_t topindex = -1;

    // aug_path_[++topindex] = unmatchedcofacet;
    aug_path_.push_back(unmatchedcofacet);

    size_t currenttopcofacet = aug_path_.back();

    //backtrack through the cofacet stack
    for (auto it = cofacet_stack_.rbegin(); it != cofacet_stack_.rend(); ++it)
    {
        size_t matchedfacet = bi_graph.match_list[*it];

        //check if the matchedfacet is adjacent to the current top cofacet on the path
        int64_t facetbindex = facet_list[matchedfacet - bi_graph.unodes].first;
        SimplexUtility::getCofacetListIndicesInPlace(binomial_table, cofacet_hash_table, cofacet_indices_, vertex_workspace_, facetbindex, npts, interfacedimension - 1);

        if (std::find(cofacet_indices_.begin(), cofacet_indices_.end(), currenttopcofacet) != cofacet_indices_.end())
        {
            aug_path_.push_back(matchedfacet);
            aug_path_.push_back(*it);
            currenttopcofacet = *it;
        }
    }

    aug_path_.push_back(facetgraphindex);    //initial facet index

    //move the container back
    pq_workspace_ = std::move(cofacet_queue.getContainer());

    return aug_path_.size();
}


size_t MaximumMorseMatching::match(MatchingContext& matching_context)
{
    // homology
    // parallelMaxFacetInit(matching_context);
    // return serialCofacetMatch(matching_context.graph);

    // cohomology
    parallelMinCofacetInit(matching_context);
    
    return serialFacetMatch(matching_context.graph);
}

size_t MaximumMorseMatching::matchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    // homology
    // parallelMaxFacetInit(matching_context);
    // return serialCofacetMatchWithPersistence(matching_context, dim_persistent_pair);

    // cohomology
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistence(matching_context, dim_persistent_pair);
}

size_t MaximumMorseMatching::matchWithPersistenceBackup(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    // homology
    parallelMaxFacetInit(matching_context);
    return serialCofacetMatchWithPersistence(matching_context, dim_persistent_pair);

    // cohomology
    // parallelMinCofacetInit(matching_context);
    // return serialFacetMatchWithPersistence(matching_context, dim_persistent_pair);
}

int64_t MaximumMorseMatching::matchWithPersistenceReturnMinCriticalIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    // cohomology
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistenceReturnMinCriticalIndex(matching_context, dim_persistent_pair);
}

std::vector< std::vector<size_t> > MaximumMorseMatching::matchWithPersistenceReturnAugPath(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    //get gradient path for quotient from augmenting path
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchWithPersistenceReturnAugPath(matching_context, dim_persistent_pair);
}


void MaximumMorseMatching::parallelMaxFacetInit(MatchingContext& matching_context)
{
    auto graph = matching_context.graph;
    //convert index
    auto u = matching_context.graph.unodes;
    auto v = matching_context.graph.vnodes;

    omp_set_num_threads(threadnum_);

//apparent pair
#pragma omp parallel for schedule(dynamic)
    for (size_t i = u; i < u + v; i++)
    {
        for (const auto& uidx: graph.adj_list[i])    //uidx in v adj is in ascending order
        {
            auto vit = graph.adj_list[uidx].begin();    //vidx in u adj is in descending order

            double cofacetweight = matching_context.sorted_cofacets[uidx].second;
            double facetweight = matching_context.sorted_facets[*vit - u].second;

            if (i == *vit)    //i == max facet of uidx
            {
                graph.match_list[i] = uidx;
                graph.match_list[uidx] = i;
                break;
            }
        }        
    }

    return;

}

void MaximumMorseMatching::parallelMinCofacetInit(MatchingContext& matching_context)
{
    auto graph = matching_context.graph;
    //convert index
    auto u = matching_context.graph.unodes;
    auto v = matching_context.graph.vnodes;

    omp_set_num_threads(threadnum_);

    int k = 0;
    
    //apparent pair, need to check the weight/diameter here?
#pragma omp parallel for schedule(dynamic)
    for (size_t i = u; i < u + v; i++)
    {
        if (graph.adj_list[i].empty()) continue;    //skip empty

        auto uidx = *(graph.adj_list[i].begin());    //uidx in v adj is in ascending order

        auto maxfacet = *(graph.adj_list[uidx].begin());    //vidx in u adj is in descending order

        double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        double facetweight = matching_context.sorted_facets[maxfacet - u].second;

        // if (i == maxfacet && cofacetweight == facetweight)    //i == min cofacet of uidx, and same weight/diameter
        if (i == maxfacet)
        {
            k+=1;

            graph.match_list[i] = uidx;
            graph.match_list[uidx] = i;
        }
    }

    std::cout<<"explicit apparent pair count = "<<k<<"\n";

    return;
}


int64_t MaximumMorseMatching::serialCofacetAugPath(BipartiteGraph& graph, const size_t cofacetindex, std::vector<size_t>& aug_path, std::vector<size_t>& facet_stack)
{
    facet_stack.clear();

    auto cmp_lamdab = [] (const size_t lhs, const size_t rhs) { return lhs < rhs; };
    
    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp_lamdab)> facet_queue(cmp_lamdab);

    int64_t topindex = -1;

    int64_t unmatchedfacet = -1;


    for (auto& vidx: graph.adj_list[cofacetindex]) 
    {
        // if (graph.match_list[vidx] != cofacetindex) facet_queue.push(vidx);
        facet_queue.push(vidx);
    }

    //facet queue may have duplicates. skip duplicates
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

                for (auto vidx: graph.adj_list[nextcofacet])
                {
                    // if (graph.match_list[vidx] != nextcofacet) facet_queue.push(vidx);
                    if (vidx != topfacet) facet_queue.push(vidx);    //push all facets except the top facet
                }
                facet_stack.push_back(topfacet);
            }
        }
        else facet_queue.pop();     //pop/skip the duplicates

        if (unmatchedfacet > 0) break;
    }


    //find augmenting path from unmatchedfacet to cofacetindex
    if (unmatchedfacet > 0)
    {
        // std::cout<<"max faccet rel idx = "<<maxfacet - u<<"  ending cofacet idx = "<<cofacetindex<<"  first cofacet of max facet = "<<adj_list[maxfacet][0]<<'\n';

        aug_path[++topindex] = unmatchedfacet;
        size_t topfacet = aug_path[topindex];
        for (auto vit = facet_stack.rbegin(); vit != facet_stack.rend(); ++vit)
        {
            size_t matchedcofacet = graph.match_list[*vit];
            for (auto facet: graph.adj_list[matchedcofacet])
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


int64_t MaximumMorseMatching::serialFacetAugPath(BipartiteGraph& graph, const size_t facetindex, std::vector<size_t>& aug_path, std::vector<size_t>& cofacet_stack)
{
    cofacet_stack.clear();

    auto cmp_lamdab = [] (const size_t lhs, const size_t rhs) { return lhs > rhs; };

    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp_lamdab)> cofacet_queue(cmp_lamdab);

    int64_t topindex = -1;

    int64_t unmatchedcofacet = -1;

    for (auto& uidx : graph.adj_list[facetindex])
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
                for (auto uidx: graph.adj_list[nextfacet])
                {
                    // if (match_list[ui] != cofacetmate) cofacet_queue.push(ui);
                    if (uidx != topcofacet) cofacet_queue.push(uidx);    //push all cofacets except the topcofacet
                }
                cofacet_stack.push_back(topcofacet);
            }
        }
        else cofacet_queue.pop();     //pop/skip the duplicates

        if (unmatchedcofacet >= 0) break;
    }

    //find augmenting path from mincofacet to facetindex
    if (unmatchedcofacet >= 0)
    {
        aug_path[++topindex] = unmatchedcofacet;
        size_t topcofacet = aug_path[topindex];
        for (auto uit = cofacet_stack.rbegin(); uit != cofacet_stack.rend(); ++uit)
        {
            size_t matchedfacet = graph.match_list[*uit];
            for (auto uidx: graph.adj_list[matchedfacet])
            //if (std::find(adj_list[vidx].begin(), adj_list[vidx].end(), topcofacet) != adj_list[vidx].end())
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


int64_t MaximumMorseMatching::serialFacetAugPathTest(BipartiteGraph& graph, const size_t facetindex, std::vector<size_t>& aug_path, std::vector<size_t>& cofacet_stack)
{
    cofacet_stack.clear();

    auto cmp_lamdab = [] (const size_t lhs, const size_t rhs) { return lhs > rhs; };

    std::priority_queue<size_t, std::vector<size_t>, decltype(cmp_lamdab)> cofacet_queue(cmp_lamdab);

    int64_t topindex = -1;

    int64_t unmatchedcofacet = -1;

    for (auto& uidx : graph.adj_list[facetindex])
    {
        // if (match_list[uidx] != facetindex) cofacet_queue.push(uidx);
        cofacet_queue.push(uidx);
    }

    while (!cofacet_queue.empty())
    {
        size_t topcofacet = cofacet_queue.top();
        cofacet_queue.pop();

        //skip duplicates
        bool duplicate = false;
        while (!cofacet_queue.empty() && cofacet_queue.top() == topcofacet)
        {
            cofacet_queue.pop();
            duplicate = true;
        }

        if (duplicate) continue;

        if (graph.match_list[topcofacet] < 0)
        {
            unmatchedcofacet = topcofacet;
            break;
        }

        size_t nextfacet = graph.match_list[topcofacet];
        for (auto uidx: graph.adj_list[nextfacet])
        {
            // if (match_list[ui] != cofacetmate) cofacet_queue.push(ui);
            if (uidx != topcofacet) cofacet_queue.push(uidx);    //push all cofacets except the topcofacet
        }
        cofacet_stack.push_back(topcofacet);
    }

    //find augmenting path from mincofacet to facetindex
    if (unmatchedcofacet >= 0)
    {
        aug_path[++topindex] = unmatchedcofacet;
        size_t topcofacet = aug_path[topindex];
        for (auto uit = cofacet_stack.rbegin(); uit != cofacet_stack.rend(); ++uit)
        {
            size_t matchedfacet = graph.match_list[*uit];
            for (auto uidx: graph.adj_list[matchedfacet])
            //if (std::find(adj_list[vidx].begin(), adj_list[vidx].end(), topcofacet) != adj_list[vidx].end())
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


size_t MaximumMorseMatching::serialCofacetMatch(BipartiteGraph& graph)
{
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> facet_stack;
    facet_stack.reserve(graph.vnodes);

    size_t count = 0;

    for (size_t i = 0; i < graph.unodes; i++)
    {
        if (graph.match_list[i] >= 0) continue;    //skip matched 
        
        int64_t augpathlen = serialCofacetAugPath(graph, i, aug_path, facet_stack);

        for (int64_t j = 0; j < augpathlen; j += 2) 
        {           
            graph.match_list[aug_path[j]] = aug_path[j + 1];
            graph.match_list[aug_path[j + 1]] = aug_path[j];
        }
    }

    for (size_t i = graph.unodes; i < graph.unodes + graph.vnodes; i++)
    {
        if (graph.match_list[i] < 0 && graph.adj_list[i].size() > 0) count += 1;    //count unmatched facets
    }

    return count;
}


size_t MaximumMorseMatching::serialFacetMatch(BipartiteGraph& graph)
{
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    size_t count = 0;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0) continue;    //skip matched or not active facet

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


size_t MaximumMorseMatching::serialCofacetMatchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    auto& graph = matching_context.graph;

    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> facet_stack;
    facet_stack.reserve(graph.vnodes);

    size_t count = 0;

    for (size_t i = 0; i < graph.unodes; i++)
    {
        if (graph.match_list[i] >= 0) continue;    //skip matched 
        
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

        //the essential feature (live feature) is not pushed to the persistent pair
        if (facetweight != cofacetweight)
        {
            dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout<<"facetweight = "<<facetweight<<"  cofacetweight = "<<cofacetweight<<'\n';
        }
    }

    for (size_t i = graph.unodes; i < graph.unodes + graph.vnodes; i++)
    {
        if (graph.match_list[i] < 0 && graph.adj_list[i].size() > 0) count += 1;    //count unmatched facets
    }

    return count;
}


size_t MaximumMorseMatching::serialFacetMatchWithPersistence(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    auto& graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    size_t count = 0;

    int k = 0;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        // if (k < 10)
        // {
        //     k += 1;
        //     std::cout<<"facet list idx = "<<i<<"   "<<"cofacet list size = "<<graph.adj_list[vidx].size()<<'\n';
        //     auto mincofacet = graph.adj_list[vidx][0];
        //     std::cout<<"mincofacet's facets idx = ";
        //     for (auto t : graph.adj_list[mincofacet]) std::cout<<t<<"  ";
        //     std::cout<<'\n';
        // }

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0) continue;    //skip matched or not active facet

        int64_t augpathlen = serialFacetAugPathTest(graph, vidx, aug_path, cofacet_stack);

        // int64_t augpathlen = -1;

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
            std::cout<<"facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';
        }

    }

    return count;
}


int64_t MaximumMorseMatching::serialFacetMatchWithPersistenceReturnMinCriticalIndex(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    auto& graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    // std::cout<<"last dim matching function is called. cofacets num = "<<graph.unodes<<'\n';

    int64_t minindex = -1;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0) continue;    //skip matched or not active facet

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
            std::cout<<"maxdim pair:  facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';
        }

    }

    return minindex;
}

std::vector< std::vector<size_t> > MaximumMorseMatching::serialFacetMatchWithPersistenceReturnAugPath(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    auto& graph = matching_context.graph;
    std::vector<size_t> aug_path(graph.unodes + graph.vnodes, 0);

    std::vector<size_t> cofacet_stack;
    cofacet_stack.reserve(graph.unodes);

    // std::cout<<"last dim matching function is called. cofacets num = "<<graph.unodes<<'\n';

    std::vector< std::vector<size_t> > aug_path_cofacet_vec;

    for (int64_t i = graph.vnodes - 1; i >= 0; i--)
    {
        auto vidx = graph.unodes + i;

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0) continue;    //skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        if (augpathlen < 0) continue;

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
            std::cout<<"maxdim pair:  facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';
        }
    }

    return aug_path_cofacet_vec;
}