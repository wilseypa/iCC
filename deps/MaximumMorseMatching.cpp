#include <execution>
#include <iostream>

#include <cassert>

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
    cofacet_indices_.reserve(dim < 5 ? 32 : 64);

    facet_indices_.reserve(dim + 1);

    aug_path_.reserve(u + v);

    cofacet_stack_.reserve(u);

    vertex_workspace_.reserve(dim + 1);

    pq_workspace_.reserve(u);

    //counters
    size_t count = 0;    //critial/unmatched count

    size_t ct = 0;    //apparent pair count

    //process facet in reverse order
    for (int64_t i = static_cast<int64_t>(v) - 1; i >= 0; --i)
    {
        const int64_t facetbindex = facet_list[static_cast<size_t>(i)].first;

        //skip non active facets
        auto fit = facet_hash.find(facetbindex);
        if (fit == facet_hash.end()) continue;

        //covert list index to graph index
        size_t facetgraphidx = static_cast<size_t>(i) + u;

        SimplexUtility::getCofacetListIndicesInPlace(binom_table, cofacet_hash, cofacet_indices_, vertex_workspace_, facetbindex, npts, dim-1);

        if (cofacet_indices_.empty())
        {
            const double facetweight = matching_context.sorted_facets[static_cast<size_t>(i)].second;
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            count += 1;
            continue;
        }

        const size_t mincofacetidx = *std::min_element(cofacet_indices_.begin(), cofacet_indices_.end());
        //mincofacetidx == list index == graph index for cofacet

        if (bi_graph.match_list[mincofacetidx] < 0)
        {
            SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices_, vertex_workspace_, cofacet_list[mincofacetidx].first, npts, dim);
            
            if (!facet_indices_.empty())
            {
                const size_t maxfacetlistidx = *std::max_element(facet_indices_.begin(), facet_indices_.end());

                //compare the list index
                if (static_cast<size_t>(i) == maxfacetlistidx)
                {
                    const double cofacetweight = cofacet_list[mincofacetidx].second;
                    const double facetweight = facet_list[static_cast<size_t>(i)].second;

                    if (cofacetweight == facetweight)
                    {
                        //apparent pair
                        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(mincofacetidx);
                        bi_graph.match_list[mincofacetidx] = static_cast<int64_t>(facetgraphidx);

                        ct+=1;

                        continue;
                    }
                }
            }
        }

        //not an apparent pair, need to find augmenting path
        // std::cout<<"will start aug path"<<'\n';
        int64_t pathlen = implicitFacetAugPath(binom_table, bi_graph, facet_list, cofacet_hash, facetgraphidx, npts, dim);
        // std::cout<<"aug path done, path len = "<<pathlen<<'\n';

        if (pathlen <= 0)
        {
            const double facetweight = matching_context.sorted_facets[static_cast<size_t>(i)].second;
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            count += 1;
            continue;
        }

        const size_t uidx = aug_path_[0];
        const double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        const size_t vidx = aug_path_[static_cast<size_t>(pathlen) - 1];
        const double facetweight = matching_context.sorted_facets[vidx - u].second;

        for (int64_t j = 0; j < pathlen; j += 2) 
        {
            bi_graph.match_list[aug_path_[static_cast<size_t>(j)]] = static_cast<int64_t>(aug_path_[static_cast<size_t>(j) + 1]);
            bi_graph.match_list[aug_path_[static_cast<size_t>(j) + 1]] = static_cast<int64_t>(aug_path_[static_cast<size_t>(j)]);
        }

        if (facetweight != cofacetweight)
        {
            // dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }

    }

    std::cout<<"interface dim = "<<dim<<" implicit apparent pair count = "<<ct<<'\n';

    return count;
}

std::vector<std::vector<size_t>> MaximumMorseMatching::implicitMatchAndCollectPVSupports(MatchingContext& matching_context)
{
    auto& binom_table = matching_context.binomial_table;
    auto& cofacet_list = matching_context.sorted_cofacets;
    auto& facet_list = matching_context.sorted_facets;
    auto& cofacet_hash = matching_context.cofacet_bindex_to_list_index;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& bi_graph = matching_context.graph;
    auto u = bi_graph.unodes;
    auto v = bi_graph.vnodes;
    auto dim = matching_context.dim;    //interface dim (cofacet dim)
    auto npts = matching_context.npts;  //total vertex count

    // init workspace (same pattern as implicitMatch)
    cofacet_indices_.reserve(dim < 5 ? 32 : 64);
    facet_indices_.reserve(dim + 1);
    aug_path_.reserve(u + v);
    cofacet_stack_.reserve(u);
    vertex_workspace_.reserve(dim + 1);
    pq_workspace_.reserve(u);

    std::vector<std::vector<size_t>> pv_support_cofacets;
    pv_support_cofacets.reserve(dim < 5 ? 64 : 128);    //estimated numbers

    // process facets in reverse order
    for (int64_t i = static_cast<int64_t>(v) - 1; i >= 0; --i)
    {
        const int64_t facetbindex = facet_list[static_cast<size_t>(i)].first;

        // skip non-active facets
        auto fit = facet_hash.find(facetbindex);
        if (fit == facet_hash.end()) continue;

        const size_t facetgraphidx = static_cast<size_t>(i) + u;

        // compute immediate cofacets of this facet (list indices)
        SimplexUtility::getCofacetListIndicesInPlace(binom_table, cofacet_hash, cofacet_indices_, vertex_workspace_,
                                                     facetbindex, npts, dim - 1);

        //if active and have no cofacets 
        if (cofacet_indices_.empty())
        {
            const double facetweight = matching_context.sorted_facets[static_cast<size_t>(i)].second;
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            continue;
        }

        const size_t mincofacetidx = *std::min_element(cofacet_indices_.begin(), cofacet_indices_.end());

        // apparent pair check (same as implicitMatch)
        if (bi_graph.match_list[mincofacetidx] < 0)
        {
            SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices_, vertex_workspace_,
                                                       cofacet_list[mincofacetidx].first, npts, dim);

            if (!facet_indices_.empty())
            {
                const size_t maxfacetlistidx = *std::max_element(facet_indices_.begin(), facet_indices_.end());

                if (static_cast<size_t>(i) == maxfacetlistidx)
                {
                    const double cofacetweight = cofacet_list[mincofacetidx].second;
                    const double facetweight = facet_list[static_cast<size_t>(i)].second;

                    if (cofacetweight == facetweight)
                    {
                        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(mincofacetidx);
                        bi_graph.match_list[mincofacetidx] = static_cast<int64_t>(facetgraphidx);
                        continue;
                    }
                }
            }
        }

        // find augmenting path
        const int64_t pathlen = implicitFacetAugPath(binom_table, bi_graph, facet_list, cofacet_hash, facetgraphidx, npts, dim);

        if (pathlen <= 0)
        {
            const double facetweight = matching_context.sorted_facets[static_cast<size_t>(i)].second;
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            continue;
        }

        const size_t terminalcofacetidx = aug_path_[0];
        pv_support_cofacets.push_back(collectReducedColumnSupport(matching_context, terminalcofacetidx, facetgraphidx));


        const size_t uidx = aug_path_[0];
        const double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        const size_t vidx = aug_path_[static_cast<size_t>(pathlen) - 1];
        const double facetweight = matching_context.sorted_facets[vidx - u].second;

        for (int64_t j = 0; j < pathlen; j += 2) 
        {
            bi_graph.match_list[aug_path_[static_cast<size_t>(j)]] = static_cast<int64_t>(aug_path_[static_cast<size_t>(j) + 1]);
            bi_graph.match_list[aug_path_[static_cast<size_t>(j) + 1]] = static_cast<int64_t>(aug_path_[static_cast<size_t>(j)]);
        }

        if (facetweight != cofacetweight)
        {
            // dim_persistent_pair.push_back(std::make_pair(facetweight, cofacetweight));
            std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }

    }

    return pv_support_cofacets;
}


int64_t MaximumMorseMatching::implicitFacetAugPath(const std::vector<std::vector<int64_t>>& binomial_table, BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& facet_list,
                                                   const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, const size_t facetgraphindex, size_t npts, size_t interfacedimension)
{
    aug_path_.clear();
    cofacet_stack_.clear();
    pq_workspace_.clear();

    MinIndexQueue cofacet_queue(std::greater<size_t>{}, std::move(pq_workspace_));

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
    if (unmatchedcofacet < 0) 
    {
        auto& pq_buffer = cofacet_queue.getContainer();
        pq_buffer.clear();
        pq_workspace_ = std::move(pq_buffer);
        return -1;
    }

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
    auto& pq_buffer = cofacet_queue.getContainer();
    pq_buffer.clear();
    pq_workspace_ = std::move(cofacet_queue.getContainer());

    return static_cast<int64_t>(aug_path_.size());
}

std::vector<size_t> MaximumMorseMatching::collectReducedColumnSupport(const MatchingContext& matching_context, const size_t terminalcofacet, const size_t expectedfacet)
{
    auto& binom_table = matching_context.binomial_table;
    auto& cofacet_list = matching_context.sorted_cofacets;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& bi_graph = matching_context.graph;
    const size_t dim = matching_context.dim;
    const size_t npts = matching_context.npts;

    facet_indices_.clear();

    MaxIndexQueue facet_queue(std::less<size_t>{}, std::move(pq_workspace_));

    std::vector<size_t> cofacet_trace;
    cofacet_trace.reserve(dim < 5 ? 16 : 32);
    cofacet_trace.push_back(terminalcofacet);

    // initialize with boundary of terminal cofacet
    SimplexUtility::getFacetListIndicesInPlace(
        binom_table, facet_hash, facet_indices_, vertex_workspace_,
        cofacet_list[terminalcofacet].first, npts, dim);

    for (size_t fi : facet_indices_) facet_queue.push(fi);

    while (!facet_queue.empty())
    {
        size_t topfacet = facet_queue.top();    //bipartite graph index 
        facet_queue.pop();

        if (!facet_queue.empty() && topfacet == facet_queue.top())
        {
            facet_queue.pop();   // cancel duplicates
            continue;
        }

        const int64_t nextcofacet = bi_graph.match_list[topfacet];
        if (nextcofacet < 0) 
        {
            assert(expectedfacet == topfacet);    //debug purpose
            break;
        }

        const size_t cfi = static_cast<size_t>(nextcofacet);
        cofacet_trace.push_back(cfi);

        SimplexUtility::getFacetListIndicesInPlace(binom_table, facet_hash, facet_indices_, vertex_workspace_, cofacet_list[cfi].first, npts, dim);

        for (size_t fi : facet_indices_)
        {
            if (fi != topfacet) facet_queue.push(fi);
        }
    }

    auto& pq_buffer = facet_queue.getContainer();
    pq_buffer.clear();
    pq_workspace_ = std::move(pq_buffer);

    return cofacet_trace;
}



//legacy explicit matching routines

size_t MaximumMorseMatching::match(MatchingContext& matching_context)
{

    // cohomology
    parallelMinCofacetInit(matching_context);
    
    return serialFacetMatch(matching_context);
}

int64_t MaximumMorseMatching::matchAndGetMinCriticalIndex(MatchingContext& matching_context)
{
    // cohomology
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchAndGetMinCriticalIndex(matching_context);
}

std::vector< std::vector<size_t> > MaximumMorseMatching::matchAndGetAugPath(MatchingContext& matching_context, std::vector<std::pair<double, double>>& dim_persistent_pair)
{
    //get gradient path for quotient from augmenting path
    parallelMinCofacetInit(matching_context);
    return serialFacetMatchAndGetAugPath(matching_context);
}


void MaximumMorseMatching::parallelMinCofacetInit(MatchingContext& matching_context)
{
    auto graph = matching_context.graph;
    //convert index
    auto u = matching_context.graph.unodes;
    auto v = matching_context.graph.vnodes;

    omp_set_num_threads(threadnum_);
    
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
            graph.match_list[i] = uidx;
            graph.match_list[uidx] = i;
        }
    }

    return;
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


size_t MaximumMorseMatching::serialFacetMatch(MatchingContext& matching_context)
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

        if (graph.match_list[vidx] >= 0 || graph.adj_list[vidx].size() == 0) continue;    //skip matched or not active facet

        int64_t augpathlen = serialFacetAugPath(graph, vidx, aug_path, cofacet_stack);

        // int64_t augpathlen = -1;        
        


        if (augpathlen < 0)
        {
            double facetweight = matching_context.sorted_facets[i].second;
            std::cout<<"facet weight = "<<facetweight<<"  cofacet weight = -1"<<'\n';
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
            std::cout<<"facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';

            if (facetweight < 1.1) continue;

            auto facetidx = aug_path[augpathlen - 1];
            auto fbindex = matching_context.sorted_facets[facetidx - graph.unodes].first;
            auto cofacetidx = aug_path[0];
            auto cfbindex = matching_context.sorted_cofacets[cofacetidx].first;
            
            auto& btable = matching_context.binomial_table;
            auto npts = btable.size();
            auto fvts = SimplexUtility::getSimplexVertices(btable, fbindex, npts, 2);
            auto cfvts = SimplexUtility::getSimplexVertices(btable, cfbindex, npts, 3);

            std::cout<<"birth facet vertices =  ";
            for(auto i : fvts) std::cout<<i<<"  ";
            std::cout<<'\n';

            std::cout<<"death cofacet vertices = ";
            for (auto i : cfvts) std::cout<<i<<"  ";
            std::cout<<'\n';
        }

    }

    return count;
}


int64_t MaximumMorseMatching::serialFacetMatchAndGetMinCriticalIndex(MatchingContext& matching_context)
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
            double facetweight = matching_context.sorted_facets[i].second;
            std::cout<<"facet weight = "<<facetweight<<"  cofacet weight = -1"<<'\n';
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
            std::cout<<"maxdim pair:  facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';
        }

    }

    return minindex;
}

std::vector< std::vector<size_t> > MaximumMorseMatching::serialFacetMatchAndGetAugPath(MatchingContext& matching_context)
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

        if (augpathlen < 0)
        {
            double facetweight = matching_context.sorted_facets[i].second;
            std::cout<<"facet weight = "<<facetweight<<"  cofacet weight = -1"<<'\n';
            continue;
        }

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
            std::cout<<"maxdim pair:  facet weight = "<<facetweight<<"  cofacet weight = "<<cofacetweight<<'\n';
        }
    }

    return aug_path_cofacet_vec;
}