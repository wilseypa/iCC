#include <execution>
#include <iostream>

#include "omp.h"

#include "MaximumMorseMatching.hpp"


void MaximumMorseMatching::parallelMaxFacetInit(MatchingContext& matching_context)
{
    auto graph = matching_context.graph;
    //convert index
    auto u = matching_context.graph.u_nodes;
    auto v = matching_context.graph.v_nodes;

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
    auto u = matching_context.graph.u_nodes;
    auto v = matching_context.graph.v_nodes;

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


