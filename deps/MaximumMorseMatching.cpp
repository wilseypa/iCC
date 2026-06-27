#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <cassert>

#include "omp.h"

#include "MaximumMorseMatching.hpp"
#include "SimplexUtility.hpp"

namespace
{
    void printSimplexVertices(std::ostream& os,
                            const std::vector<std::vector<int64_t>>& binomial_table,
                            const int64_t bindex,
                            const size_t npts,
                            const size_t simplex_dim)
    {
        if (bindex < 0)
        {
            os << "{?}";
            return;
        }

        const auto simplex_vertices =
            SimplexUtility::getSimplexVertices(binomial_table, bindex, npts, simplex_dim);

        os << '{';
        for (size_t i = 0; i < simplex_vertices.size(); ++i)
        {
            if (i > 0) os << ',';
            os << simplex_vertices[i];
        }
        os << '}';
    }

    void printGhostPair(const std::vector<std::vector<int64_t>>& binomial_table,
                        const size_t npts,
                        const size_t cofacet_dim,
                        const double facetweight,
                        const double cofacetweight,
                        const int64_t facetbindex,
                        const int64_t cofacetbindex)
    {
        std::cerr << "[ghost pair] cofacet weight != facet weight\n"
                << "  cofacet dim = " << cofacet_dim << '\n'
                << "  facet weight = " << facetweight
                << ", bindex = " << facetbindex
                << ", vertices = ";
        printSimplexVertices(std::cerr,
                            binomial_table,
                            facetbindex,
                            npts,
                            cofacet_dim == 0 ? 0 : cofacet_dim - 1);
        std::cerr << '\n'
                << "  cofacet weight = " << cofacetweight
                << ", bindex = " << cofacetbindex
                << ", vertices = ";
        printSimplexVertices(std::cerr, binomial_table, cofacetbindex, npts, cofacet_dim);
        std::cerr << '\n';
    }
}    // helper namespace

size_t MaximumMorseMatching::implicitMatch(MatchingContext& matching_context,
                                           std::vector<std::pair<double, double>>& dim_persistent_pair,
                                           std::vector<PersistentPairInfo>* persistent_pair_info)
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

    const auto appendPersistentPair =
        [&](const double facetweight, const double cofacetweight, const int64_t facetbindex, const int64_t cofacetbindex)
    {
        dim_persistent_pair.emplace_back(facetweight, cofacetweight);

        if (persistent_pair_info != nullptr)
        {
            persistent_pair_info->push_back(PersistentPairInfo{facetweight, cofacetweight, facetbindex, cofacetbindex});
        }
    };

    //process facet in reverse order
    for (int64_t i = static_cast<int64_t>(v) - 1; i >= 0; --i)
    {
        const size_t facet_list_index = static_cast<size_t>(i);
        const int64_t facetbindex = facet_list[facet_list_index].first;

        //skip non active facets
        auto fit = facet_hash.find(facetbindex);
        if (fit == facet_hash.end()) continue;

        //covert list index to graph index
        size_t facetgraphidx = facet_list_index + u;

        SimplexUtility::getCofacetListIndicesInPlace(binom_table, cofacet_hash, cofacet_indices_, vertex_workspace_, facetbindex, npts, dim-1);

        if (cofacet_indices_.empty())
        {
            const double facetweight = matching_context.sorted_facets[facet_list_index].second;
            // std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            appendPersistentPair(facetweight, -1.0, facetbindex, -1);
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
                if (facet_list_index == maxfacetlistidx)
                {
                    const double cofacetweight = cofacet_list[mincofacetidx].second;
                    const double facetweight = facet_list[facet_list_index].second;

                    if (cofacetweight == facetweight)
                    {
                        //apparent pair
                        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(mincofacetidx);
                        bi_graph.match_list[mincofacetidx] = static_cast<int64_t>(facetgraphidx);

                        ct+=1;

                        continue;
                    }

                    // printGhostPair(binom_table, npts, dim, facetweight, cofacetweight, facetbindex, cofacet_list[mincofacetidx].first);
                }
            }
        }

        const int64_t terminalcofacet = implicitFacetCompressedAugPath(binom_table, bi_graph, facet_list,
                                                                       cofacet_hash, cofacet_indices_, npts, dim);

        if (terminalcofacet < 0)
        {
            const double facetweight = matching_context.sorted_facets[facet_list_index].second;
            // std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = -1 "<<'\n';
            appendPersistentPair(facetweight, -1.0, facetbindex, -1);
            count += 1;
            continue;
        }

        const size_t uidx = static_cast<size_t>(terminalcofacet);
        const double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        const double facetweight = matching_context.sorted_facets[facet_list_index].second;
        const int64_t reduced_facet_bindex = matching_context.sorted_facets[facet_list_index].first;
        const int64_t reduced_cofacet_bindex = matching_context.sorted_cofacets[uidx].first;

        // compressed path: pair the current facet directly with the terminal cofacet.
        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(uidx);
        bi_graph.match_list[uidx] = static_cast<int64_t>(facetgraphidx);

        if (facetweight != cofacetweight)
        {
            appendPersistentPair(facetweight, cofacetweight, reduced_facet_bindex, reduced_cofacet_bindex);
            // std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }

    }

    return count;
}

MaximumMorseMatching::MatchSupportInfo MaximumMorseMatching::implicitMatchAndCollectSupportInfo(
    MatchingContext& matching_context,
    std::vector<std::pair<double, double>>& dim_persistent_pair,
    const bool collect_pv_support,
    std::vector<PersistentPairInfo>* persistent_pair_info)
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

    MatchSupportInfo match_support_info;
    if (collect_pv_support)
        match_support_info.raw_pv_support_cofacet_indices.reserve(dim < 5 ? 64 : 128);    //estimated numbers
    match_support_info.protected_facet_list_indices.reserve(dim < 5 ? 32 : 64);

    // size_t ct = 0;    //apparent pair count

    const auto appendPersistentPair =
        [&](const double facetweight, const double cofacetweight, const int64_t facetbindex, const int64_t cofacetbindex)
    {
        dim_persistent_pair.emplace_back(facetweight, cofacetweight);

        if (persistent_pair_info != nullptr)
        {
            persistent_pair_info->push_back(PersistentPairInfo{
                facetweight,
                cofacetweight,
                facetbindex,
                cofacetbindex,
            });
        }
    };

    const auto appendProtectedFacet =
        [&](const size_t facet_list_index)
    {
        match_support_info.protected_facet_list_indices.push_back(facet_list_index);
    };
    
    // process facets in reverse order
    for (int64_t i = static_cast<int64_t>(v) - 1; i >= 0; --i)
    {
        const size_t facet_list_index = static_cast<size_t>(i);
        const int64_t facetbindex = facet_list[facet_list_index].first;

        // skip non-active facets
        auto fit = facet_hash.find(facetbindex);
        if (fit == facet_hash.end()) continue;

        const size_t facetgraphidx = facet_list_index + u;

        // compute immediate cofacets of this facet (list indices)
        SimplexUtility::getCofacetListIndicesInPlace(binom_table, cofacet_hash, cofacet_indices_, vertex_workspace_,
                                                     facetbindex, npts, dim - 1);

        //if active and have no cofacets 
        if (cofacet_indices_.empty())
        {
            const double facetweight = matching_context.sorted_facets[facet_list_index].second;
            appendPersistentPair(facetweight, -1.0, facetbindex, -1);
            appendProtectedFacet(facet_list_index);

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

                if (facet_list_index == maxfacetlistidx)
                {
                    const double cofacetweight = cofacet_list[mincofacetidx].second;
                    const double facetweight = facet_list[facet_list_index].second;

                    if (cofacetweight == facetweight)
                    {
                        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(mincofacetidx);
                        bi_graph.match_list[mincofacetidx] = static_cast<int64_t>(facetgraphidx);

                        continue;
                    }

                    // printGhostPair(binom_table, npts, dim, facetweight, cofacetweight, facetbindex, cofacet_list[mincofacetidx].first);
                }
            }
        }

        const int64_t terminalcofacet = implicitFacetCompressedAugPath(binom_table, bi_graph, facet_list,
                                                                       cofacet_hash, cofacet_indices_, npts, dim);

        if (terminalcofacet < 0)
        {
            const double facetweight = matching_context.sorted_facets[facet_list_index].second;
            appendPersistentPair(facetweight, -1.0, facetbindex, -1);
            appendProtectedFacet(facet_list_index);

            continue;
        }

        const size_t terminalcofacetidx = static_cast<size_t>(terminalcofacet);
        if (collect_pv_support)
        {
            // PV detection needs the reduced column before the compressed pair is installed.
            match_support_info.raw_pv_support_cofacet_indices.push_back(collectReducedColumnSupport(matching_context, terminalcofacetidx, facetgraphidx));
        }

        const size_t uidx = terminalcofacetidx;
        const double cofacetweight = matching_context.sorted_cofacets[uidx].second;
        const double facetweight = matching_context.sorted_facets[facet_list_index].second;
        const int64_t reduced_facet_bindex = matching_context.sorted_facets[facet_list_index].first;
        const int64_t reduced_cofacet_bindex = matching_context.sorted_cofacets[uidx].first;

        // Compressed install: pair the current facet directly with the terminal cofacet.
        bi_graph.match_list[facetgraphidx] = static_cast<int64_t>(uidx);
        bi_graph.match_list[uidx] = static_cast<int64_t>(facetgraphidx);

        if (facetweight != cofacetweight)
        {
            appendPersistentPair(facetweight, cofacetweight, reduced_facet_bindex, reduced_cofacet_bindex);
            // std::cout <<"interface dim = "<<dim<< "  facet weight = " << facetweight << "  cofacet weight = " << cofacetweight << '\n';
        }

    }

    // std::cout<<"interface dim = "<<dim<<" implicit apparent pair count = "<<ct<<'\n';
    
    return match_support_info;
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
        size_t topcofacet = cofacet_queue.top();  //cofacet's graph index == list index
        cofacet_queue.pop();

        // std::cout<<"current top cofacet of queue = "<<topcofacet<<'\n';

        // count multiplicity of the same cofacet at the top of the queue
        size_t multiplicity = 1;

        while (!cofacet_queue.empty() && cofacet_queue.top() == topcofacet)
        {
            cofacet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL) continue;

        if (bi_graph.match_list[topcofacet] < 0)
        {
            unmatchedcofacet = topcofacet;
            break;
        }

        size_t nextfacet = bi_graph.match_list[topcofacet];  //graph index

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

    size_t currenttopcofacet = aug_path_.front();

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

int64_t MaximumMorseMatching::implicitFacetCompressedAugPath(const std::vector<std::vector<int64_t>>& binomial_table,
                                                             const BipartiteGraph& bi_graph,
                                                             const std::vector<std::pair<int64_t, double>>& facet_list,
                                                             const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table,
                                                             const std::vector<size_t>& start_cofacet_indices,
                                                             size_t npts,
                                                             size_t interfacedimension)
{
    pq_workspace_.clear();

    MinIndexQueue cofacet_queue(std::greater<size_t>{}, std::move(pq_workspace_));

    for (auto cofidx : start_cofacet_indices)
    {
        cofacet_queue.push(cofidx);
    }

    while (!cofacet_queue.empty())
    {
        const size_t topcofacet = cofacet_queue.top();  // cofacet graph index == list index
        cofacet_queue.pop();

        size_t multiplicity = 1;
        while (!cofacet_queue.empty() && cofacet_queue.top() == topcofacet)
        {
            cofacet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL) continue;

        if (bi_graph.match_list[topcofacet] < 0)
        {
            auto& pq_buffer = cofacet_queue.getContainer();
            pq_buffer.clear();
            pq_workspace_ = std::move(pq_buffer);
            return static_cast<int64_t>(topcofacet);
        }

        const size_t nextfacet = static_cast<size_t>(bi_graph.match_list[topcofacet]);
        const size_t nextfacet_list_idx = nextfacet - bi_graph.unodes;
        enqueueReducedCompressedColumnTail(binomial_table,
                                           bi_graph,
                                           facet_list,
                                           cofacet_hash_table,
                                           nextfacet_list_idx,
                                           topcofacet,
                                           npts,
                                           interfacedimension,
                                           cofacet_queue);
    }

    auto& pq_buffer = cofacet_queue.getContainer();
    pq_buffer.clear();
    pq_workspace_ = std::move(pq_buffer);

    return -1;
}

void MaximumMorseMatching::enqueueReducedCompressedColumnTail(const std::vector<std::vector<int64_t>>& binomial_table,
                                                              const BipartiteGraph& bi_graph,
                                                              const std::vector<std::pair<int64_t, double>>& facet_list,
                                                              const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table,
                                                              const size_t facet_list_index,
                                                              const size_t expected_pivot_cofacet,
                                                              size_t npts,
                                                              size_t interfacedimension,
                                                              MinIndexQueue& target_queue)
{
    const int64_t facetbindex = facet_list[facet_list_index].first;
    SimplexUtility::getCofacetListIndicesInPlace(binomial_table,
                                                 cofacet_hash_table,
                                                 cofacet_indices_,
                                                 vertex_workspace_,
                                                 facetbindex,
                                                 npts,
                                                 interfacedimension - 1);

    const auto min_cofacet_iter = std::min_element(cofacet_indices_.begin(), cofacet_indices_.end());
    if (min_cofacet_iter != cofacet_indices_.end() && *min_cofacet_iter == expected_pivot_cofacet)
    {
        for (const auto cofidx : cofacet_indices_)
        {
            if (cofidx != expected_pivot_cofacet)
                target_queue.push(cofidx);
        }
        return;
    }

    MinIndexQueue column_queue(std::greater<size_t>{});

    for (auto cofidx : cofacet_indices_)
    {
        column_queue.push(cofidx);
    }

    while (!column_queue.empty())
    {
        const size_t topcofacet = column_queue.top();
        column_queue.pop();

        size_t multiplicity = 1;
        while (!column_queue.empty() && column_queue.top() == topcofacet)
        {
            column_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL) continue;

        if (topcofacet == expected_pivot_cofacet)
        {
            // Expected pivot exposed. Do not enqueue it; this is the term being canceled.
            break;
        }

        if (topcofacet > expected_pivot_cofacet)
        {
            // Since the queue is a min-queue, we have already passed the expected pivot.
            // The pivot cannot reappear under monotone reduced-tail expansion.
            throw std::logic_error("enqueueReducedCompressedColumnTail: passed expected pivot");
        }

        // Now topcofacet < expected_pivot_cofacet. It must be reducible by an already processed matched facet.
        const int64_t matchedfacet = bi_graph.match_list[topcofacet];
        if (matchedfacet < 0)
        {
            throw std::logic_error("enqueueReducedCompressedColumnTail: unmatched term before expected pivot");
        }

        const size_t matched_facet_list_index = static_cast<size_t>(matchedfacet) - bi_graph.unodes;
        if (matched_facet_list_index <= facet_list_index)
        {
            throw std::logic_error("enqueueReducedCompressedColumnTail: reducer is not available in processing order");
        }

        enqueueReducedCompressedColumnTail(binomial_table,
                                           bi_graph,
                                           facet_list,
                                           cofacet_hash_table,
                                           matched_facet_list_index,
                                           topcofacet,
                                           npts,
                                           interfacedimension,
                                           column_queue);
        continue;
    }

    while (!column_queue.empty())
    {
        const size_t queued_cofacet = column_queue.top();
        column_queue.pop();

        size_t multiplicity = 1;
        while (!column_queue.empty() && column_queue.top() == queued_cofacet)
        {
            column_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) != 0ULL && queued_cofacet != expected_pivot_cofacet)
            target_queue.push(queued_cofacet);
    }
}

std::vector<size_t> MaximumMorseMatching::collectReducedColumnSupport(const MatchingContext& matching_context, const size_t terminalcofacet, const size_t expected_facet_graph_idx)
{
    auto& binom_table = matching_context.binomial_table;
    auto& cofacet_list = matching_context.sorted_cofacets;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& bi_graph = matching_context.graph;
    const size_t u = bi_graph.unodes;
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

    for (size_t fi : facet_indices_) facet_queue.push(fi);  //list indices

    while (!facet_queue.empty())
    {
        size_t top_facet_list_idx = facet_queue.top();  //list index
        facet_queue.pop();

        size_t multiplicity = 1;
        while (!facet_queue.empty() && top_facet_list_idx == facet_queue.top())
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL) continue;

        const int64_t nextcofacet = bi_graph.match_list[top_facet_list_idx + u];
        if (nextcofacet < 0) 
        {
            assert(expected_facet_graph_idx == (top_facet_list_idx + u));    //debug purpose
            break;
        }

        const size_t cfi = static_cast<size_t>(nextcofacet);
        cofacet_trace.push_back(cfi);

        enqueueReducedCofacetBoundaryTail(matching_context,
                                          cfi,
                                          top_facet_list_idx,
                                          facet_queue,
                                          cofacet_trace);
    }

    auto& pq_buffer = facet_queue.getContainer();
    pq_buffer.clear();
    pq_workspace_ = std::move(pq_buffer);

    return cofacet_trace;
}

void MaximumMorseMatching::enqueueReducedCofacetBoundaryTail(const MatchingContext& matching_context,
                                                             const size_t cofacet_list_index,
                                                             const size_t pivot_facet_list_index,
                                                             MaxIndexQueue& target_queue,
                                                             std::vector<size_t>& cofacet_trace)
{
    auto& binom_table = matching_context.binomial_table;
    auto& cofacet_list = matching_context.sorted_cofacets;
    auto& facet_hash = matching_context.facet_bindex_to_list_index;
    auto& bi_graph = matching_context.graph;
    const size_t u = bi_graph.unodes;
    const size_t dim = matching_context.dim;
    const size_t npts = matching_context.npts;

    SimplexUtility::getFacetListIndicesInPlace(binom_table,
                                               facet_hash,
                                               facet_indices_,
                                               vertex_workspace_,
                                               cofacet_list[cofacet_list_index].first,
                                               npts,
                                               dim);

    const auto max_facet_iter = std::max_element(facet_indices_.begin(), facet_indices_.end());
    if (max_facet_iter != facet_indices_.end() && *max_facet_iter == pivot_facet_list_index)
    {
        for (const auto fi : facet_indices_)
        {
            if (fi != pivot_facet_list_index)
                target_queue.push(fi);
        }
        return;
    }

    MaxIndexQueue facet_queue(std::less<size_t>{});

    for (size_t fi : facet_indices_)
    {
        facet_queue.push(fi);
    }

    while (!facet_queue.empty())
    {
        const size_t top_facet_list_idx = facet_queue.top();
        facet_queue.pop();

        size_t multiplicity = 1;
        while (!facet_queue.empty() && top_facet_list_idx == facet_queue.top())
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL) continue;

        const int64_t nextcofacet = bi_graph.match_list[top_facet_list_idx + u];
        if (nextcofacet >= 0 && top_facet_list_idx > pivot_facet_list_index)
        {
            const size_t cfi = static_cast<size_t>(nextcofacet);
            cofacet_trace.push_back(cfi);

            enqueueReducedCofacetBoundaryTail(matching_context,
                                              cfi,
                                              top_facet_list_idx,
                                              facet_queue,
                                              cofacet_trace);
            continue;
        }

        if (top_facet_list_idx != pivot_facet_list_index)
        {
            assert(false && "cofacet boundary replay did not stop at its pivot facet");
            target_queue.push(top_facet_list_idx);
        }

        break;
    }

    while (!facet_queue.empty())
    {
        const size_t queued_facet = facet_queue.top();
        facet_queue.pop();

        size_t multiplicity = 1;
        while (!facet_queue.empty() && facet_queue.top() == queued_facet)
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) != 0ULL && queued_facet != pivot_facet_list_index)
            target_queue.push(queued_facet);
    }
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


int64_t MaximumMorseMatching::implicitFacetAugPathDebug(
    const std::vector<std::vector<int64_t>>& binomial_table,
    BipartiteGraph& bi_graph,
    const std::vector<std::pair<int64_t, double>>& facet_list,
    const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table,
    const size_t facetgraphindex,
    size_t npts,
    size_t interfacedimension)
{
    aug_path_.clear();
    cofacet_stack_.clear();
    pq_workspace_.clear();

    MinIndexQueue cofacet_queue(std::greater<size_t>{}, std::move(pq_workspace_));

    int64_t unmatchedcofacet = -1;

    // set to false when you want to silence the debug prints
    // bool debug_multiplicity = (interfacedimension == 3); 
    bool debug_multiplicity = false;

    const size_t start_facet_list_idx = facetgraphindex - bi_graph.unodes;
    const int64_t start_facet_bindex = facet_list[start_facet_list_idx].first;
    const double start_facet_weight = facet_list[start_facet_list_idx].second;

    // histogram for this call
    size_t mult1_count = 0;
    size_t mult2_count = 0;
    size_t mult3_count = 0;
    size_t mult4plus_count = 0;

    auto printSimplexVertices = [&](int64_t bindex, size_t simplex_dim)
    {
        if (bindex < 0)
        {
            std::cerr << "{?}";
            return;
        }

        const auto simplex_vertices =
            SimplexUtility::getSimplexVertices(binomial_table, bindex, npts, simplex_dim);

        std::cerr << '{';
        for (size_t i = 0; i < simplex_vertices.size(); ++i)
        {
            if (i > 0) std::cerr << ',';
            std::cerr << simplex_vertices[i];
        }
        std::cerr << '}';
    };

    // reverse lookup only for debug printing; O(n), but used only on rare events
    auto findCofacetBindex = [&](size_t target_list_idx) -> int64_t
    {
        for (const auto& kv : cofacet_hash_table)
        {
            if (kv.second == target_list_idx) return kv.first;
        }
        return -1;
    };

    auto printMultiplicitySummary = [&](bool path_found)
    {
        if (!debug_multiplicity) return;
        if (mult3_count == 0 && mult4plus_count == 0) return;

        std::cerr << "\n[DMT debug] implicitFacetAugPath summary\n"
                  << "  interface dim = " << interfacedimension << '\n'
                  << "  start facet list idx = " << start_facet_list_idx
                  << ", graph idx = " << facetgraphindex
                  << ", bindex = " << start_facet_bindex
                  << ", weight = " << start_facet_weight
                  << ", vertices = ";
        printSimplexVertices(start_facet_bindex, interfacedimension - 1);
        std::cerr << '\n'
                  << "  multiplicity histogram: "
                  << "1->" << mult1_count
                  << ", 2->" << mult2_count
                  << ", 3->" << mult3_count
                  << ", >=4->" << mult4plus_count << '\n'
                  << "  augmenting path found = " << (path_found ? "yes" : "no") << '\n';
    };

    // reuse the immediate cofacets computed for apparent pair check
    for (auto cofidx : cofacet_indices_)
    {
        cofacet_queue.push(cofidx);
    }

    while (!cofacet_queue.empty())
    {
        size_t topcofacet = cofacet_queue.top();  // cofacet graph index == list index
        cofacet_queue.pop();

        // count multiplicity of the same cofacet at the top of the queue
        size_t multiplicity = 1;
        while (!cofacet_queue.empty() && cofacet_queue.top() == topcofacet)
        {
            cofacet_queue.pop();
            ++multiplicity;
        }

        if (multiplicity == 1) ++mult1_count;
        else if (multiplicity == 2) ++mult2_count;
        else if (multiplicity == 3) ++mult3_count;
        else ++mult4plus_count;

        if (debug_multiplicity && multiplicity >= 3)
        {
            const int64_t topcofacet_bindex = findCofacetBindex(topcofacet);

            std::cerr << "\n[DMT debug] implicitFacetAugPath saw repeated cofacet\n"
                      << "  interface dim = " << interfacedimension << '\n'
                      << "  multiplicity = " << multiplicity << '\n'
                      << "  start facet list idx = " << start_facet_list_idx
                      << ", graph idx = " << facetgraphindex
                      << ", bindex = " << start_facet_bindex
                      << ", weight = " << start_facet_weight
                      << ", vertices = ";
            printSimplexVertices(start_facet_bindex, interfacedimension - 1);
            std::cerr << '\n'
                      << "  repeated cofacet list idx = " << topcofacet
                      << ", bindex = " << topcofacet_bindex
                      << ", vertices = ";
            printSimplexVertices(topcofacet_bindex, interfacedimension);
            std::cerr << '\n'
                      << "  current match(topcofacet) = " << bi_graph.match_list[topcofacet] << '\n';

            if (bi_graph.match_list[topcofacet] >= 0)
            {
                const size_t matchedfacet_graph_idx =
                    static_cast<size_t>(bi_graph.match_list[topcofacet]);
                const size_t matchedfacet_list_idx =
                    matchedfacet_graph_idx - bi_graph.unodes;

                std::cerr << "  matched facet list idx = " << matchedfacet_list_idx
                          << ", graph idx = " << matchedfacet_graph_idx
                          << ", bindex = " << facet_list[matchedfacet_list_idx].first
                          << ", weight = " << facet_list[matchedfacet_list_idx].second
                          << ", vertices = ";
                printSimplexVertices(facet_list[matchedfacet_list_idx].first,
                                     interfacedimension - 1);
                std::cerr << '\n';
            }

            std::cerr << "  queue size after consolidation = "
                      << cofacet_queue.size() << '\n';
        }

        // current behavior: skip all duplicates to preserve single-path traversal
        // if (multiplicity > 1) continue;

        // parity-test alternative:
        if ((multiplicity & 1ULL) == 0ULL) continue;

        if (bi_graph.match_list[topcofacet] < 0)
        {
            unmatchedcofacet = static_cast<int64_t>(topcofacet);
            break;
        }

        size_t nextfacet = static_cast<size_t>(bi_graph.match_list[topcofacet]);  // graph index

        // get cofacets of the nextfacet
        int64_t facetbindex = facet_list[nextfacet - bi_graph.unodes].first;
        SimplexUtility::getCofacetListIndicesInPlace(
            binomial_table,
            cofacet_hash_table,
            cofacet_indices_,
            vertex_workspace_,
            facetbindex,
            npts,
            interfacedimension - 1);

        for (auto cofidx : cofacet_indices_)
        {
            if (cofidx != topcofacet) cofacet_queue.push(cofidx);
        }

        cofacet_stack_.push_back(topcofacet);
    }

    // no aug path found
    if (unmatchedcofacet < 0)
    {
        printMultiplicitySummary(false);

        auto& pq_buffer = cofacet_queue.getContainer();
        pq_buffer.clear();
        pq_workspace_ = std::move(pq_buffer);
        return -1;
    }

    // reconstruct the actual path
    aug_path_.push_back(static_cast<size_t>(unmatchedcofacet));

    size_t currenttopcofacet = aug_path_.front();

    // backtrack through the cofacet stack
    for (auto it = cofacet_stack_.rbegin(); it != cofacet_stack_.rend(); ++it)
    {
        size_t matchedfacet = static_cast<size_t>(bi_graph.match_list[*it]);

        // check if the matchedfacet is adjacent to the current top cofacet on the path
        int64_t facetbindex = facet_list[matchedfacet - bi_graph.unodes].first;
        SimplexUtility::getCofacetListIndicesInPlace(
            binomial_table,
            cofacet_hash_table,
            cofacet_indices_,
            vertex_workspace_,
            facetbindex,
            npts,
            interfacedimension - 1);

        if (std::find(cofacet_indices_.begin(), cofacet_indices_.end(), currenttopcofacet)
            != cofacet_indices_.end())
        {
            aug_path_.push_back(matchedfacet);
            aug_path_.push_back(*it);
            currenttopcofacet = *it;
        }
    }

    aug_path_.push_back(facetgraphindex);    // initial facet index

    printMultiplicitySummary(true);

    // move the container back
    auto& pq_buffer = cofacet_queue.getContainer();
    pq_buffer.clear();
    pq_workspace_ = std::move(pq_buffer);

    return static_cast<int64_t>(aug_path_.size());
}
