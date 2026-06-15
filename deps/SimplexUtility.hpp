#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <numeric> // For std::iota
#include <utility> // For std::swap

#include "robin_hood.h"

namespace SimplexUtility
{
    inline std::vector<std::vector<int64_t>> getBinomialTable(const size_t vtnum, const size_t maxdim)
    {
        // binom table by rows. same as lhf
        std::vector<std::vector<int64_t>> binom_table(vtnum + 1, std::vector<int64_t>(maxdim + 1 + 1, 0));

        binom_table[0][0] = 1;
        for (size_t i = 1; i < vtnum + 1; i++)
        {
            binom_table[i][0] = 1;
            for (size_t j = 1; j < (maxdim + 2); j++)
            {
                binom_table[i][j] = binom_table[i - 1][j - 1] + binom_table[i - 1][j];
                if (binom_table[i][j] < 0)
                    throw std::overflow_error("Binomial Overflow Error!");
            }
        }

        return binom_table;
    }

    inline int64_t getBinomialIndex(const std::vector<std::vector<int64_t>> &binomial_table, const std::vector<size_t> &simplex_vt, const size_t shift) // shift is used for cofacet enum
    {
        int64_t bindex = 0;
        size_t k = simplex_vt.size();

        // simplex_vt is sorted in descending order
        for (auto it = simplex_vt.begin(); it != simplex_vt.end(); it++)
        {
            bindex += binomial_table[*it][k + shift];
            if (bindex < 0)
                throw std::overflow_error("Binomial overflow in get simplex binomial index.");
            k--;
        }

        return bindex;
    }

    inline int64_t getEdgeBinomialIndex(const std::vector<std::vector<int64_t>> &binomial_table, size_t j, size_t i)
    {
        //simplex vts are sorted in descending order, edge (j, i),  j > i 
        return (j > i) ? (binomial_table[j][2] + i) : (binomial_table[i][2] + j);
    }

    inline size_t getSimplexMaxVertex(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t bindex, const size_t vtnum, const size_t dim)
    {
        // top = vtnum - 1, bottom = dim
        size_t top = vtnum - 1;
        if (binomial_table[top][dim + 1] > bindex)
        {
            size_t span = top - dim;
            while (span > 0)
            {
                //mid = top - floor(span/2) = top - floor((top - dim)/2) = ceil((dim + top)/2)
                //mid = dim + ((top - dim + 1) >> 1) or = top - (span >> 1)  avoid overflow
                size_t half = (span >> 1);
                size_t mid = top - half;
                if (binomial_table[mid][dim + 1] > bindex)
                {
                    top = mid - 1;
                    span -= (half + 1);
                }
                else
                {
                    span = half;
                }
            }
        }
        return top;
    }

    inline std::vector<size_t> getSimplexVertices(const std::vector<std::vector<int64_t>> &binomial_table, int64_t bindex, size_t vtnum, const size_t dim)
    {
        // max vertex index = vtnum - 1
        std::vector<size_t> simplex_vt;
        simplex_vt.reserve(dim + 1);

        for (size_t k = dim + 1; k > 0; k--) // size_t is unsigned. cannot use it to check against negative number
        {
            vtnum = getSimplexMaxVertex(binomial_table, bindex, vtnum, k - 1);

            bindex -= binomial_table[vtnum][k];

            simplex_vt.push_back(vtnum);
        }
        return simplex_vt;
    }


    //get simplex vertices in pre-allocated space
    inline void getSimplexVerticesInPlace(const std::vector<std::vector<int64_t>>& binomial_table, std::vector<size_t>& output_workspace, int64_t bindex, size_t vtnum, const size_t dim)
    {
        output_workspace.clear();

        for (size_t k = dim + 1; k > 0; k--)  //size_t is unsigned. cannot use it to check against negative number
        {
            vtnum = getSimplexMaxVertex(binomial_table, bindex, vtnum, k - 1);

            bindex -= binomial_table[vtnum][k];

            output_workspace.push_back(vtnum);
        }

        return;
    }


    inline void sortSimplexByWeightThenIndex(std::vector<std::pair<int64_t, double>> &simplex_list)
    {
        auto sort_lambda = [](const auto &lhs, const auto &rhs)
        { return (lhs.second < rhs.second) || ((lhs.second == rhs.second) && (lhs.first < rhs.first)); };

        std::sort(simplex_list.begin(), simplex_list.end(), sort_lambda);

        return;
    }

    inline std::vector<int64_t> getFacetBinomialIndices(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t bindex, const size_t dim)
    {
        //{i, j, k} -> {j, k}, {i, k}, {i, j}. i > j > k
        std::vector<int64_t> facet_bindex;
        facet_bindex.reserve(dim + 1);

        size_t npts = binomial_table.size() - 1;

        std::vector<size_t> simplex_vt = getSimplexVertices(binomial_table, bindex, npts, dim);

        int64_t above = 0;
        int64_t below = bindex;
        size_t k = dim;

        for (auto i = 0; i < simplex_vt.size(); i++)
        {
            size_t vt = simplex_vt[i];
            below -= binomial_table[vt][k + 1];

            facet_bindex.push_back(above + below);

            above += binomial_table[vt][k];

            k--;
        }

        return facet_bindex;
    }

    inline size_t mstFindRoot(std::vector<size_t> &parent_idx, size_t x)
    {
        while (parent_idx[x] != x)
        {
            parent_idx[x] = parent_idx[parent_idx[x]];
            x = parent_idx[x];
        }

        return x;
    }

    inline void mstSetUnion(std::vector<size_t> &parent_idx, size_t x, size_t y)
    {
        size_t p = mstFindRoot(parent_idx, x);
        size_t q = mstFindRoot(parent_idx, y);
        parent_idx[p] = parent_idx[q];
        return;
    }

    inline robin_hood::unordered_map<int64_t, size_t> getActiveEdgeIndexHashTable(const std::vector<std::vector<int64_t>> &binomial_table, const std::vector<std::pair<int64_t, double>> &sorted_edge, const size_t npts)
    {
        robin_hood::unordered_map<int64_t, size_t> active_edge_index_hash;
        active_edge_index_hash.reserve(sorted_edge.size());

        std::vector<size_t> parent_idx(npts);
        std::iota(parent_idx.begin(), parent_idx.end(), 0);

        size_t x, y;
        int64_t bindex;
        size_t m = sorted_edge.size();
        std::vector<size_t> edge_vt;
        for (auto i = 0; i < m; i++)
        {
            bindex = sorted_edge[i].first;
            edge_vt = getSimplexVertices(binomial_table, bindex, npts, 1);
            x = edge_vt[0];
            y = edge_vt[1];
            if (mstFindRoot(parent_idx, x) != mstFindRoot(parent_idx, y))
            {
                mstSetUnion(parent_idx, x, y);
                continue;
            }
            active_edge_index_hash.emplace(bindex, i);
        }

        return active_edge_index_hash;
    }

    inline robin_hood::unordered_map<int64_t, size_t> getActiveSimplexIndexHashTable(const std::vector<int64_t> &graph_match_list, const std::vector<std::pair<int64_t, double>> &facet_list)
    {
        robin_hood::unordered_map<int64_t, size_t> active_facet_index_hash;
        active_facet_index_hash.reserve(facet_list.size());

        auto facetnum = facet_list.size();

        for (auto i = 0; i < facetnum; ++i)
        {
            if (graph_match_list[i] < 0)
            {
                auto bindex = facet_list[i].first;
                active_facet_index_hash.emplace(bindex, i);
            }
        }

        return active_facet_index_hash;
    }

    inline robin_hood::unordered_map<int64_t, size_t> getSimplexIndexHashTable(const std::vector<std::pair<int64_t, double>> &facet_list)
    {
        robin_hood::unordered_map<int64_t, size_t> facet_index_hash;
        facet_index_hash.reserve(facet_list.size());

        for (auto i = 0; i < facet_list.size(); ++i)
        {
            auto bindex = facet_list[i].first;
            facet_index_hash.emplace(bindex, i);
        }

        return facet_index_hash;
    }

    inline double getVirtualLabelDistance(const robin_hood::unordered_map<uint64_t, double> &virtual_distance_hash_table, size_t i, size_t j)
    {
        if (i > j)
            std::swap(i, j);

        const uint64_t key = (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
        const auto it = virtual_distance_hash_table.find(key);
        return (it != virtual_distance_hash_table.end()) ? it->second : -1.0;
    }

    inline void updateBinomialTable(std::vector<std::vector<int64_t>> &binomial_table, const size_t originalvtnum, const size_t virtualvtnum, const size_t maxdim)
    {
        size_t currentvtnum = binomial_table.size() - 1; // current vertex number in binomial table

        if (currentvtnum >= originalvtnum + virtualvtnum)
            return;

        for (size_t i = 0; i < virtualvtnum; i++)
        {
            binomial_table.emplace_back(maxdim + 1 + 1, 0);
        }

        // std::cout<<"binomial table size = "<<binomial_table.size()<<'\n';

        for (size_t i = 0; i < virtualvtnum; i++)
        {
            binomial_table[currentvtnum + 1 + i][0] = 1;
            for (size_t j = 1; j < (maxdim + 1 + 1); j++)
            {
                binomial_table[currentvtnum + 1 + i][j] = binomial_table[currentvtnum + 1 + i - 1][j - 1] + binomial_table[currentvtnum + 1 + i - 1][j];
                if (binomial_table[i][j] < 0)
                    throw std::overflow_error("Binomial Overflow Error!");
            }
        }

        return;
    }

    //helper function to compute cofacet binomial index
    inline int64_t computeCofacetBindex(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<size_t>& facet_vertices, size_t covt, size_t cofacetdim)
    {
        int64_t bindex = 0;
        size_t k = cofacetdim + 1;
        bool inserted = false;

        for (auto i = 0; i < facet_vertices.size(); ++i)
        {
            //facet_vertices are in descending order
            //insert co-vertex
            if (!inserted && covt > facet_vertices[i])
            {
                bindex += binomial_table[covt][k];
                inserted = true;
                k--;
            }
            bindex += binomial_table[facet_vertices[i]][k];
            k--;
        }

        //co-vertex is the smallest one
        if (!inserted) bindex += binomial_table[covt][k];

        return bindex;
    }


    inline void getCofacetListIndicesInPlace(const std::vector<std::vector<int64_t>>& binomial_table, const robin_hood::unordered_map<int64_t, size_t>& cofacet_hash_table, 
                                  std::vector<size_t>& cofacet_indices, std::vector<size_t>& facet_vertices_workspace, int64_t facetbidx,  size_t npts, size_t facetdim)
    {
        cofacet_indices.clear();

        getSimplexVerticesInPlace(binomial_table, facet_vertices_workspace, facetbidx, npts, facetdim);

        const size_t maxvt = facet_vertices_workspace.front();
        const size_t minvt = facet_vertices_workspace.back();

        //range skipping
        for (size_t covt = 0; covt < minvt; ++covt)
        {
            int64_t cofacetbindex = computeCofacetBindex(binomial_table, facet_vertices_workspace, covt, facetdim + 1);
            auto it = cofacet_hash_table.find(cofacetbindex);
            if (it != cofacet_hash_table.end()) cofacet_indices.push_back(it -> second);
        }

        for (size_t covt = minvt + 1; covt < maxvt; ++covt)
        {
            if (!std::binary_search(facet_vertices_workspace.begin(), facet_vertices_workspace.end(), covt, std::greater<size_t>()))    //descending order
            {
                int64_t cofacetbindex = computeCofacetBindex(binomial_table, facet_vertices_workspace, covt, facetdim + 1);
                auto it = cofacet_hash_table.find(cofacetbindex);
                if (it != cofacet_hash_table.end()) cofacet_indices.push_back(it->second);
            }
        }

        for (size_t covt = maxvt + 1; covt < npts; ++covt)
        {
            int64_t cofacetbindex = computeCofacetBindex(binomial_table, facet_vertices_workspace, covt, facetdim + 1);
            auto it = cofacet_hash_table.find(cofacetbindex);
            if (it != cofacet_hash_table.end()) cofacet_indices.push_back(it -> second);
        }

        return;
    }

    inline void getFacetListIndicesInPlace(const std::vector<std::vector<int64_t>>& binomial_table, const robin_hood::unordered_map<int64_t, size_t>& facet_hash_table, 
                                       std::vector<size_t>& facet_indices, std::vector<size_t>& cofacet_vertices_workspace, int64_t cofacetbindex, size_t npts, size_t cofacetdim)
    {
        facet_indices.clear();

        getSimplexVerticesInPlace(binomial_table, cofacet_vertices_workspace, cofacetbindex, npts, cofacetdim);

        int64_t above = 0;
        int64_t below = cofacetbindex;
        size_t k = cofacetdim;

        for (auto i = 0; i < cofacet_vertices_workspace.size(); i++)
        {
            size_t vt = cofacet_vertices_workspace[i];
            below -= binomial_table[vt][k + 1];

            auto facetbindex = above + below;

            auto it = facet_hash_table.find(facetbindex);
            if (it != facet_hash_table.end()) facet_indices.push_back(it->second);

            above += binomial_table[vt][k];

            k--;
        }

        return;
    }


}
