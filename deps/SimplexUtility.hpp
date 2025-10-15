#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <numeric> // For std::iota

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
        // assume j > i
        if (j < i)
            std::swap(i, j);

        return (binomial_table[j][2] + i);
    }

    inline size_t getSimplexMaxVertex(const std::vector<std::vector<int64_t>> &binomial_table, const int64_t bindex, size_t vtnum, const size_t dim)
    {
        // top = vtnum - 1, bottom = dim
        size_t top = vtnum - 1;
        if (binomial_table[top][dim + 1] > bindex)
        {
            size_t span = top - dim;
            while (span > 0)
            {
                size_t half = (span >> 1);
                size_t mid = top - half;
                if (binomial_table[mid][dim + 1] > bindex)
                {
                    top = mid - 1;
                    span -= (half + 1);
                }
                else
                    span = half;
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

    inline void sortSimplexByWeightThenIndex(std::vector<std::pair<int64_t, double>> &simplex_list)
    {
        auto sort_lambda = [](const auto &lhs, const auto &rhs)
        { return (lhs.second < rhs.second) || ((lhs.second == rhs.second) && (lhs.first < rhs.first)); };

        std::sort(simplex_list.begin(), simplex_list.end(), sort_lambda);

        return;
    }

    inline std::vector<int64_t> getFacetBinomialIndex(const std::vector<std::vector<int64_t>> &binomial_table, const int64_t binomindex, const size_t dim)
    {
        //{i, j, k} -> {j, k}, {i, k}, {i, j}. i > j > k
        std::vector<int64_t> facet_bindex;
        facet_bindex.reserve(dim + 1);

        size_t npt = binomial_table.size() - 1;

        std::vector<size_t> simplex_vt = getSimplexVertices(binomial_table, binomindex, npt, dim);

        int64_t above = 0;
        int64_t below = binomindex;
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

}