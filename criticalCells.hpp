#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <chrono>
#include <Eigen/Sparse>

#include <set>
#include <unordered_set>

#include<cstdint>

#include"robin_hood.h"
#include "biGraphMorseMatch.hpp"

#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>


std::ostream &operator<<(std::ostream &os, const std::map<double, std::vector<std::vector<int>>> &bins)
{
    for (const auto &[weight, simplexes] : bins)
    {
        if (!simplexes.empty())
        {
            os << weight << " ";
            for (const auto &simplex : simplexes)
            {
                os << "(";
                for (const auto &i : simplex)
                    os << i << " ";
                os << ")";
            }
            os << std::endl;
        }
    }
    return os;
}

struct SparseDistMat
{
    Eigen::SparseMatrix<double> distMatrix;
    inline double distance(size_t i, size_t j) const { return distMatrix.coeff(i, j); };
};
struct NormallDistMat
{
    std::vector<std::vector<double>> distMatrix;
    inline double distance(size_t i, size_t j) const { return this->distMatrix[i][j]; };
};

/****************************************************** */
//temp dummy struct for alpha complex
//crit cell should inherit this
struct Alpha
{};

/****************************************************** */


struct VR // Creates a constructor for combinatorial simplexes n_pts = n, dim = r
{
public:
    std::vector<int> simplex;
    bool is_initialized; // Marks if the object is initialized
    bool active;         // Marks end of simplexes
    VR() : active(true), is_initialized(false) {};

    void initialize(size_t n_pts, size_t dim)
    {
        n_pts_ = n_pts;
        dim_ = dim;
        simplex.resize(dim); // Resize the vector to the appropriate size
        std::iota(simplex.begin(), simplex.end(), 0);
        is_initialized = true;
    };

    void next_simplex()
    {
        if (simplex[0] == n_pts_ - dim_)
        {
            active = false;
            return; // No more simplices can be generated.
        }

        size_t j = dim_ - 1;

        while (simplex[j] - j == n_pts_ - dim_)
            j--;

        simplex[j]++;

        for (size_t i = j + 1; i < dim_; i++)
            simplex[i] = simplex[j] + i - j;
    };

private:
    size_t n_pts_;
    size_t dim_;
};

template <typename ComplexType, typename DistMatType>
class CritCells : public DistMatType
{
public:
    CritCells(const std::string &fileName);               // FileName to read InputData from
    CritCells(std::vector<std::vector<double>> &distMat); // Input normal Distance matrix
    CritCells(Eigen::SparseMatrix<double> &distMat);
    void run_Compute(int maxDim, int batchsize = 0);
    


    int findRoot(std::vector<int>& parent_idx, int x);
    void setUnion(std::vector<int>& parent_idx, int x, int y);
    std::vector<int> getMSTEdgeIndices(std::vector<std::vector<int>>& sorted_edges);

    //rewrite
    std::vector<std::vector<int64_t>> getBinomialTable(const size_t vtnum, const size_t maxdim);

    int64_t getBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<size_t>& simplex_vt, const size_t shift = 0);

    int64_t getEdgeBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, size_t j, size_t i);

    size_t getSimplexMaxVertex(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t bindex, size_t vtnum, const size_t dim);

    std::vector<size_t> getSimplexVertices(const std::vector<std::vector<int64_t>>& binomial_table, int64_t bindex, size_t vtnum, const size_t dim);

    void sortSimplexByWeightThenIndex(std::vector<std::pair<int64_t, double>>& simplex_list);
    
    std::vector<std::pair<int64_t, double>> getSortedEdge(const std::vector<std::vector<int64_t>>& binomial_table, const double maxeps);

    size_t mstFindRoot(std::vector<size_t>& parent_idx, size_t x);

    void mstSetUnion(std::vector<size_t>& parent_idx, size_t x, size_t y);

    robin_hood::unordered_map<int64_t, size_t> getActiveEdgeIndexHashTable(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<std::pair<int64_t, double>>& sorted_edge);

    // robin_hood::unordered_map<int64_t, size_t> getActiveFacetIndexHashTable(const std::vector<std::pair<int64_t, double>>& facet_list, const std::unordered_set<size_t>& active_facet_index_set);

    std::vector<std::pair<int64_t, double>> getSortedCofacetList(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum);

    std::vector<int64_t> getFacetBinomialIndex(const std::vector<std::vector<int64_t>>& binomial_table, const int64_t binomindex, const size_t dim);

    void buildInterface(Bi_Graph_Match& bi_graph, const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<std::pair<int64_t, double>>& cofacet_list, const size_t dim, const robin_hood::unordered_map<int64_t, size_t>& active_facet_index);
    
    void runMorseTest(size_t maxdim, double maxeps, int threadnumber);

    std::vector<std::pair<int64_t, double>> getSortedDimCells(const std::vector<std::vector<int64_t>>& binomial_table, std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t>& vertex_handle_index, 
                                                              CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>& delaunay_d, const size_t dim, double maxeps);

    double getAlphaSimplexWeight(const std::vector<size_t>& alpha_simplex);

    void runAlphaTest(const std::string &fileName, double maxeps, int threadnumber);



    std::vector<std::unordered_set<size_t>> getVirtualVertexIndices(const std::vector<std::vector<int64_t>>& binomial_table, const size_t maxdim, double& initeps, const int threadnumber);

    std::vector<std::unordered_set<size_t>> mergeIntersectingVirtualVertexIndices(std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);

    void updateBinomialTable(std::vector<std::vector<int64_t>>& binomial_table, const size_t vvtnum, const size_t maxdim);

    std::vector<size_t> getActiveVertex(const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);

    double computeVirtualDistance(const size_t i, const size_t j, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);

    robin_hood::unordered_map<uint64_t, double> getVirtualDistanceHashTable(const std::vector<size_t>& active_vertex, const std::vector<std::unordered_set<size_t>>& virtual_vertex_indices);

    double getVirtualDistance(size_t i, size_t j, const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table);

    std::vector<std::pair<int64_t, double>> getVirtualSortedEdge(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<size_t>& active_vertex, 
                                                                 const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table, const double maxeps);

    robin_hood::unordered_map<int64_t, size_t> getVirtualActiveEdgeIndexHashTable(const std::vector<std::vector<int64_t>>& binomial_table, 
                                                                                  const std::vector<std::pair<int64_t, double>>& sorted_virtual_edge, const size_t vvtnum);
    
    std::vector<std::pair<int64_t, double>> getVirtualSortedCofacetList(const std::vector<std::vector<int64_t>>& binomial_table, const std::vector<std::pair<int64_t, double>>& sorted_virtual_simplex, 
                                                                        const robin_hood::unordered_map<uint64_t, double>& virtual_distance_hash_table, const std::vector<size_t>& active_vertex, const size_t dim, const double maxeps, const int threadnum);

    void runMorseReductionTest(size_t maxdim, double initeps, double maxeps, int threadnumber);
private:
    std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                             // Direct creation of edgebins to a map
    void binByWeights(std::map<double, std::vector<std::vector<int>>> &weighted_simplicies, std::map<double, std::vector<std::vector<int>>> &bins); // Merged higher dim feature to bins
    std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType &simplex_const, size_t dim, size_t batch_size);                  // Worker is invokation counter
    std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>> &simplexes, size_t dim, bool final);
};