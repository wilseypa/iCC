#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <chrono>
#include <numeric>

#include <set>
#include <unordered_set>

#include <cstdint>

#include "robin_hood.h"

#ifdef BUILD_ALPHA_COMPLEX
#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#endif

#include "BipartiteGraph.hpp"
#include "DistanceMatrix.hpp"

/****************************************************** */
// temp dummy struct for alpha complex
struct Alpha
{
};

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

/****************************************************** */

template <typename ComplexType, typename DistMatType>
class CritCells : public DistMatType
{
public:
    CritCells(const std::string& filename); // FileName to read InputData from

    CritCells(std::vector<std::vector<double>> point_coordinates); // primary ctor

    // CritCells(Eigen::SparseMatrix<double>& distMat);

    /*legacy*/
#ifdef BUILD_LEGACY_ICC

    struct ICCLegacy
    {
        CritCells<ComplexType, DistMatType>& parent_cc;

        LegacyRunner(CritCells<ComplexType, DistMatType>& critical_cells) : parent_cc(critical_cells) {}

        void run_Compute(int maxDim, int batchsize = 0);

    private:
        std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                             // Direct creation of edgebins to a map
        void binByWeights(std::map<double, std::vector<std::vector<int>>>& weighted_simplicies, std::map<double, std::vector<std::vector<int>>>& bins); // Merged higher dim feature to bins
        std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType& simplex_const, size_t dim, size_t batch_size);                  // Worker is invokation counter
        std::vector<std::vector<int>> dimMatching(std::vector<std::vector<int>>& simplexes, size_t dim, bool final);
    };

    ICCLegacy icc_legacy_runner{*this};

#endif
    /*legacy end*/

    void buildInterface(BipartiteGraph& bi_graph, const std::vector<std::pair<int64_t, double>>& cofacet_list,
                        const robin_hood::unordered_map<int64_t, size_t>& active_facet_index, const std::vector<std::vector<int64_t>>& binomial_table, const size_t dim);

    void morseVRTest(size_t maxdim, double maxeps, int threadnumber);

    void morseVRPH(size_t maxdim, double maxeps, int threadnumber);

#ifdef BUILD_ALPHA_COMPLEX
    void morseAlphaTest(double maxeps, int threadnumber);
#endif

    void morseQuotientAndExpand(const size_t maxdim, const double initeps, const double maxeps, const int threadnumber);

    void morsePiecewisePH(const size_t maxdim, const std::vector<double>& eps_breaks, const int thread_number, const double pv_cap_scale, const bool verbose = false);

private:
    std::vector<std::vector<double>> point_cloud_;
};

extern template class CritCells<VR, NormalDistMat>;
#ifdef BUILD_ALPHA_COMPLEX
extern template class CritCells<Alpha, NormalDistMat>;
#endif
