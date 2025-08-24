#include "omp.h"

#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"


template class SimplexEnumerator<NormalDistMat>;
template class SimplexEnumerator<SparseDistMat>;

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVREdges(const double maxeps)
{
    //********************use openmp later************************//
    
    std::vector<std::pair<int64_t, double>> sorted_edge;

    size_t npt = dist_mat_.dist_mat.size();

    for (size_t i = 0; i < npt - 1; i++)
    {
        for (size_t j = i + 1; j < npt; j++)
        {
            double weight = dist_mat_.getDistance(i, j);
            if (weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(binomial_table, j, i);
                sorted_edge.push_back(std::make_pair(bindex, weight));
            }
        }
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_edge);

    return sorted_edge;
}

template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex_list, const size_t dim, const double maxeps, const int threadnum)
{
    //dim == simplex dimension == cofacet dimension - 1
    std::vector<std::pair<int64_t, double>> cofacet_list;

    size_t npts = dist_mat_.dist_mat.size();

    std::vector< std::vector< std::pair<size_t, double> > > thread_workspace(threadnum);

    omp_set_num_threads(threadnum);

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_simplex_list.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];

        const int64_t bindex = sorted_simplex_list[i].first;
        const double weight = sorted_simplex_list[i].second;

        std::vector<size_t> simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        const size_t minvt = simplex_vertices.back();    //sorted in descending order
        for (size_t covt = 0; covt < minvt; ++covt)
        {
            double newweight = 0.0;
            for (const auto& vt : simplex_vertices)
            {
                newweight = std::max(newweight, dist_mat_.getDistance(covt, vt));
            }

            double cofacetweight = std::max(newweight, weight);

            if (newweight < maxeps)
            {
                int64_t shiftedbindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vertices, 1);
                int64_t cofacetbindex = shiftedbindex + covt;
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    for (const auto& thread_cofacets : thread_workspace)
    {
        cofacet_list.insert(cofacet_list.end(), thread_cofacets.begin(), thread_cofacets.end());
    }

    SimplexUtility::sortSimplexByWeightThenIndex(cofacet_list);

    return cofacet_list;
}


template <typename DistMatType>
std::vector<std::pair<int64_t, double>> SimplexEnumerator<DistMatType>::getSortedAlphaCells(const std::vector<std::vector<int64_t>>& binomial_table, 
                                                                std::unordered_map<CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>::Vertex_handle, size_t>& vertex_handle_index, 
                                                                CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>& delaunay_d, const size_t dim, double maxeps)
{
    std::vector<std::pair<int64_t, double>> sortd_d_cell;

    if (dim == 0) return sortd_d_cell;

    int maxdim = delaunay_d.maximal_dimension();

    std::vector<size_t> simplex_pt;
    simplex_pt.reserve(dim + 1);

    if (dim == maxdim - 1)
    {
        for (auto facetit = delaunay_d.finite_facets_begin(); facetit != delaunay_d.finite_facets_end(); facetit++)
        {
            auto fullcellit = facetit->full_cell();
            auto covt = facetit->index_of_covertex();

            auto neighborit = fullcellit->neighbor(covt);
            if (neighborit < fullcellit) continue;    //skip to avoid double counting the same facet

            for (size_t i = 0; i <= maxdim; i++)
            {
                if (i == covt) continue;    //skip the co_vertex of the facet
                auto vhandle = fullcellit->vertex(i);    //need to use the full cell iter. facet iter does not have vertex method
                simplex_pt.push_back(vertex_handle_index[vhandle]);
            }
            std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());
            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, simplex_pt, 0);
            double weight = getAlphaSimplexWeight(simplex_pt);
            if (weight < maxeps) sortd_d_cell.emplace_back(bindex, weight);
            simplex_pt.clear();
        }
        SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
        return sortd_d_cell;
    }

    if (dim == maxdim)
    {
        for (auto fullcellit = delaunay_d.finite_full_cells_begin(); fullcellit != delaunay_d.finite_full_cells_end(); fullcellit++)
        {
            for (size_t i = 0; i <= dim; i++)
            {
                auto vhandle = fullcellit->vertex(i);
                simplex_pt.push_back(vertex_handle_index[vhandle]);
            }
            std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());
            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, simplex_pt, 0);
            double weight = getAlphaSimplexWeight(simplex_pt);
            if (weight < maxeps) sortd_d_cell.emplace_back(bindex, weight);
            simplex_pt.clear();
        }
        SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
        return sortd_d_cell;
    }

    //for the rest of the dim
    std::vector<size_t> cell_pt;
    cell_pt.reserve(dim + 1);
    std::unordered_set<int64_t> cell_bindex_lookup;
    for (auto fullcellit = delaunay_d.finite_full_cells_begin(); fullcellit != delaunay_d.finite_full_cells_end(); fullcellit++)
    {
        for (size_t i = 0; i <= maxdim; i++)
        {
            auto vhandle = fullcellit->vertex(i);
            simplex_pt.push_back(vertex_handle_index[vhandle]);
        }
        std::sort(simplex_pt.begin(), simplex_pt.end(), std::greater<size_t>());
        
        //subroutine from lhf getDimEdges(int dim)
        size_t powsize = pow(2, maxdim + 1);
        for (size_t counter = 1; counter < powsize; counter++)
        {
            //count the number of 1 in binary form of counter
            if (__builtin_popcount(counter) != dim + 1) continue;

            //collect the corresponding pt(subset) of full cell
            for (size_t i = 0; i < maxdim + 1; i++)
            {
                if (counter & (1 << i)) cell_pt.push_back(simplex_pt[i]);
            }

            std::sort(cell_pt.begin(), cell_pt.end(), std::greater<size_t>());

            int64_t bindex = SimplexUtility::getBinomialIndex(binomial_table, cell_pt, 0);

            // if (bindex == 0)
            // {
            //     std::cout<<"zero bindex edge = "<<cell_pt[0]<<"  "<<cell_pt[1]<<'\n';
            // }

            //already found
            if (cell_bindex_lookup.find(bindex) != cell_bindex_lookup.end()) 
            {
                cell_pt.clear();
                continue;
            }

            cell_bindex_lookup.insert(bindex);
            double weight = getAlphaSimplexWeight(cell_pt);
            if (weight < maxeps) sortd_d_cell.emplace_back(bindex, weight);
            cell_pt.clear();    //clean up the pt array for d cell
        }
        simplex_pt.clear();    //clean up the pt array for full cell
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sortd_d_cell);
    return sortd_d_cell;
}


template <typename DistMatType>
double SimplexEnumerator<DistMatType>::getAlphaSimplexWeight(const std::vector<size_t>& alpha_simplex)
{
    double weight = 0;
    //simplex pt is sorted in descending order
    for (auto rfirst = alpha_simplex.rbegin(); rfirst != alpha_simplex.rend() - 1; rfirst++)
    {
        for (auto rsecond = rfirst + 1; rsecond != alpha_simplex.rend(); rsecond++)
        {
            weight = std::max(weight, dist_mat_.getDistance(*rfirst, *rsecond));
        }
    }
    return weight;
}



