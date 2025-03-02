#include <algorithm>
#include <vector>
#include <ranges>
#include <numeric>
#include <execution>
#include <cstdlib>
#include "omp.h"

#include "biGraphMorseMatch.hpp"

#include <iostream>

void Bi_Graph_Match::parallelKarpSipserInit(const int threadnum) 
{
    std::vector<int> node_deg(u, 0);
    // std::vector<int> node_deg(v, 0);
    std::vector<int> visit_flag(u + v, 0);
    // int nodecount = 0;

    omp_set_num_threads(threadnum);

    std::transform(std::execution::par, adj_list.begin(), adj_list.begin() + u, node_deg.begin(), [](const auto& u_adj) { return int(u_adj.size()); });

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < u; i++)
    // for (int i = u -1; i >= 0; i--)
    {

        if (node_deg[i] == 1)
        {
            // std::cout<<"init non recur/iter start = "<<i<<'\n';
            pairDegreeOne(i, visit_flag, node_deg);
        }
            
    }

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < u; i++)
    //for (int i = u - 1; i >= 0; i--)
    // for (int i = u; i < (u + v); i++)
    {
        if (visit_flag[i] == 0 && node_deg[i] > 0) {
            pairUnmatched(i, visit_flag);
        }
    }

    // auto unmatched = std::count_if(std::execution::par, match_list.begin() + u, match_list.end(), [](int value) { return value < 0; });
    // std::cout<<"unmatched dim - 1 after init = "<<unmatched<<'\n';

    return;
}

void Bi_Graph_Match::pairDegreeOne(const size_t uidx, std::vector<int>& visit_flag, std::vector<int>& node_deg) 
{
    // if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
    //     return;
    // }

    // std::vector<size_t> dfs_stack;

    // for (const auto& vidx : adj_list[uidx]) {
    //     if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
    //         // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
    //         match_list[uidx] = vidx;
    //         match_list[vidx] = uidx;
    //         //update degree of the neighbor of v
    //         for (const auto& index : adj_list[vidx]) {
    //             //found new node with degree == 1
    //             if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) dfs_stack.push_back(index);
    //         }
    //         break;
    //     }
    // }
    // while(!dfs_stack.empty()) {
    //     size_t top = dfs_stack.back();
    //     dfs_stack.pop_back();

    //     if (__sync_fetch_and_add(&(visit_flag[top]), 1) != 0) continue;

    //     for (const auto& vidx : adj_list[top]) {
    //         if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) 
    //         {
    //             // std::cout<<"subsequent uidx = "<<top<<"  match v = "<<vidx<<'\n';
    //             match_list[top] = vidx;
    //             match_list[vidx] = top;
    //             for (const auto& index : adj_list[vidx]) {
    //                 if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) dfs_stack.push_back(index);
    //             }
    //             break;
    //         }
    //     }
    // }
    // return;
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }

    for (const auto& vidx : adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            match_list[uidx] = vidx;
            match_list[vidx] = uidx;
            //update degree of the neighbor of v
            for (const auto& index : adj_list[vidx]) {
                //found new node with degree == 1
                if (__sync_fetch_and_sub(&(node_deg[index]), 1) == 2) {
                    pairDegreeOne(index, visit_flag, node_deg);
                }
            }
            break;
        }
    }
}

void Bi_Graph_Match::pairUnmatched(const size_t uidx, std::vector<int>& visit_flag) 
{
    if (__sync_fetch_and_add(&(visit_flag[uidx]), 1) != 0) {
        return;
    }


    for (const auto& vidx : adj_list[uidx]) {
        if (__sync_fetch_and_add(&(visit_flag[vidx]), 1) == 0) {
            // std::cout<<"pair "<<uidx<<"  "<<vidx<<'\n';
            match_list[uidx] = vidx;
            match_list[vidx] = uidx;
            return;
        }
    }
    return;
}

// void Bi_Graph_Match::elementaryCollapse(int vidx, int umin, int vmin, std::vector<int>& visit_flag_u, std::vector<int>& visit_flag_v, std::vector<int>& node_deg)
// {
    // if (__sync_fetch_and_add(&(visit_flag_v[vidx - vmin]), 1) != 0) return;

    // std::vector<int> dfs_stack;

    // for (const auto& uidx: adj_list[vidx])
    // {
    //     if (__sync_fetch_and_add(&(visit_flag_u[uidx - umin]), 1) == 0)
    //     {
    //         match_list[uidx] = vidx;
    //         match_list[vidx] = uidx;
    //         //update uidx neighbor deg
    //         for (const auto& vindex: adj_list[uidx])
    //         {
    //             if (__sync_fetch_and_sub(&(node_deg[vindex - vmin]), 1) == 2) dfs_stack.push_back(vindex); 
    //         }
    //     }
    // }

    // while(!dfs_stack.empty())
    // {
    //     int vtop = std::move(dfs_stack.back());
    //     dfs_stack.pop_back();

    //     if (__sync_fetch_and_add(&(visit_flag_v[vtop - vmin]), 1) != 0) continue;
        
    //     for (const auto& uidx: adj_list[vtop])
    //     {
    //         if (__sync_fetch_and_add(&(visit_flag_u[uidx - umin]), 1) == 0)
    //         {
    //             match_list[uidx] = vtop;
    //             match_list[vtop] = uidx;
    //             for(const auto& vindex: adj_list[uidx])
    //             {
    //                 if (__sync_fetch_and_sub(&(node_deg[vindex - vmin]), 1) == 2) dfs_stack.push_back(vindex);
    //             }
    //         }
    //     }
    // }

    // return;
// }

void Bi_Graph_Match::parallelMaxFacetInit(const size_t cofacet_index_min, const size_t cofacet_index_max, const size_t facet_index_min, const size_t facet_index_max, const int threadnum)
{
    //convert index
    size_t umin = cofacet_index_min;
    size_t umax = cofacet_index_max;
    size_t vmin = facet_index_min + u;
    size_t vmax = facet_index_max + u;

    std::vector<int> visit_flag_u(umax - umin, 0);
    std::vector<int> visit_flag_v(vmax - vmin, 0);

    // std::vector<int> node_deg(vmax - vmin, 0);
    // std::transform(std::execution::par, adj_list.begin() + vmin, adj_list.begin() + vmax, node_deg.begin(), [](const auto& u_adj) { return u_adj.size(); });

    omp_set_num_threads(threadnum);

//match max facet
#pragma omp parallel for schedule(dynamic)
    for (size_t i = vmin; i < vmax; i++)
    {
        for (const auto& uidx: adj_list[i])    //uidx in v adj is in ascending order
        {
            auto vit = adj_list[uidx].begin();    //vidx in u adj is in descending order
            if (i == *vit)    //i == max facet of uidx
            {
                //no race here
                visit_flag_u[uidx - umin] += 1;
                visit_flag_v[i - vmin] += 1;
                match_list[uidx] = i;
                match_list[i] = uidx;

                break;
            }
        }
    }

//para min cofacet
#pragma omp parallel for schedule(static)
    for (size_t i = vmin; i < vmax; i++)
    {
        if ((match_list[i] < 0) && !(adj_list[i].empty()))
        {
            //if i's first cofacet is available
            if (__sync_fetch_and_add(&(visit_flag_u[adj_list[i][0] - umin]), 1) == 0)
            {
                match_list[adj_list[i][0]] = i;
                match_list[i] = adj_list[i][0];
            }
        }
    }

    // auto unmatched = std::count_if(std::execution::par, match_list.begin() + u, match_list.end(), [](int value) { return value < 0; });
    // std::cout<<"unmatched dim - 1 after init = "<<unmatched<<'\n';

    return;
}


int Bi_Graph_Match::facetDfsAugPath(const size_t startnode, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<size_t>& aug_path_tid) 
{
    int64_t topindex = -1;
    size_t firstlookahead = u + 1;
    aug_path_tid[++topindex] = startnode;

    while (topindex >= 0) 
    {
        size_t vidx = aug_path_tid[topindex];
        int endflag = 0;
        //look ahead, look for v's unmatched neighbor
        for (const auto& uidx : adj_list[vidx]) 
        {
            if (__sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0) 
            {
                if (match_list[uidx] < 0) 
                {   
                    __sync_fetch_and_add(&(dfs_flag[uidx]), 1);
                    aug_path_tid[++topindex] = uidx;
                    return topindex + 1;    //path length
                }
            }
        }
        //dfs
        for (const auto& uidx : adj_list[vidx]) 
        {
            // if (match_list[uidx] < vidx) continue;

            if (__sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0) 
            {
                if (match_list[uidx] >= 0) 
                {
                    aug_path_tid[++topindex] = uidx;
                    aug_path_tid[++topindex] = match_list[uidx];
                    endflag = 1;
                    break;
                }
            }
        }
        //dfs cannot augment, pop
        if (!endflag) {
            topindex -= 2;
        }
    }
    return topindex + 1;
}

// bool Bi_Graph_Match::isRightSinglePath(int cofacetindex)
// {
//     //return true if cofacet has only one right neighbor
//     //return false if multi or no right neighbor
//     if (match_list[cofacetindex] < 0) return false;

//     int vidx = match_list[cofacetindex];
//     //iter to the second to the last item
//     //facet adj list is in ascending order
//     std::vector<int>::reverse_iterator rit = adj_list[vidx].rbegin() - 1;
//     if (*rit == cofacetindex)
//     {
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }

// int Bi_Graph_Match::facetRightSinglePathDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid)
// {
//     //if not single, only lookahead
//     int topindex = -1;
//     aug_path_tid[++topindex] = facetindex;
//     int endflag = 0;
//     int cycleflag = 0;
//     while (topindex >= 0)
//     {
//         endflag = 0;
//         int vidx = aug_path_tid[topindex];

//         //try look ahead. look for unmatched
//         if (topindex != 0)
//         {
//             for (const auto& uidx: adj_list[vidx])
//             {
//                 //already matched. pass
//                 if (match_list[uidx] >= 0) continue;

//                 cycleflag = 0;
//                 //check acyclicity
//                 for (int i = 0; i < topindex; i += 2)
//                 {
//                     if (std::find(adj_list[uidx].begin(), adj_list[uidx].end(), aug_path_tid[i]) != adj_list[uidx].end())
//                     {
//                         cycleflag = 1;
//                         break;
//                     }
//                 }

//                 if (!cycleflag && __sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//                 // if (uidx > aug_path_tid[topindex - 1] && __sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0)
//                 {
//                     aug_path_tid[++topindex] = uidx;
//                     return topindex + 1;    //aug path length
//                 }
//             }
//         }

//         //try aug
//         for (const auto& uidx: adj_list[vidx])
//         {
//             if (match_list[uidx] < vidx || !isRightSinglePath(uidx)) continue;

//             cycleflag = 0;
//             //check acyclicity
//             for (int i = 0; i < topindex; i += 2)
//             {
//                 if (std::find(adj_list[uidx].begin(), adj_list[uidx].end(), aug_path_tid[i]) != adj_list[uidx].end())
//                 {
//                     cycleflag = 1;
//                     break;
//                 }
//             }

//             if (!cycleflag && __sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//             {
//                 aug_path_tid[++topindex] = uidx;
//                 aug_path_tid[++topindex] = match_list[uidx];
//                 endflag = 1;
//                 break;
//             }
//         }
//         //cannot right augment, pop
//         if (!endflag) topindex -= 2;
//     }
//     return topindex + 1;
// }

// int Bi_Graph_Match::facetRightDFSAugPathWithCheck(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid)
// {
//     int topindex = -1;
//     aug_path_tid[++topindex] = facetindex;
//     int endflag = 0;
//     int cycleflag = 0;
//     while (topindex >= 0)
//     {
//         endflag = 0;
//         int vidx = aug_path_tid[topindex];
        
//         //try look ahead. look for unmatched
//         if (topindex != 0)
//         {
//             //int vidx = aug_path_tid[topindex];
//             for (const auto& uidx: adj_list[vidx])
//             {
//                 //already matched. pass
//                 if (match_list[uidx] >= 0) continue;

//                 cycleflag = 0;
//                 //check acyclicity
//                 for (int i = 0; i < topindex; i += 2)
//                 {
//                     if (std::find(adj_list[uidx].begin(), adj_list[uidx].end(), aug_path_tid[i]) != adj_list[uidx].end())
//                     {
//                         cycleflag = 1;
//                         break;
//                     }
//                 }

//                 if (!cycleflag && __sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//                 // if (uidx > aug_path_tid[topindex - 1] && __sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0)
//                 {
//                     aug_path_tid[++topindex] = uidx;
//                     return topindex + 1;    //aug path length
//                 }
//             }
//         }
//         //try aug
//         for (const auto& uidx: adj_list[vidx])
//         {   
//             if (match_list[uidx] < vidx) continue;

//             cycleflag = 0;
//             //check acyclicity
//             for (int i = 0; i < topindex; i += 2)
//             {
//                 if (std::find(adj_list[uidx].begin(), adj_list[uidx].end(), aug_path_tid[i]) != adj_list[uidx].end())
//                 {
//                     cycleflag = 1;
//                     break;
//                 }
//             }

//             if (!cycleflag && __sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//             {
//                 aug_path_tid[++topindex] = uidx;
//                 aug_path_tid[++topindex] = match_list[uidx];
//                 endflag = 1;
//                 break;
//             }
//         }
//         //cannot right augment, pop
//         if (!endflag) topindex -= 2;
//     }
//     return topindex + 1;
// }

// int Bi_Graph_Match::facetLeftDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid)
// {
//     int topindex = -1;
//     aug_path_tid[++topindex] = facetindex;
//     while (topindex >= 0)
//     {
//         int vidx = aug_path_tid[topindex];
//         int endflag = 0;
//         //if augmented, try look ahead to the left
//         if (topindex != 0)
//         {
//             int vidx = aug_path_tid[topindex];
//             for (const auto& uidx: adj_list[vidx])
//             {
//                 if (uidx < aug_path_tid[topindex - 1] && __sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0)
//                 {
//                     if (match_list[uidx] < 0)
//                     {
//                         __sync_fetch_and_add(&(dfs_flag[uidx]), 1);
//                         aug_path_tid[++topindex] = uidx;
//                         return topindex + 1;    //aug path length
//                     }
//                 }
//             }
//         }
//         //try left aug
//         for (const auto& uidx: adj_list[vidx])
//         {
//             if (match_list[uidx] >= 0 && match_list[uidx] < vidx)
//             {
//                 if (__sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//                 {
//                     aug_path_tid[++topindex] = uidx;
//                     aug_path_tid[++topindex] = match_list[uidx];
//                     endflag = 1;
//                     break;
//                 }
//             }
//         }
//         //cannot left augment, pop
//         if (!endflag) topindex -= 2;
//     }
//     return topindex + 1;
// }


// int Bi_Graph_Match::facetRightDFSAugPath(int facetindex, std::vector<int>& dfs_flag, std::vector<int>& look_ahead_flag, std::vector<int>& aug_path_tid)
// {
//     int topindex = -1;
//     aug_path_tid[++topindex] = facetindex;
//     while (topindex >= 0)
//     {
//         int vidx = aug_path_tid[topindex];
//         int endflag = 0;
//         //if augmented, try look ahead to the right
//         if (topindex != 0)
//         {
//             // int vidx = aug_path_tid[topindex];
//             for (const auto& uidx: adj_list[vidx])
//             {
//                 if (__sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0)
//                 // if (uidx > aug_path_tid[topindex - 1] && __sync_fetch_and_add(&(look_ahead_flag[uidx]), 1) == 0)
//                 {
//                     if (match_list[uidx] < 0)
//                     {
//                         __sync_fetch_and_add(&(dfs_flag[uidx]), 1);
//                         aug_path_tid[++topindex] = uidx;
//                         return topindex + 1;    //aug path length
//                     }
//                 }
//             }
//         }
//         //try right aug
//         for (const auto& uidx: adj_list[vidx])
//         {
//             // if (match_list[uidx] > vidx)
//             if (match_list[uidx] >= 0)
//             {
//                 if (__sync_fetch_and_add(&(dfs_flag[uidx]), 1) == 0)
//                 {
//                     aug_path_tid[++topindex] = uidx;
//                     aug_path_tid[++topindex] = match_list[uidx];
//                     endflag = 1;
//                     break;
//                 }
//             }
//         }
//         //cannot right augment, pop
//         if (!endflag) topindex -= 2;
//     }
//     return topindex + 1;
// }


void Bi_Graph_Match::add2SingleOrRemove(const size_t index, std::vector<uint64_t>& removed_flag, const size_t round)
{
    int64_t flagabs = abs(removed_flag[index]);

    //if index is not visited in this round
    if (flagabs < round)
    {
        removed_flag[index] = round;
        return;
    }

    //if index has been visited once in this round
    if (removed_flag[index] == round)
    {
        removed_flag[index] = -round;
        return;
    }

    //if index is removed in current round. removed_flag[index] == -round
    if (flagabs == round) return;
}

int Bi_Graph_Match::serialCofacetLeftDFSAugPath(const size_t cofacetindex, std::vector<uint8_t>& cofacet_dfs_flag, std::vector<size_t>& aug_path, std::vector<uint64_t>& removed_flag, const size_t round)
{
    std::vector<size_t> dfs_stack;
    dfs_stack.reserve(u);
    dfs_stack.push_back(cofacetindex);

    size_t maxindex = u + v;

    cofacet_dfs_flag[cofacetindex] += 1;

    int64_t topindex;

    // std::cout<<"ustart = "<<cofacetindex<<"  "<<dfs_stack.back()<<'\n';

    while (!dfs_stack.empty())
    {
        topindex = dfs_stack.back();
        dfs_stack.pop_back();

        // std::cout<<"stack topindex = "<<topindex<<'\n';

        for (auto& vidx: adj_list[topindex])
        {

            if (match_list[vidx] < 0) 
            {
                add2SingleOrRemove(vidx, removed_flag, round);
                continue;
            }

            size_t uidx = match_list[vidx];
            if (uidx != topindex)
            {
                dfs_stack.push_back(uidx);
                add2SingleOrRemove(uidx, removed_flag, round);
            }
        }
    }

    int64_t maxfacetindex = -1;
    for (int64_t i = maxindex - 1; i >= u; i--)
    {
        if (removed_flag[i] == round)
        {
            maxfacetindex = i;
            break;
        }
    }
    if (maxfacetindex < 0) return -1;    //if no avaiable facet

    int endflag = 0;
    //pick the largest availble facet. do backward dfs to find the aug path

    topindex = -1;
    aug_path[++topindex] = maxfacetindex;

    std::cout<<"cofacet idx = "<<cofacetindex<<"  max facet = "<<maxfacetindex<<"  rel max facet = "<<maxfacetindex - u<<"  adj list = ";
    for (auto t: adj_list[cofacetindex]) std::cout<<t<<"  ";
    std::cout<<'\n';
        
    while (!endflag)
    {
        endflag = 1;
        size_t vidx = aug_path[topindex];
        for (auto& uidx: adj_list[vidx])
        {
            if (uidx == match_list[vidx]) continue;

            //end of backward dfs
            if (uidx == cofacetindex)
            {
                aug_path[++topindex] = uidx;
                return topindex + 1;
            }

            //found eligible backward path
            // if (single_cofacet.contains(uidx) && cofacet_dfs_flag[uidx] < 1)
            if (removed_flag[uidx] == round)
            {   
                cofacet_dfs_flag[uidx] += 1;
                aug_path[++topindex] = uidx;
                aug_path[++topindex] = match_list[uidx];
                endflag = 0;
                break;
            }
        }

        //cannot go backward to the cofacetindex with this topindex
        if (endflag)
        {
            for (int64_t j = topindex - 1; j > 0; j -= 2) 
            {
                cofacet_dfs_flag[aug_path[j]] -= 1;
            }
        }
    } 
    
    return -1;
}

// std::set<int> Bi_Graph_Match::criticalFacetBackwardSearch(int facetindex)
// {
//     //facet index is the index in simplex/facet bin
//     //facetindex + u = index in bi graph
//     std::vector<int> dfs_stack;
//     dfs_stack.push_back(facetindex + u);

//     std::set<int> single_facet;
//     std::set<int> removed_index;

//     single_facet.insert(facetindex);

//     while(!dfs_stack.empty())
//     {
//         int top = dfs_stack.back();
//         dfs_stack.pop_back();

//         for (auto uidx: adj_list[top])
//         { 
//             if (match_list[uidx] < 0) continue;

//             int vidx = match_list[uidx];

//             if (vidx != top)
//             {
//                 dfs_stack.push_back(vidx);
//                 add2SingleOrRemove(vidx - u, single_facet, removed_index);
//             }
//         }
//     }

//     return single_facet;
// }


// std::vector<std::set<int>> Bi_Graph_Match::getBackwardSingleFacetIndex(std::vector<int> critical_facet_index)
// {
//     std::vector<std::set<int>> backward_single_facet_index;

//     for (auto ci : critical_facet_index)
//     {
//         std::set<int> ci_backward_index = criticalFacetBackwardSearch(ci);

//         //remove duplicates with previous ci's
//         for (auto& facet_index_set: backward_single_facet_index)
//         {
//             for (auto i: facet_index_set)
//             {
//                 auto it = ci_backward_index.find(i);
//                 if (it != ci_backward_index.end()) ci_backward_index.erase(it);
//             }
//         }

//         backward_single_facet_index.push_back(ci_backward_index);

        
//     }

//     return backward_single_facet_index;
// }


void Bi_Graph_Match::serialCofacetDFSMatch()
{
    std::vector<size_t> unmatched_u_init(u, 0);
    std::vector<size_t> unmatched_u_final(u, 0);
    std::vector<uint8_t> cofacet_dfs_flag(u, 0);

    size_t initialunmatched = 0;
    size_t finalunmatched = 0;

    for(size_t i = 0; i < u; i++)
    {
        if (match_list[i] < 0 && adj_list[i].size() > 0) unmatched_u_init[initialunmatched++] = i;
    }

    std::vector<size_t> aug_path(u, 0);

    std::vector<uint64_t> removed_flag(u + v, 0);

    size_t round = 0;

    while (true)
    {
        //do not reset dfs flag

        // std::cout<<"dfs round "<<'\n';
        std::fill(std::execution::par, removed_flag.begin(), removed_flag.end(), 0);

        for (size_t i = 0; i < initialunmatched; i++)
        {
            round = i + 1;

            size_t ustart = unmatched_u_init[i];

            // std::cout<<"ustart = "<<ustart<<"  match = "<<match_list[ustart]<<'\n';
            
            int augpathlen = serialCofacetLeftDFSAugPath(ustart, cofacet_dfs_flag, aug_path, removed_flag, round);

            // std::cout<<"returned from ith init unmatched = "<<i<<'\n';

            for (int j = 0; j < augpathlen; j += 2) 
            {
                match_list[aug_path[j]] = aug_path[j + 1];
                match_list[aug_path[j + 1]] = aug_path[j];
            }

            if (augpathlen <= 0) unmatched_u_final[finalunmatched++] = ustart;  
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched)) break;

        std::swap(unmatched_u_init, unmatched_u_final);

        initialunmatched = finalunmatched;

        finalunmatched = 0;
    }

    return;
}


// int Bi_Graph_Match::findRoot() {
//     auto root_lambda = [&](const auto i) {return (state_flag[i] == 0);};
//     std::vector<int> range(u);
//     std::iota(range.begin(), range.end(), 0);
//     auto it = std::find_if(range.begin(), range.end(), root_lambda);
//     if (it != range.end()) {
//         return *it;
//     } else {
//         return -1;
//     }
// }

// int Bi_Graph_Match::getParent(std::vector<int>& parent_workspace, int uidx) {
//     parent_workspace.clear();
//     int vidx = match_list[uidx];
//     //uidx is a source node has no parent (umatched)
//     if (vidx == -1) return 0;
//      //i is d simp
//     for (auto& i: adj_list[vidx]) {
//         int imate = match_list[i];
//         if (imate != -1 && imate != vidx) parent_workspace.push_back(i);
//     }
//     return parent_workspace.size();
// }

size_t Bi_Graph_Match::getChild(std::vector<size_t>& child_workspace, const size_t uidx) 
{
    child_workspace.clear();
    //i is d-1 simp
    for (auto& i: adj_list[uidx]) {
        int64_t imate = match_list[i];
        if (imate != -1 && imate != uidx) child_workspace.push_back(imate); 
    }

    return child_workspace.size();
}

// void Bi_Graph_Match::getAncestor(std::vector<int>& ancestor_workspace, int rootnum, int uidx) {
//     //use rootnum to confine the ancestor search in current "bfs component"
//     ancestor_workspace.clear();

//     std::deque<int> bfs_queue;
//     std::vector<int> parent_workspace;
//     parent_workspace.reserve(udegree);

//     //temp change uidx's rootflag to avoid cycle infi loop
//     int rootflag = state_flag[uidx];
//     state_flag[uidx] = -1;
    
//     getParent(parent_workspace, uidx);
//     for (auto& i: parent_workspace) {
//         if(state_flag[i] == rootnum) bfs_queue.push_back(i);
//     }

//     //cycle can be found during this process but need ancestors for look up
//     while (!bfs_queue.empty()) {
//         int front = std::move(bfs_queue.front());
//         bfs_queue.pop_front();
//         if (std::find(ancestor_workspace.begin(), ancestor_workspace.end(), front) == ancestor_workspace.end()) ancestor_workspace.push_back(front);
        
//         //parent of queue front
//         getParent(parent_workspace, front);
//         for(auto& i: parent_workspace) {
//             if (state_flag[i] == rootnum && std::find(bfs_queue.begin(), bfs_queue.end(), i) == bfs_queue.end()) {
//                 bfs_queue.push_back(i);
//             }
//         }
//     }

//     //recover rootflag
//     state_flag[uidx] = rootflag;

//     return;
// }

// bool Bi_Graph_Match::isBackwardAcyclic(std::vector<int>& ancestor_d_simp, std::vector<int>& u_child) {
    //par
//     for (auto& i: u_child) {
//         if (state_flag[i] == 0) continue;
//         if (std::find(ancestor_d_simp.begin(), ancestor_d_simp.end(), i) != ancestor_d_simp.end()) {
//             std::cout<<"found cycle at = "<<i<<"   ancestor simp = ";
//             for(auto i: ancestor_d_simp) std::cout<<i<<"  ";
//             std::cout<<"  child simp = ";
//             for(auto j: u_child) std::cout<<j<<"  ";
//             std::cout<<'\n';
//             return false;
//         }
//     }

//     return true;
// }

// int Bi_Graph_Match::lookAheadDFS(std::deque<int>& graph_bfs_queue, std::vector<int>& ancestor_d_simp, std::vector<int>& child_d_simp, int rootnum, int uidx) {
//     if (child_d_simp.size() == 0) return 0;

//     std::vector<int> lookahead_vec;

//     //par ch workspace do not overlap. can be reused
//     std::vector<int> parent_workspace;
//     parent_workspace.reserve(udegree);

//     std::vector<int> child_workspace;
//     child_workspace.reserve(udegree);

//     for (auto& i: child_d_simp) {
//         int parentnum = getParent(parent_workspace, i);
//         //if uidx is the only parent of i (in the whole graph). use look ahead shortcut
//         if (parentnum == 1) {
//             lookahead_vec.push_back(i);
//         } else if (state_flag[i] == 0 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
//             //if i has more than 1 par and has not been searched yet (or not in current bfs component). push to global bfs queue
//             graph_bfs_queue.push_back(i);
//         }
//     }

//     //start look ahead op on children of uidx

//     ancestor_d_simp.push_back(uidx);
//     int reverted = 0;
//     bool flag;    //working var
//     //the nodes appeared in lookahead stack at the same time are independent 
//     while (!lookahead_vec.empty()) {

//         flag = false;

//         int top = std::move(lookahead_vec.back());
//         lookahead_vec.pop_back();
//         state_flag[top] = rootnum;

//         //child_workspace contains the children of top
//         int childnum = getChild(child_workspace, top);
//         if (childnum == 0) continue;

//         //child in unsearched part of current bfs component
//         auto iter = std::remove_if(child_workspace.begin(), child_workspace.end(), [&](const auto& i) {return (state_flag[i] != 0 && state_flag[i] != rootnum);});
//         child_workspace.erase(iter, child_workspace.end());

//         //if top is not acyclic. the lookahead of top is over
//         //push its unsearched children or reencountered children to bfs queue
//         if (!isBackwardAcyclic(ancestor_d_simp, child_workspace)) {
            
//             int temp = match_list[top];

//             match_list[top] = -1;
//             match_list[temp] = -1;
//             reverted += 1;

//         } else {
//         //check if the child of top can be pushed to look ahead stack
//             for (auto& i: child_workspace) {
//                 if (getParent(parent_workspace, i) == 1) {
//                     //if top is the only parent of i. keep look ahead op
//                     lookahead_vec.push_back(i);
//                     //raise flag for adding top to ancestor
//                     flag = true;
//                 }else if (state_flag[i] == 0 && std::find(graph_bfs_queue.begin(), graph_bfs_queue.end(), i) == graph_bfs_queue.end()) {
//                     graph_bfs_queue.push_back(i);
//                 } 
//             }
//         }
//         //if any of top's child is pushed to look ahead stack, update ancestor
//         if (flag) ancestor_d_simp.push_back(top);
//         //not all nodes in the stack have ances-des relationship with top
//         //but this does not bring extra cycle rm. simplfied implementation
//     }
//     return reverted;
//  }

//  int Bi_Graph_Match::serialBFS(std::vector<int>& ancestor_d_simp, int rootnum, int root) {
//     int reverted = 0;

//     std::deque<int> bfs_queue;
//     bfs_queue.push_back(root);

//     std::vector<int> front_child;
//     front_child.reserve(udegree);

//     while (!bfs_queue.empty()) {
//         int front = std::move(bfs_queue.front());
//         bfs_queue.pop_front();
//         state_flag[front] = rootnum;
        
//         int childnum = getChild(front_child, front);
//         if (childnum == 0) continue;

//         //child in unsearched part of current bfs component
//         auto iter = std::remove_if(front_child.begin(), front_child.end(), [&](const auto& i) {return (state_flag[i] != 0 && state_flag[i] != rootnum);});
//         front_child.erase(iter, front_child.end());

//         getAncestor(ancestor_d_simp, rootnum, front);

//         if (!isBackwardAcyclic(ancestor_d_simp, front_child)) {
//             int temp = match_list[front];

//             match_list[front] = -1;
//             match_list[temp] = -1;
//             reverted += 1;
//         } else {
//             //look ahead and push descendant to bfs queue
//             reverted += lookAheadDFS(bfs_queue, ancestor_d_simp, front_child, rootnum, front);
//         }
//     }
    
//     return reverted;
// }


void Bi_Graph_Match::parallelFacetDFSMatch(const int threadnum) 
{
    std::vector<size_t> unmatched_v_init(v, 0);
    std::vector<size_t> unmatched_v_final(v, 0);
    std::vector<int> dfs_flag(u + v, 0);
    std::vector<int> look_ahead_flag(u, 0);

    omp_set_num_threads(threadnum);

    size_t initialunmatched = 0;
    size_t finalunmatched = 0;
    //find unmatched left nodes from initialized matching
#pragma omp parallel
    {
        size_t ct = 0;
        std::vector<size_t> thread_buff;    //overflow? limit its size?
#pragma omp for
        for (size_t i = u; i < u + v; i++)
        // for (size_t i = u + v - 1; i >= u; i--)
        {
            if (match_list[i] < 0 && adj_list[i].size() > 0) {
                thread_buff.push_back(i);
                ct += 1;
            }
        }
        if (ct > 0) {
            int offset = __sync_fetch_and_add(&initialunmatched, ct);
            std::copy_n(thread_buff.begin(), ct, unmatched_v_init.begin() + offset);
        }
    }

    //storing aug path one path per thread
    std::vector<std::vector<size_t>> aug_path(threadnum, std::vector<size_t>(u + v, 0));

    while (true) {

        // std::cout<<"aug match round "<<'\n';

        //shared among threads
        finalunmatched = 0;

        std::fill(std::execution::par, dfs_flag.begin(), dfs_flag.end(), 0);

#pragma omp parallel for schedule(dynamic)
        // for (int i = 0; i < initialunmatched; i++)
        for (int i = initialunmatched - 1; i >= 0; i--)
        {
            auto& aug_path_tid = aug_path[omp_get_thread_num()];

            int vstart = unmatched_v_init[i];
            int augpathlen = facetDfsAugPath(vstart, dfs_flag, look_ahead_flag, aug_path_tid);

            //augmentation
            for (int j = 0; j < augpathlen; j += 2) {
                match_list[aug_path_tid[j]] = aug_path_tid[j + 1];
                match_list[aug_path_tid[j + 1]] = aug_path_tid[j];
            }
            //store the unmatched node, need atomic op on the shared count var
            if (augpathlen <= 0)
                unmatched_v_final[__sync_fetch_and_add(&finalunmatched, 1)] = vstart;
        }

        if ((finalunmatched == 0) || (initialunmatched == finalunmatched)) {
            break;
        }
        //try more aug paths, let final_unmatched be the init_unmatched
        std::swap(unmatched_v_init, unmatched_v_final);

        initialunmatched = finalunmatched;
    }

    //non-isolated unmatched
    return;   
}

// void Bi_Graph_Match::parallelDirectionalFacetDFSMatch()
// {
//     std::vector<int> unmatched_v_init(v, -1);
//     std::vector<int> unmatched_v_final(v, 0);
//     std::vector<int> dfs_flag(u + v, 0);

//     //size?????
//     std::vector<int> look_ahead_flag(u + v, 0);  

//     omp_set_num_threads(maxthreadnum);

//     int initialunmatched = 0;
//     int finalunmatched = 0;
//     //find unmatched left nodes from initialized matching
// #pragma omp parallel
//     {
//         int ct = 0;
//         std::vector<int> thread_buff;    //overflow? limit its size?
// #pragma omp for
//         for (int i = u; i < u + v; i++) 
//         {
//             if (match_list[i] < 0 && adj_list[i].size() > 0) {
//                 thread_buff.push_back(i);
//                 ct += 1;
//             }
//         }
//         if (ct > 0) {
//             int offset = __sync_fetch_and_add(&initialunmatched, ct);
//             std::copy_n(thread_buff.begin(), ct, unmatched_v_init.begin() + offset);
//         }
//     }


//     //storing aug path one path per thread
//     std::vector<std::vector<int>> aug_path(maxthreadnum, std::vector<int>(u + v, 0));

//     while (true) {
//         //shared among threads
//         finalunmatched = 0;

//         std::fill(std::execution::par, dfs_flag.begin(), dfs_flag.end(), 0);

 
// #pragma omp parallel for schedule(dynamic)
//         for (int i = 0; i < initialunmatched; i++)
//         {
//             auto& aug_path_tid = aug_path[omp_get_thread_num()];

//             int vstart = unmatched_v_init[i];

//             // int augpathlen = facetRightDFSAugPath(vstart, dfs_flag, look_ahead_flag, aug_path_tid);
//             int augpathlen = facetRightDFSAugPathWithCheck(vstart, dfs_flag, look_ahead_flag, aug_path_tid);

//             //augmentation
//             for (int j = 0; j < augpathlen; j += 2) {
//                 match_list[aug_path_tid[j]] = aug_path_tid[j + 1];
//                 match_list[aug_path_tid[j + 1]] = aug_path_tid[j];
//             }
//             //store the unmatched node, need atomic op on the shared count var
//             if (augpathlen <= 0)
//                 unmatched_v_final[__sync_fetch_and_add(&finalunmatched, 1)] = vstart;
//         }
        
//         if ((finalunmatched == 0) || (initialunmatched == finalunmatched)) {
//             break;
//         }
//         //try aug again, let final_unmatched be the init_unmatched
//         std::swap(unmatched_v_init, unmatched_v_final);
//         initialunmatched = finalunmatched;
//     }

//     //non-isolated unmatched
//     return;  
// }


// // int Bi_Graph_Match::serialCycleRemoval() {
//     //assume u is the d simplex of d interface
//     state_flag.resize(u, 0);

//     int reverted = 0;

//     int root = findRoot();
//     int rootnum = 0;

//     //ancestor workspace. avoid multiple reallocation
//     std::vector<int> ancestor_workspace;
//     ancestor_workspace.reserve(u);

//     while (root != -1) {

//         // std::cout << "current root = " << root << '\n';

//         rootnum += 1;

//         reverted += serialBFS(ancestor_workspace, rootnum, root);

//         root = findRoot();
//     }
//     return reverted;
// }

int Bi_Graph_Match::dfsCycleRemoval()
{
    int reverted = 0;

    //0: not visited, 1: visiting, 2: visited
    std::vector<uint8_t> state_flag(u, 0);

    std::vector<size_t> dfs_stack;
    dfs_stack.reserve(v);

    std::vector<size_t> child_workspace;
    child_workspace.reserve(udegree + 1);


    for (size_t i = 0; i < u; i++)
    {

        if (state_flag[i] != 0) continue;


        dfs_stack.push_back(i);

        while (!dfs_stack.empty())
        {
            size_t top = dfs_stack.back();
            dfs_stack.pop_back();

            if (state_flag[top] == 0)
            {
                state_flag[top] = 1;
                dfs_stack.push_back(top);

                getChild(child_workspace, top);

                for (size_t child: child_workspace)
                {

                    if (state_flag[child] == 0) 
                    {
                        dfs_stack.push_back(child);
                    } else if (state_flag[child] == 1)
                    {
                        //found back edge
                        int64_t temp = match_list[child];

                        // if (temp < 0)
                        // {
                        //     std::cout<<"found negative dfs child match. top = "<<top<<"    adj list = ";
                        //     for (auto k : adj_list[top]) std::cout<<k<<"("<<k-u<<")"<<"  ";
                        //     std::cout<<'\n';
                        //     continue;
                        // }
                        match_list[child] = -1;
                        match_list[temp] = -1;
                        reverted += 1;
                    }
                }
            } else if (state_flag[top] == 1) state_flag[top] = 2;
        }
        
    }

    return reverted;
}


void Bi_Graph_Match::addEdge(int leftnode, int rightnode) {
    adj_list[leftnode].push_back(u + rightnode);
    adj_list[u + rightnode].push_back(leftnode);
    return;
}

Bi_Graph_Match::Bi_Graph_Match(size_t leftnum, size_t rightnum, size_t leftdimension) : u(leftnum), v(rightnum), udegree(leftdimension + 1) 
{
    adj_list.resize(u + v);
    match_list.resize(u + v, -1);
}

void Bi_Graph_Match::updateDimension(size_t newleftnum, size_t newrightnum) {
    u = newleftnum;
    v = newrightnum;
    udegree += 1;
    adj_list.clear();
    adj_list.resize(u + v);
    match_list.clear();
    match_list.resize(u + v, -1);
}



// void Bi_Graph_Match::buildInterface(const std::vector<std::vector<int>>& cofacet_bin, int cofacet_index_min, int cofacet_index_max, const std::vector<std::vector<int>>& simplex_bin, int simplex_index_max, const std::vector<int>& active_index) 
// {    
//     std::fill(match_list.begin(), match_list.end(), -1);

//     for(int i = cofacet_index_min; i < cofacet_bin.size(); i++)
//     {
//         for (int j = active_index.size() - 1; j >= 0; j--)
//         {
//             //if (active_index[j] >= simplex_index_max) continue;

//             if (std::includes(cofacet_bin[i].begin(), cofacet_bin[i].end(), simplex_bin[active_index[j]].begin(), simplex_bin[active_index[j]].end())) 
//             {   
//                 this->addEdge(i, active_index[j]);
//             }
            
//         }
//     }

// }

std::vector<size_t> Bi_Graph_Match::getActiveIndex() 
{
    std::vector<size_t> active_index;
    active_index.reserve(u);
    for (size_t i = 0; i < u; i++) 
    {
        if (match_list[i] < 0) active_index.push_back(i);
    }
    return active_index;
}

std::vector<size_t> Bi_Graph_Match::getCriticalIndex(const std::vector<size_t>& dim_active_index, const size_t simplex_index_max) 
{
    std::vector<size_t> critical_index;

    for (size_t i = u; i < (u + simplex_index_max); i++) 
    {
        //for kth simplex, its index in graph i == k + u
        if (match_list[i] < 0)
        {
            if (std::find(dim_active_index.begin(), dim_active_index.end(), (i - u)) != dim_active_index.end()) 
            {
                critical_index.push_back(i - u);
            }
        }            
    }
    return critical_index;
}


void Bi_Graph_Match::checkSimplex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_simplex)
{
    for (auto& simplex: target_simplex)
    {
        auto it = std::find(simplex_bin.begin(), simplex_bin.end(), simplex);
        int index = std::distance(simplex_bin.begin(), it);
        std::cout<<"target simplex = ";
        for(auto pt: simplex) std::cout<<pt<<" ";
        std::cout<<"  index = "<<index;
        std::cout<<"  adj index = ";
        for(auto idx: adj_list[u + index]) std::cout<<idx<<" ";
        std::cout<<"  match = "<<match_list[u + index]<<'\n';
    }
}

void Bi_Graph_Match::checkCofacet(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<std::vector<int>>& target_cofacet)
{
    for (auto& cofacet: target_cofacet)
    {
        auto it = std::find(cofacet_bin.begin(), cofacet_bin.end(), cofacet);
        int index = std::distance(cofacet_bin.begin(), it);
        std::cout<<"target cofacet = ";
        for(auto pt: cofacet) std::cout<<pt<<" ";
        std::cout<<"  index = "<<index;
        std::cout<<"  adj index = ";
        for(auto idx: adj_list[index]) std::cout<<idx - u<<" ";
        std::cout<<"  match = "<<match_list[index] - u<<'\n';
    }
}

void Bi_Graph_Match::checkSimplexByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_simplex_index)
{
    for (auto& i: target_simplex_index)
    {
        auto simplex = simplex_bin[i];
        std::cout<<"target simplex = ";
        for(auto pt: simplex) std::cout<<pt<<" ";
        std::cout<<"  index = "<<i;
        std::cout<<"  adj index = ";
        for(auto idx: adj_list[u + i]) std::cout<<idx<<" ";
        std::cout<<"  match = "<<match_list[u + i]<<'\n';
    }
}

void Bi_Graph_Match::checkCofacetByIndex(std::vector<std::vector<int>>& cofacet_bin, std::vector<std::vector<int>>& simplex_bin, std::vector<int>& target_cofacet_index)
{
    for (auto& i: target_cofacet_index)
    {
        auto cofacet = cofacet_bin[i];
        std::cout<<"target cofacet = ";
        for(auto pt: cofacet) std::cout<<pt<<" ";
        std::cout<<"  index = "<<i;
        std::cout<<"  adj index = ";
        for(auto idx: adj_list[i]) std::cout<<idx - u<<" ";
        std::cout<<"  match = "<<match_list[i] - u<<'\n';
    }
}