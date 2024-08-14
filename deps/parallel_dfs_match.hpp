#pragma once

#include <vector>
#include "bi_graph.hpp"

int parallelDFSMatch(Bi_Graph*, int);
int dfsAugPath(Bi_Graph*, int, std::vector<int>&, std::vector<int>&, std::vector<int>&);