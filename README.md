# iCC
Incremental Critical Cells (iCC)

This repository contains the c++ code for studies with iCC that are being performed in the HPC lab at the University of Cincinnati.

## REQUIREMENTS 

- C++23
  
- CMake
  
- OpenMP

- Eigen3

---
			  
## COMPILING 
  
  	##    mkdir build && cd build
	##    cmake .. -DCMAKE_BUILD_TYPE=Release
	##    make

---

##  RUNNING 

	##    ./CritCells <input_filename> <maxDim> [<batch_size>]

## CUSTOMIZATION

1. **Implementation Parameters:**
   - The implementation requires two template parameters to control the input complex type and the distance matrix implementation to be used.
   - Sample implementations can be found in `criticalCells.hpp` files in the form of structs.

2. **Customization Process:**
   - Customization can be achieved by specializing predefined functions.
   - The `DistMat` must provide a distance function to obtain the distance between two indices.

3. **Custom Complex Handling:**
   - If a custom complex is used, there are two options:
     - Store the current simplex in the `simplex` field and implement the `next_simplex` function for seeking.
     - Provide a custom implementation, such as:
       ```cpp
       std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();                                                                        
       std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType &simplex_const, size_t dim, size_t batch_size);
       ```
   - These functions are crucial for handling custom complex types and ensuring proper functionality.

