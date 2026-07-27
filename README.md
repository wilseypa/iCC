# iCC

Incremental Critical Cells (iCC)

This repository contains the C++ code for studies with iCC that are being performed in the HPC lab at the University of Cincinnati.

## Requirements

- C++23 compiler
- CMake 3.22 or newer
- OpenMP
- CLI11, fetched automatically by CMake if it is not already installed

## Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable is written to `build/iCC`.

## Running

Interactive mode:

```bash
./build/iCC
```

Command-line mode for `morsePH`:

```bash
./build/iCC --tool ph --file-name <input.csv> --max-dim <maxDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

Command-line mode for `morsePiecewisePH` with explicit epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-breaks <eps1> <eps2> ... --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v]
```

Command-line mode for `morsePiecewisePH` with automatically generated epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-interval-count <count> [--eps-interval-scale <scale>] --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v]
```

For automatic epsilon intervals, the program sorts the distinct positive pairwise distances from the distance matrix and chooses `<count>` interval upper bounds from those ranks. `--eps-interval-scale` must be at least `1.0`; `1.0` gives linear distance ranks, and larger values make earlier intervals contain more distinct distances. The default scale is `1.0`.

`--eps-breaks` and `--eps-interval-count` are mutually exclusive. The aliases `--eps-break-count`, `--eps-window-count`, `--eps-break-scale`, and `--eps-window-scale` are also accepted.

`--pv-cap-scale` bounds a new PV's diameter by `<scale> * eps_hi`, where `eps_hi` is the cut scale at which the PV is formed. `--pv-min-separation` requires every pair of PVs to have minimum cross distance at least `<delta> * eps_max`; its default is `0`, which disables the separation constraint.

The executable accepts any finite positive PV cap scale for experiments. The quotient and single-coordinate safety arguments assume the cap does not exceed the formation scale (`scale <= 1`); because filtration admission uses a strict upper bound, use `scale < 1` to avoid the equality boundary.

For all CLI options:

```bash
./build/iCC --help
```

## Timing

The program prints `run time = ...` in both verbose and non-verbose command-line modes. When automatic epsilon intervals are used, the timer starts after the interval breaks have been generated.

Use `-v` to show diagnostic output, including the generated epsilon breaks for automatic piecewise runs. For piecewise runs with PVs in dimensions 3 through 7, verbose mode also reports canonical-witness false-facet-identification statistics for the top interface. Each report describes the complex rebuilt at that window's upper scale; it is not filtered to cofacets whose weights lie strictly inside the recording interval.

## Reviewer Guide

For paper-review use, see [README_REVIEWER.md](README_REVIEWER.md). That file focuses on building and running the current command-line workflow.

## Customization

1. **Implementation Parameters**
   - The implementation requires two template parameters to control the input complex type and the distance matrix implementation to be used.
   - Sample implementations can be found in `criticalCells.hpp` in the form of structs.

2. **Customization Process**
   - Customization can be achieved by specializing predefined functions.
   - The `DistMat` must provide a distance function to obtain the distance between two indices.

3. **Custom Complex Handling**
   - If a custom complex is used, there are two options:
     - Store the current simplex in the `simplex` field and implement the `next_simplex` function for seeking.
     - Provide a custom implementation, such as:

```cpp
std::map<double, std::vector<std::vector<int>>> binEdgeSimplexes();
std::map<double, std::vector<std::vector<int>>> dsimplices_batches(ComplexType& simplex_const, size_t dim, size_t batch_size);
```

   - These functions are crucial for handling custom complex types and ensuring proper functionality.
