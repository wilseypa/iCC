# iCC Reviewer README

This guide is intended for paper reviewers who want to build and run the current iCC command-line workflow.

## Requirements

- C++23 compiler
- CMake 3.22 or newer
- OpenMP
- Internet access during the first CMake configuration if CLI11 is not already installed locally

## Build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable is:

```bash
./build/iCC
```

## Input

The command-line driver expects a CSV point-cloud input file. Pass the input path with `--file-name`.

## Command-Line Runs

Run standard `morsePH`:

```bash
./build/iCC --tool ph --file-name <input.csv> --max-dim <maxDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

Run `morsePiecewisePH` with explicit epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-breaks <eps1> <eps2> ... --pv-cap-scale <scale> [-n <thread_count>] [-v]
```

Run `morsePiecewisePH` with automatically generated epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-interval-count <count> [--eps-interval-scale <scale>] --pv-cap-scale <scale> [-n <thread_count>] [-v]
```

For complete command-line help:

```bash
./build/iCC --help
```

## Automatic Epsilon Intervals

With `--eps-interval-count`, the program computes all distinct positive pairwise distances from the distance matrix, sorts them, and chooses the requested number of epsilon interval upper bounds from the sorted ranks.

`--eps-interval-scale` controls how those ranks are selected:

- `1.0` gives linear distance ranks.
- Values greater than `1.0` make earlier intervals contain more distinct distances.
- The default is `1.0`.

`--eps-breaks` and `--eps-interval-count` are mutually exclusive.

## Timing

The executable prints `run time = ...` for both verbose and non-verbose runs. For automatic epsilon-interval runs, the timer starts after the interval breaks are generated.

Use `-v` to print diagnostics, including generated epsilon breaks for automatic piecewise runs.

## Suggested Reviewer Workflow

1. Configure and build in Release mode.
2. Run `./build/iCC --help` to verify the available options.
3. Run either `--tool ph` or `--tool piecewise` with the input files and parameters from the paper.
4. Record the printed runtime and generated output for comparison.
