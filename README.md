# iCC

Incremental Critical Cells (iCC) computes persistent-homology information for
dense Euclidean point clouds. The executable provides three concrete pipelines:

- `ph`: direct implicit Vietoris–Rips persistent homology.
- `piecewise`: piecewise persistent homology (PwPH) with pseudo-vertices (PVs).
- `apparent`: apparent-pair-only reduction, without augmenting-path matching.

This repository contains the C++ implementation used by the HPC lab at the
University of Cincinnati.

## Requirements

- A C++23 compiler
- CMake 3.22 or newer
- OpenMP
- CLI11 (CMake fetches CLI11 v2.6.2 when no installed package is found)

CGAL and Eigen are not required.

## Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable is written to `build/iCC`.

To build and run the core regression tests:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Running

Running without arguments opens the interactive runner:

```bash
./build/iCC
```

Direct persistent homology:

```bash
./build/iCC --tool ph --file-name <input.csv> --homology-dim <homologyDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

PwPH with explicit epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --homology-dim <homologyDim> --eps-breaks <eps1> <eps2> ... --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v [--ffi]]
```

PwPH with automatically generated epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --homology-dim <homologyDim> --eps-interval-count <count> [--eps-interval-scale <scale>] --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v [--ffi]]
```

Apparent-pair-only reduction:

```bash
./build/iCC --tool apparent --file-name <input.csv> --homology-dim <homologyDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

`--homology-dim H` (or `-d H`) sets the maximum homology dimension to compute,
following [Ripser's dimension convention](https://github.com/Ripser/ripser#usage).
The accepted range is `0` through `62`. Internally, `maxdim` remains the maximum
interface dimension and is set to `H + 1`, since deaths in homology dimension
`H` require cofacets of dimension `H + 1`.

For example, `--homology-dim 2` computes through H2 using interfaces through
dimension 3. It replaces the previous `--max-dim 3`; the old long option is no
longer accepted. The interactive runner also asks for the homology dimension.
Existing processing and CSV dimension fields retain their interface numbering.

H0 is still handled by the initial minimum spanning forest and its intervals
are not printed. Consequently, `--homology-dim 0` performs the edge-only run
without reporting positive-dimensional PH intervals.

The apparent pipeline prints deterministic `dimension, unmatched weight` rows.
It includes unmatched top-dimensional cofacets and does not print H0 cells.

The numeric tool aliases `1`, `2`, and `3` select `piecewise`, `ph`, and
`apparent`, respectively. Existing PH aliases remain accepted; the apparent
pipeline also accepts `apparent-pairs` and `apparentpairs`.

For every option and alias:

```bash
./build/iCC --help
```

## Epsilon schedules

`--eps-breaks` and `--eps-interval-count` are mutually exclusive. The aliases
`--eps-break-count` and `--eps-window-count` are accepted for the interval
count. `--eps-break-scale` and `--eps-window-scale` are accepted for the
interval scale.

Automatic scheduling sorts the distinct positive pairwise distances and chooses
the requested number of upper bounds from those ranks. The interval scale must
be finite and at least `1.0`; `1.0` gives linear distance ranks, while larger
values place more distinct distances in earlier windows. Selected distances are
advanced by one representable floating-point value because filtration admission
uses a strict upper bound.

`--pv-cap-scale` limits a new PV's diameter to `scale * eps_hi`, where `eps_hi`
is its formation-window upper bound. `--pv-min-separation` requires every pair
of PVs to have minimum cross-distance at least `delta * eps_max`; its default is
`0`, which disables this constraint.

The executable accepts any finite positive PV cap scale for experiments. The
quotient and single-coordinate safety arguments assume the cap does not exceed
the formation scale. Because admission is strict, use a scale below `1.0` to
avoid the equality boundary.

## Output and timing

Direct PH retains the existing per-dimension persistent-pair report. Nonverbose
PwPH writes a CSV-like window header followed by
`dimension, birth weight, death weight` rows. It also reports the top interface's
cofacet and facet list sizes in each window, without lower-interface counts or
PV details. Verbose direct PH includes decoded
birth-facet and death-cofacet labels for every persistence pair.

Verbose PwPH starts each window with its epsilon bounds. Interfaces are reported
starting at dimension 2, with cofacet, facet, and active-facet list sizes before
their PH results (homology dimension 1 and up). Only the top interface includes
birth/death simplex labels and geometric realizations as original vertices;
weights appear in the PH intervals. A surviving interval has no death simplex.
PV counts, diameters, and original representatives follow the PH results in
non-final windows. The final
window omits PV information because no new PVs are formed. The edge/MST
interface is omitted.

Add `--ffi` to `--verbose` to enable canonical-witness
false-facet-identification statistics for the top interface in PwPH windows
containing PVs in top interface dimensions 3 through 7 (`--homology-dim` 2
through 6). Plain verbose mode does not collect
FFI realization maps or run the FFI report. `--ffi` requires `--verbose`.
Each report describes the complex rebuilt at that window's upper scale; it is
not restricted to cofacets whose weights lie strictly inside the recording
interval.

Before the top-dimensional matching, verbose `ph` and `apparent` runs with
homology dimension at least 1 report the total and active facet counts, cofacet
count, and allocated list storage. They also report the active-facet map
allocation and the estimated saving from an indirect 32-bit rank hash at the
same bucket count. This estimate includes hash metadata and overflow slots;
it excludes allocator overhead and
does not predict whole-run peak RSS.

Every command-line run prints `run time = ...` in milliseconds. Input reading,
distance-matrix construction, and automatic epsilon-schedule generation happen
before this timer starts.
