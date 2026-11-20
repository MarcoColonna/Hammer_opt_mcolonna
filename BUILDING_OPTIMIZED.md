# Building and validating the optimized fork

The main `Hammer::Hammer` C++ interface retains the upstream `std::vector`
arguments and return types. Conversion to the optimized containers happens
inside the implementation. The Cython wrapper retains the upstream histogram
interface, including explicit bin edges and its original overflow default.
Internal tensor types still use Boost small vectors; code that directly uses
those types may require additional compatibility work.

## Build on the analysis machine

Load the compiler, CMake, Boost and yaml-cpp environment used for the official
HAMMER build. Use a separate build directory and installation prefix to keep
the official installation available for comparison.

```sh
cmake -S . -B build/optimized \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_INSTALL_PREFIX="$HOME/opt/hammer-optimized" \
  -DWITH_PYTHON=OFF \
  -DWITH_ROOT=OFF \
  -DENABLE_TESTS=ON \
  -DINSTALL_HEADERS=ON
cmake --build build/optimized --parallel 4
ctest --test-dir build/optimized --output-on-failure
```

Enable `WITH_ROOT` if your analysis needs HAMMER's ROOT integration. Configure
`CMAKE_PREFIX_PATH` as needed for dependencies outside standard search paths.
The dependency libraries must match the compiler and target architecture.

After the build and tests pass:

```sh
cmake --install build/optimized
```

Rebuild the analysis against this installation, using the same method used to
select the official HAMMER installation. For a CMake consumer, point
`Hammer_DIR` at the installed directory containing `HammerConfig.cmake`.
Existing compiled analysis binaries must be rebuilt because internal container
types changed. Make sure the runtime also loads the optimized libraries.

## Validate before production use

Run the official and optimized versions on identical input with the same
settings. Compare event weights, histogram contents and uncertainties with
appropriate numerical tolerances, then compare elapsed times using equivalent
Release builds. Passing the indexing tests alone does not establish numerical
equivalence or a speedup for an analysis.

Local checks so far: the bundled Google Test compiles; 10 iterator and
sequential-indexing tests pass; compiler syntax checks pass for `Hammer.cc`,
the C++ Hammer tests, outer-dot tests, and the changed iterator, indexing,
outer-container and dot-product sources. These checks used a temporary config
header with ROOT disabled. A full CMake build, installation test, Python build,
physics validation and performance comparison remain to be run on the analysis
machine.
