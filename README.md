# SpinWave Library (Refactored)

A object-oriented library for calculating spin-wave excitations, self-energy corrections, and dynamic structure factors in quantum magnetic systems. This is a complete refactoring of the original codebase with improved safety, modularity, and maintainability.

## Key Features

*   **Modern C++ Design:** Uses C++17 features, smart pointers, RAII, and a clean class hierarchy to manage memory safely and reduce leaks.
*   **Modular Architecture:** Functionality is separated into distinct, focused classes:
    *   **`SpinWave`:** Core spin-wave calculations including Holstein-Primakoff transformation and Bogoliubov diagonalization.
    *   **`MeanField`:** Hartree-Fock decoupling and self-consistent field (SCF) calculations.
    *   **`SelfEnergy`:** One-loop and two-loop self-energy corrections, spectral functions, and dynamic structure factors.
    *   **`SimulatedAnnealing`:** Classical ground state search using simulated annealing with heat-bath sampling.
*   **Operator Algebra:** Classes (`boson`, `boson_`, `Bosons`, `Bosons_`) for symbolic manipulation of bosonic creation/annihilation operators in both real and momentum space.
*   **Linear Spin-Wave Theory:**
    *   Automatic generation of bosonic Hamiltonians from spin interactions.
    *   Bogoliubov transformation with automatic handling of non-positive definite matrices (adds a small gap).
*   **Self-Energy Corrections:**
    *   One-loop (`oneLoop`) and two-loop (`twoLoop`) self-energy calculations.
    *   Optimized vertex calculation with reduced redundant multiplications.
    *   OpenMP parallelization for performance.
*   **Hartree-Fock Self-Consistent Field (SCF):**
    *   Self-consistent calculation of mean-field parameters using the `MeanField` class.
    *   Thermal effects via Bose-Einstein distribution.
*   **Classical Ground State Search:**
    *   Simulated annealing with heat-bath sampling algorithm.
    *   Configurable lattice size, temperature schedule, and Monte Carlo steps.
*   **Structure Factor Calculation:**
    *   Dynamic structure factor `S(q, ω)` including nonlinear corrections.
*   **Path Interpolation:**
    *   Bilinear interpolation of data on the Brillouin zone grid.

## Dependencies

*   **C++ Compiler:** C++20 compatible compiler (GCC, Clang, MSVC).
*   **LAPACK & BLAS:** Required for linear algebra operations.
*   **OpenMP:** Used for parallelization.

## Building the Project

A `CMakeLists.txt` file is provided. To build:

1.  Create a build directory:
    mkdir build && cd build

2.  Configure and build:
    cmake .. && make


## Code Structure

*   `geometry.h`: Defines `Vec2`, `Vec3`, and `Term` template classes for geometric and algebraic structures.
*   `spin.h` / `spin.cpp`: Defines `spin` and `Spins` classes for spin operators and interactions.
*   `Boson.h` / `Boson.cpp`: Core boson operator classes:
    *   `boson_`: Real-space boson operators.
    *   `boson`: Momentum-space boson operators with phase factors.
    *   `Bosons_`: Polynomial of real-space boson operators.
    *   `Bosons`: Polynomial of momentum-space boson operators.
*   `lapack.h` / `lapack.cpp`: LAPACK wrapper with RAII-style memory management.
*   `SpinWave.h` / `SpinWave.cpp`: Main spin-wave class handling:
    *   Holstein-Primakoff transformation.
    *   Bogoliubov transformation.
    *   HMF (mean-field Hamiltonian) management.
*   `MeanField.h` / `MeanField.cpp`: Hartree-Fock decoupling and self-consistent calculations.
*   `SelfEnergy.h` / `SelfEnergy.cpp`: Self-energy calculations, Green's functions, and structure factors with optimized vertex reduction.
*   `SimulatedAnnealing.h` / `SimulatedAnnealing.cpp`: Classical ground state search via simulated annealing.
*   `FMKH.cpp`: Example main program for a ferromagnetic Kitaev-Heisenberg model.

## Performance Notes

This refactored version prioritizes code safety and maintainability over raw performance. Recent optimizations to the vertex calculation logic have significantly improved efficiency:
- Execution time for the Kitaev-Heisenberg test case is now approximately 2.1s.
- Memory management is safer with smart pointers, reducing leak risks.
- The code is more modular and easier to extend or debug.

## Limitations

- Momentum-space self-consistent Hartree-Fock has been removed in this refactoring.
- Some features from the original code (e.g., real-space SCF) are preserved but may be slower.

## New Features

- **Simulated Annealing:** Added classical ground state search using heat-bath Monte Carlo sampling with configurable temperature schedule.
- **Optimized Vertex Calculation:** The `Bosons::reduce()` method and `Bosons::swap()` operations reduce redundant multiplications in self-energy calculations, improving performance.

## Contributing

This is a refactored version of a legacy codebase. Contributions for performance optimizations, bug fixes, and feature additions are welcome.
