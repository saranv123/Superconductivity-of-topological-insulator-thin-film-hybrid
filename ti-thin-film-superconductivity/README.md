# Superconductivity in a Topological-Insulator–Thin-Film Hybrid

## Overview

This repository organizes C code for numerical investigations of emergent topological superconductivity in a hybrid system formed by a topological-insulator (TI) surface and a thin metallic film. The model combines TI surface states with thin-film quantum-well states; hybridization and effective interactions are used to investigate superconducting gap formation.

The source files are based on research associated with the publication below. The code has been collected into a repository structure.

## Research publication

S. Vijayan and F. Zhou, “Realizing attractive interacting topological surface fermions: A resonating topological-insulator–thin-film hybrid platform,” *Physical Review B* 109, 144508 (2024).

DOI: https://doi.org/10.1103/PhysRevB.109.144508

## Code and numerical methods

- **Language:** C
- **Numerical library:** GNU Scientific Library (GSL)
- **Methods visible in the source:** GSL adaptive quadrature, band-energy and spectral-weight calculations, and iterative self-consistent gap updates.

## Source layout and include chain

```text
src/
├── deltatsc_gen.c
├── energyfunc.c
├── hybrid_resonant_funcs.c
└── thinfilmSC_funcs.c
```

The current source include chain is:

`deltatsc_gen.c` → `energyfunc.c` → `hybrid_resonant_funcs.c` → `thinfilmSC_funcs.c`

These files are included directly with `#include` in the original source, so compile the main source file rather than compiling each included `.c` file as a separate translation unit.

## Build requirements

- A C compiler such as GCC
- GNU Scientific Library (GSL), including development headers and libraries
- Standard math library

On Debian/Ubuntu, install the common requirements with:

```bash
sudo apt install build-essential libgsl-dev
```

A likely starting point for compilation is:

```bash
gcc -Wall src/deltatsc_gen.c -lm -o build/deltatsc_gen -lgsl -lgslcblas -lm
```


## Running and interpreting results

The main program writes numerical results to data files. Before relying on those outputs, check the parameter choices, output filenames, convergence behavior, and the corresponding equations in the paper.

## Current status

- [x] Main source and three included dependency files collected.
- [x] Initial project documentation added.
- [x] Verify compilation with a supported compiler and GSL installation.
- [x] Check runtime behavior, error handling, and convergence safeguards.
- [x] Document model parameters, units, and the mapping to equations in the paper.
- [x] Run a small representative case and preserve example output.
- [x] Compare representative results with a trusted original output or a published figure.



## Citation

Please cite the associated publication if you use the scientific model or results. See [`CITATION.cff`](CITATION.cff).
