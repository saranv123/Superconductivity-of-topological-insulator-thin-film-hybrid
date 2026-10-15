# Numerical methodology

## Scientific context

The associated publication studies a hybrid formed from a topological-insulator surface and a thin metallic film. The source code evaluates hybrid-band properties and superconducting gap quantities in this setting.

Reference: S. Vijayan and F. Zhou, *Physical Review B* **109**, 144508 (2024), https://doi.org/10.1103/PhysRevB.109.144508.

## Code organization

- `deltatsc_gen.c`: main numerical calculation; defines integrands and calls GSL integration routines.
- `energyfunc.c`: band-energy functions, spectral weights, interaction terms, and integration-limit helpers.
- `hybrid_resonant_funcs.c`: hybrid-system quantities including detuning, effective spin-orbit coupling, spectral weights, and density of states.
- `thinfilmSC_funcs.c`: thin-film interaction and gap-related functions, including iterative updates of gap values.

The source files use nested `#include` directives to include one another. The main source should therefore be compiled as the translation unit; do not separately compile all four `.c` files and link them together without first refactoring the include structure.

## Numerical workflow

1. Calculate hybrid-band energies and spectral weights from model parameters.
2. Determine momentum-space integration limits using helper routines.
3. Start with a trial value for superconducting gap.
4. Evaluate integrals with GSL adaptive quadrature.
5. Update superconducting gap values iteratively.
6. Use relative changes in gap values as convergence criteria.
7. Stop the iteration when successive gap values fall within the convergence criteria. 
8. Write calculated quantities to data file.


This is a high-level description of the implementation visible in the source. It is not a substitute for mapping each routine to the exact equations and conventions in the paper.

## Verification checklist

- [ ] Build with a documented compiler and GSL version.
- [ ] Resolve compiler warnings and confirm standard-compliant program entry point.
- [ ] Check all allocations, file operations, GSL status codes, and numerical edge cases.
- [ ] Add iteration limits and clear failure reporting where needed.
- [ ] Explain all model parameters and units.
- [ ] Map routines to equations in the paper.
- [ ] Run a small test case and record parameters and output.
- [ ] Compare representative results with a trusted original output or published figure.

