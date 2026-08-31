# FVM-RG-Solver (Methane Shockwave Branch)

This branch contains a numerical solver and the associated models for
the paper:
**"Simulation of Shock Waves in Methane: A Self-Consistent Continuum
Approach Enhanced Using Machine Learning"**
*Mathematics, 2024, 12(18), 2924.*

**DOI:** [10.3390/math12182924](https://doi.org/10.3390/math12182924)

## Overview

This work numerically models the internal structure of a one-dimensional
shock wave propagating through single-component methane, a polyatomic gas
whose internal (vibrational) energy relaxation is significantly more
complex than that of diatomic gases. The flow is described using
conservation equations for mass, momentum, and energy, closed with
transport coefficients (viscosity, thermal conductivity, bulk viscosity)
derived from kinetic theory rather than empirical correlations, and solved
using the finite-volume method with the Godunov scheme and HLLE
approximate Riemann solver.

Understanding methane shock structure is relevant to entry into planetary
atmospheres containing methane — particularly Titan. This single-species study is
presented as a necessary foundation before extending the model to
multi-species mixtures (e.g., methane–nitrogen) for future planetary
mission-relevant work.

The paper was presented at open-science event YSM 2024 (All-Russian Conference of Young Scientists).

## Scope of This Branch-Repository

To be precise about what is and is not modeled here:

- **What this solves:** 1D shock wave structure (pre-shock to
  post-shock relaxation zone) for single-component methane, including
  vibrational energy relaxation, using a continuum (Navier–Stokes-level)
  description with kinetic-theory transport coefficients.
- **What this does not solve:** this is not a combustion or detonation
  chemistry solver — no chemical reaction network, ignition, or
  multi-species reacting mixture is modeled in this branch. There is no
  geometry, vehicle, or engine model of any kind; the domain is a 1D flow
  behind a normal shock.
- **Machine-learning component:** part of the solver's cost lies in
  evaluating methane's vibrational energy and specific heat, each of which
  requires summing over ~2000 vibrational energy levels per grid cell,
  per timestep. This repository includes regression models (classical
  and a feedforward neural network) trained to approximate these two
  quantities, integrated into the solver as a drop-in replacement for the
  direct summation, to reduce computational cost without changing the
  underlying physical model.

## Branch-Repository Structure
- `src`: Core finite-volume solver (C++), implementing the governing
  equations, transport coefficients, and Riemann solver.
- `hdr`: Header files, including the neural-network-based approximations
  for vibrational energy and specific heat.
- `example-shockwave`: Setup and configuration for the 1D methane shock
  wave case described in the paper.
- `plot/CreatePlots`: Python scripts for post-processing and generating
  result plots.
- `FVM-RG-Solver-struct.pdf`: Diagram of the solver's code structure.

## Technologies Used
- C++ (core solver)
- CMake (build system)
- Python, PyTorch, scikit-learn (model training for vibrational
  energy/specific heat approximation)

## Notes on Methodology

Direct evaluation of methane's vibrational energy and specific heat via
summation over its vibrational spectrum is the dominant computational cost
in the solver. Several nonlinear regression approaches (k-NN, Decision
Tree, Random Forest, Gradient Boosting, and a single-hidden-layer
feedforward neural network in PyTorch) were compared for approximating
these quantities; the FNN was selected and integrated into the solver,
giving a substantial reduction in computation time while keeping relative
error low.

## How to Cite This Work

**Journal Article:**
```bibtex
@article{maksudova2024methane,
  title={Simulation of Shock Waves in Methane: A Self-Consistent Continuum Approach Enhanced Using Machine Learning},
  author={Maksudova, Z. and Shakurova, L. and Kustova, E.},
  journal={Mathematics},
  volume={12},
  number={18},
  pages={2924},
  year={2024},
  publisher={MDPI},
  doi={10.3390/math12182924}
}
```

## License
This project is licensed under the MIT License — see the LICENSE file for
details.
