# 2D-Incompressible-Navier-Stokes-Solver
here I have added the source code which i have coded myself for this project


# 2D Incompressible Navier-Stokes Solver

A high-performance 2D computational fluid dynamics (CFD) solver built from scratch in **Fortran**. This solver is designed to accurately simulate complex fluid dynamics, including the classic lid-driven cavity and the transient vortex shedding (Von Kármán street) of a square cylinder. 

Heavy emphasis is placed on hardware-accelerated computing, utilizing **OpenMP** for multi-core CPU parallelization and **OpenACC** for GPU offloading.

## Features

* **Numerical Method:** Fractional-step (projection) method for pressure-velocity coupling.
* **Grid Topology:** Staggered Cartesian grid architecture to prevent pressure checkerboarding.
* **Hardware Acceleration:** 
  * OpenMP implementation for shared-memory CPU architectures (achieving up to **3.85x speedup**).
  * OpenACC directives for highly parallelizable kernel execution on GPUs (achieving up to **3.26x speedup**).
* **Convergence:** Strict L-infinity norm criteria (> 1e-6) for momentum and Pressure Poisson solvers.
* **Validation:** Drag coefficients and flow topology rigorously validated against ANSYS Fluent.

## Governing Equations

The solver implements the incompressible Navier-Stokes equations using a fractional-step projection method. An intermediate velocity field is computed first, followed by solving the Pressure Poisson Equation (PPE) to enforce mass conservation:

∇²p = (ρ/Δt) ∇·u*

The velocities are then corrected using the updated pressure gradient.

## Repository Structure

```text
├── src/                  # Fortran source code (.f90 files)
└── README.md             # Project documentation
