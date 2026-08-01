# Monte Carlo Simulations for Coarse-Grained DNA Braiding

This repository contains MATLAB code for performing Monte Carlo (MC) simulations of coarse-grained DNA braiding systems.

## Requirements

The code was developed and tested with:

- MATLAB R2024b
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox

## Repository Structure

### `running_general_braiding_simulation.m`

Main driver script for setting up and running a DNA braiding simulation.

This script specifies:

- Worm-like chain (WLC) parameters
- Geometrical constraints
- Mechanical parameters (e.g., force, turns)
- Monte Carlo simulation parameters

---

### `braiding_init_configs.m`

Generates the initial DNA braid configuration while ensuring the correct geometry and topology of the system.

---

### `braiding_simu_general.m`

Core Monte Carlo simulation routine.

This function generates equilibrium DNA braid conformations using the Monte Carlo algorithm and evaluates the corresponding system energies.

---

### `interpolant_acos_separation_length0.64312nm_veff_5.0803epernm.mat`

Lookup table containing interpolated Debye–Hückel electrostatic interaction values used during energy calculations.

The interpolant computes the electrostatic repulsion between two rod segments based on their relative positions and orientations.

For details of the underlying method, see:

> Klenin, K. et al., *Biophysical Journal* **74**, 780–788 (1998).

## Notes

The simulations employ a coarse-grained worm-like chain (WLC) representation of DNA and use Monte Carlo sampling to generate equilibrium braid conformations under prescribed mechanical constraints.
