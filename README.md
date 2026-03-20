# Minimum Time Consensus with energy constraint for Multi-Agent Systems

## 📖 Overview
This repository contains a comprehensive MATLAB framework for analyzing and simulating minimum time consensus under energy constraint for multi-agent systems. The codebase provides algorithms for computing optimal consensus times, defining energy attainable sets (via ellipsoidal approximations), and visualizing state trajectories and control inputs across various system dynamics (e.g., double integrators, third-order systems).

This project is structured to support theoretical and algorithmic research in control theory, specifically addressing the challenges of driving a multi-agent system to consensus under stringent energy constraints and optimal time conditions.

## Author(s)
* **Akansha Rautela^, Deepak Patil^, Ameer Mulla^^, Indra Narayan Kar^** - *^Indian Institute of Technology (IIT) Delhi*, ^^*Indian Institute of Technology (IIT) Dharwad*

## 🗂️ Directory Structure

The repository is organized into a modular structure to separate core algorithms, attainable set computations, utility functions, and simulation results.

```text
organized_minimum_energy/
├── root/                   # Main execution scripts and core case evaluations
│   ├── main_energy_consensus.m
│   ├── case0.m - case4.m
├── src/                    # Source code modules
│   ├── attainable_sets/    # Scripts for computing 2D/3D intersection of reachable sets
│   ├── core_algorithms/    # Core simulation loops, consensus logic, and system dynamics
│   ├── tests/              # Unit tests for specific functional cases
│   └── utils/              # Plotting, projections, and visualization helpers
├── results/                # Generated output data
│   ├── data_logs/          # Raw text logs of coefficients, times, and states
│   └── figures/            # MATLAB figures (.fig) and publication-ready vector graphics (.eps)
├── docs/                   # Documentation and generated HTML files
└── backups/                # Auto-saved and backup script files
