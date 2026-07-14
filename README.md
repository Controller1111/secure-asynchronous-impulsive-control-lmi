# MATLAB LMI Codes for Secure Asynchronous Impulsive Control

This repository provides the MATLAB LMI codes used to verify the feasibility conditions reported in the simulation examples of the manuscript entitled:

**Secure Asynchronous Impulsive Control of Nonlinear Multi-agent Systems under Deception Attacks**

## Requirements

- MATLAB
- MATLAB LMI Toolbox

## Files

- `Example1_Case1_c0439_LMI.m`  
  LMI feasibility test for Case 1 of Example 1 with `c = 0.439`.

- `Example1_Case1_c0990_LMI.m`  
  LMI feasibility test for Case 1 of Example 1 with `c = 0.990`.

- `Example1_Case2_c0680_LMI.m`  
  LMI feasibility test for Case 2 of Example 1 with `c = 0.680`.

- `Example2_LMI.m`  
  LMI feasibility test for Example 2 with `c = 0.3`.

## Usage

1. Install MATLAB with the MATLAB LMI Toolbox.
2. Download or clone this repository.
3. Open MATLAB and set the current folder to the directory containing the `.m` files.
4. Run the file corresponding to the desired simulation case.

Each file formulates the corresponding LMI conditions, checks their feasibility, and extracts the decision matrices `P_j` and `Q_j`.

## Reproducibility

The parameter settings are consistent with those reported in the manuscript. In `Example2_LMI.m`, a fixed random seed is used to reproduce the impulse index sequence adopted in Example 2.
