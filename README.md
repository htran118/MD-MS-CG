# Molecular dynamics with Multiscale Coarse-Graining

MD-MS-CG is a object-oriented simulation package for Molecular Dynamics (MD) with support for multiscale coarse-graining (CG).

## Introduction

In MD, many accurate atomistic models have been developed. However, because of the sheer number of atoms, it is inefficient to simulate many important processes with a full atomistic system. In these situations, CG simulations, which consider the overall motions of clusters of atoms and ignore less relevant details, are more efficient and useful. Mathematically, it means one needs to average the total external force acting on the CG group.

However, it is well-known that this simple method results in incorrect dynamics, due to the lower number of degrees of freedom. Intuitively, as there are less objects in the system, they will be able to move around more easily. To compensate for the inaccuracy, one can introduce a new, stochastic force. It has been shown by [Noid et al.](https://aip.scitation.org/doi/10.1063/1.2938860) that this force is quite simple and can be easily computed from the pre-CG simulation. They have also showed their mathematical model is correct in some very elementary cases.

There are many popular packages for MD simulations, including LAMMPS, a open-source, object-oriented program in C++. Nevertheless, owing to the need for performance mentioned above, the inner working of these packages is complicated and hard to be modified. On the other hand, the implementation of multiscale CG requires the complete overhaul of some crucial machineries, where most of the users don't ever change. Thus, with more transparent design, support for basic MD techniques and barebone optimization, MD-MS-CG was developed as the solution to these diffculties, 

## Features

- Support for multiscale CG dynamics.
- Support for most of standard MD techniques, including periodic boundary conditions, Verlet's list, Andersen-Lowe thermostat.

## Prerequisites

- Unix-like OS
- G++11
- Make

## Installation

At the main folder, simply run `make all`.

## Instruction

Run the main simulation with simula.exe. For further exaplanations on the role of each file and the structure of the project, refer to README.txt in the main folder.
