# TUPS
This Code realizes an implementation of the Truncated Unity parquet equations and is thus called  **T**runcated **U**nity **P**arquet **S**olver (TUPS).
As a reference, also checkout the corresponding  [paper](https://arxiv.org/abs/1912.07469).
For more information on the Truncated Unity parquet equations in general also see  [here](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.075143).
The code can be used under the GNU General Public License v3.0.

## Installation / Compilation
This package includes a makefile.
Compilation with mpif90 has been tested.
Further required kernels are BLAS and FFT.
We include a FFT implementation in ./lib which one can simply link to, alternatively one could also link to another implementation (eg. MKL).
For BLAS we link to intel's MKL library - linking to a different implementation should also be straight forward through altering the makefile.


## Execution
The code can be executed via mpiexec.
As the first commandline argument the parameter file must be given - also see [Input](#Input).
The call would thus for example be
```shell
$mpiexec -n 5 ./TUPS prms
```
if `prms` is the parameter file and one spawns 5 mpi tasks.
The repository contains an example parameter file called `prms`.

## Input
The input to the program is defined in an extra file.
The different options and parameters are explained in the sample input file `prms`.
Some of the options will be further discussed here.
In general the readin is perform in `parquet_ini.f90` which can also be used as a reference.

### picklist
The variable 'picklist' lets one choose the used form factors by hand through a file containing booleans.
This file must have as many lines as there are k-points so for example 64 for a 8x8 grid.
The form factors can be enumerted by their distance to the origin.
If the distance of two form factors to the origin equals their order is defined through the functions in `parquet_formfactors.f90`.
Each line in the picklist-file should have an entry 'T' if one wants to include the form factor or 'F' if not.
Note that this option is not compatible with the option 'old = .false.' which realizes a more optimized implementation of the parquet equations.
It is, however, compatible with the option 'SDE_old = .false.' since the more optimized implementation of the Schwinger Dyson Equation (SDE) is compatible with taking arbitrary form factors.

### DGA
Takes a boolean/logical value and decides whether one performs calculations in parquet approximation (PA) or Dynamical Vertex Approximation (DGA).
Thus one either uses *U* as fully irreducible vertex in which case no further input file is required or one uses the fully irreducible vertex of an auxillary impurity problem in which case one needs to specify the corresponding file - which is done in the next line of the input file. 




 

