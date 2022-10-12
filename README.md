# fomc

Classical molecular force fields for atomistic modeling via Molecular Mechanics/Dynamics
or Monte Carlo approaches come with parameters that can be obtained in a number of ways.
While parameters for intramolecular force field terms (bonds, angles, dihedrals) as well
as partial charges can be obtained from first principles, VdW interactions are
often parameterized by empirical fitting to reproduce experimental data.
fomc is a program that performs this type of optimization, but, different from
most common approaches, the experimental quantities it uses for optimization
are not interaction energies, sublimation energies, or the like, but crystal
structures. fomc starts with a given set of Lennard Jones (LJ) parameters that can
be taken from a general purpose force field like GAFF. It then performs a brief
MD simulation of a sample of the crystal structure(s) (the empirical input data) using
Gromacs as external process. After comparing the input structure (the ground truth) to
the average structure from the MD simulation, one of the LJ parameters is modified
randomly within certain boundaries, and the change is accepted (or not) based on a Monte
Carlo criterion. The process is repeated for a given number of cycles or until convergence
is reached.

## Features

- Takes the crystal structure of a small organic molecule, and an initial force field parameter set
(provided in Gromacs topology format) as input, and adjusts the parameters so that the
reproduction of the crystal at a given temperature is optimized.
- The optimization algorithm is a straight forward Monte Carlo algorithm
- Molecular Dynamics simulations in the NPT ensemble allow for an optimization of parameters to
reproduce structural properties at a given temperature.
- A minimum of one crystal structure is required as input, but more than one can also be
provided (e.g. different polymorphs), which should provide parameters with improved
transferability.
- The total number of parameters is limited using Lorentz Berthelot rules.

## Prerequisites

fomc is a plain C code that requires only standard libraries, and should compile
with gcc on most Linux distributions. For execution a working installation of Gromacs
(any version newer than 2016 should work) is required. Other than that, for the generation
and manipulation of topologies and input structures, a number of external tools can be
helpful, and are listed below.

- [Gromacs](http://www.gromacs.org/)
- [openbabel](http://openbabel.org/)
- [acpype](https://github.com/alanwilter/acpype)
- [Ambertools](http://ambermd.org/AmberTools.php)
- [gdis](https://github.com/arohl/gdis)

For more details see the documented [workflow](examples/README.md) in the examples folder.

## Installation

```
git clone https://github.com/mbatgh/fomc.git
cd fomc/src
make 
make install
```
should generate a binary called fomc in the bin folder in your home directory.

## Usage

fomc is executed on the command line. Various settings can be controlled in a parameter input file,
or with command line parameters. For a concise overview consider the output of the command fomc -h,
as shown below.


```
fomc -h

This is fomc (force field parameter optimization via MC/MD
Version 1.0

usage: fomc [command-line parameters]

  -h          show this msg
  -i  string  id, basename of input files: parameter inp, xtal struc pdb
              as in: id.inp, id-01.pdb to id-0x.pdb
  -s  int     seed for random generator
  -n  int     number of xtal polymorphs
  -r  string  run id, used in combination with random seed for unique output filenames
  -t  float   MC temperature
  -d  float   RMSD threshold [nm]
  -l  float   lattice parameter deviation threshold [%]
  -z  int     timeout for a single MD run in seconds
  -y  int     total number of MC cycles
  -x  float   max dtot at which to stop (usually 2)
  -a  string  input atp filename
```

For more details on the generation of input files, their syntax, and the overall workflow, see the 
example provided [here](examples/README.md)

## Contact

Please send any questions/bug-reports/suggestions to me. My full name is rare enough, so that,
together with the place where I live (Graz) google will promptly provide you with my email address.

## Reference

If you use this code in published work, please cite the reference given below.

Michael Brunsteiner, Sten Nilsoon-Lill, Lucy M. Morgan, Adrian Davis, and Amrit Paudel
Finite temperature mechanical properties of organic molecular solids from classical molecular simulation.
To be submitted (2022).

## Acknowledgements

The work involved in the design of this tool was done at the [Research Center Pharmaceutical Engineering](http://www.rcpe.at),
in a collaborative K project, towards understanding chemical stability of small molecule drugs in the solid state, with funding
from the [Austrian COMET program](https://www.ffg.at/en/comet/programme), and a number of companies (Astra-Zeneca, Janssen,
Pfizer, and UCB).

