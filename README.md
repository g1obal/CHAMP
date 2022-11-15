# CHAMP  
[Cornell-Holland Ab-initio Materials Package](https://cyrus.lassp.cornell.edu/champ)

CHAMP is a quantum Monte Carlo suite of programs for electronic structure 
calculations on a variety of systems (atoms, molecules, clusters, solids and 
nanostructures) principally written by Cyrus Umrigar, Claudia Filippi, and 
Julien Toulouse with major contributions by Devrim Guclu and Frank Petruzielo 
and smaller contributions by other postdocs, students, and collaborators.

## Fork Information
Currently, CHAMP is maintained and distributed via a Subversion repository. 
If you wish to use the most recent version of this program, please contact one 
of the principal authors.

This fork is created for 2D systems of quantum dots starting from CHAMP version 
3.08.2 (SVN Revision: 733). It is specifically modified and tested for finite
sized artificial superlattices. 

## Installation
Installation requires *gcc*, *gfortran*, *openmpi* and *make* (or corresponding 
Intel compilers) to be installed in your system.

Download or clone this repository and navigate to *CHAMP* directory.
```
$ git clone https://github.com/g1obal/CHAMP.git
$ cd CHAMP
```

Compile the libraries (see *lib/0Info* for information).
```
$ cd lib
$ make
$ cd ..
```

Navigate to *qmc* directory and compile the program (serial or MPI).
```
$ cd qmc
```

For the serial version: 
```
$ make clean
$ make
```

For the MPI version:
```
$ make clean
$ make mpi
```

Between serial and MPI version compilations, `$ make clean` should be done so 
that MPI and serial objects are not mixed.

If a file is added/deleted or module dependencies are modified, generate a new 
makefile using: `$ make make`

For more information see *documentation* directory.

## Usage
For the serial version:
```
$ champ.exe -m vmc_mov1 < input > output
```
```
$ champ.exe -m dmc_mov1 < input > output
```

For the MPI version:
```
$ mpirun -np 2 champ_mpi.exe -m vmc_mov1_mpi -i input > output
```
```
$ mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi1 -i input > output
$ mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi2 -i input > output
$ mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi3 -i input > output
```

DMC MPI information: <br />
For mpi1 runs, the population control is on each process; the population size is
the number of walkers per process. <br />
mpi2 and mpi3 runs distribute the walkers to all processes, thus the population
size is the total number of walkers on all processes.
