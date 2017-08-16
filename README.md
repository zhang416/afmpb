## Adaptive Fast Multipole Poisson-Boltzmann Solver
The Adaptive Fast Multipole Poisson-Boltzmann (AFMPB) package computes the 
numerical solution of the linearized Poisson-Boltzmann equation that 
describes electrostatic interactions of molecular systems in ionic solutions. 

The linearized Poisson-Boltzmann equation is reformulated as a boundary 
integral equation and is subsequently discretized using the node-patch scheme. 
The resulting linear system is solved using GMRES. Within each iteration, the 
matrix-vector multiplication is accelerated using the DASHMM library. 

### Installation
DASHMM is built on top of the Asynchronous Multi-Tasking HPX-5 runtime system. 
Users must install HPX-5 on their systems before installing AFMPB. HPX-5 can 
be downloaded from [here](https://hpx.crest.iu.edu/download) and version 4.1.0
is used. Once HPX-5 is installed, the AFMPB package can be built in the 
following steps:

```
> git clone git@github.com:zhang416/afmpb.git
> mkdir afmpb-build
> cd afmpb-build
> cmake ../afmpb 
> make 
```
This put the executable `afmpb` under `afmpb/example` directory. 

### Job Examples 
The minimum input to 'afmpb' is the PQR file. This file is first processed by 
`parsePQR.sh` which extracts the atoms information. 

AFMPB can read meshes generated from [MSMS](https://www.ncbi.nlm.nih.gov/pubmed/8906967)
or [TMSMesh](http://lsec.cc.ac.cn/~lubz/Meshing.html). If no input mesh is 
provided, AFMPB will invoke the built-in surface meshing routine. Some 
examples are 

```
// use built-in meshing routine
> ./afmpb --pqr-file=GLY.pqr.ext

// use MSMS mesh
> ./afmpb --pqr-file=GLY.pqr.ext --mesh-format=1 --mesh-file=GLY.pqr-mesh.data-d20-r0.5 

// use TMSMesh 
> ./afmpb --pqr-file=fas2.pqr.ext --mesh-format=2 --mesh-file=fas2.off 
```

On a cluster with Slurm workload manager, a job script looks like this 
```
#! /bin/bash -l
#SBATCH -p queue
#SBATCH -N 2
#SBATCH -t 00:10:00

srun -n 2 -c 48 ./afmpb --pqr-file=fas2.pqr.ext --mesh-format=2 --mesh-file=fas2.off --hpx-threads=24
```
The `-c` option equals to the number of cores Slurm sees on each compute node and the `--hpx-threads`
option equals the number of physical cores available on the compute node. 
