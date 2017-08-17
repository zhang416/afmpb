The Adaptive Fast Multipole Poisson-Boltzmann (AFMPB) package computes the 
numerical solution of the linearized Poisson-Boltzmann equation that 
describes electrostatic interactions of molecular systems in ionic solutions. 

The linearized Poisson-Boltzmann equation is reformulated as a boundary 
integral equation and is subsequently discretized using the node-patch scheme. 
The resulting linear system is solved using GMRES. Within each iteration, the 
matrix-vector multiplication is accelerated using the DASHMM library. 

### Installation
AFMPB depends on two external libraries: DASHMM and HPX-5. Version 4.1.0 of 
HPX-5 can be downloaded [here](https://hpx.crest.iu.edu/download). DASHMM is 
automatically downloaded by AFMPB when the application is built. 

Users must install HPX-5 on their systems before installing the AFMPB solver. 
Instructions on building HPX-5 are available [here](http://hpx.crest.iu.edu/users_guide#overview). 
Once HPX-5 is installed, the AFMPB package can be built in the following 
steps: 

```
> git clone git@github.com:zhang416/afmpb.git
> mkdir afmpb-build
> cd afmpb-build
> cmake ../afmpb 
> make 
```
This put the executable `afmpb` under `afmpb/example` directory. 

### Job Examples 
The minimum input to `afmpb` is the PQR file. This file is first processed by 
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

### References
* B. Lu, X. Cheng, J. Huang, J. A. McCammon. **AFMPB: An Adaptive Fast Multipole
Poisson-Boltzmann Solver for Calculating Electrostatics in Biomolecular Systems.** _Comput. Phys. Commun._ 
181 (2010) 1150

* B. Zhang, B. Peng, J. Huang, N. P. Pitsianis, X. Sun, B. Lu. **Parallel AFMPB Solver 
with Automatic Surface Meshing for Calculation of Molecular Solvation Free Energy.** 
_Comput. Phys. Commun._ 190(2015) 173 


Benzhuo Lu, Xiaolin Cheng, Jingfang Huang, J. Andrew McCammon, Comput. Phys. Commun. 181(2010)1150

Parallel AFMPB solver with automatic surface meshing for calculation of molecular solvation free energy ☆

Bo Zhanga, Bo Pengb, Jingfang Huangc, Nikos P. Pitsianisd, e, Xiaobai Sune, Benzhuo Lub,
