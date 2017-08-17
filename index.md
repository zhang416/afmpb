The Adaptive Fast Multipole Poisson-Boltzmann (AFMPB) package computes the 
numerical solution of the linearized Poisson-Boltzmann equation that 
describes electrostatic interactions of molecular systems in ionic solutions. 

The linearized Poisson-Boltzmann equation is reformulated as a boundary 
integral equation and is subsequently discretized using the node-patch scheme. 
The resulting linear system is solved using GMRES. Within each iteration, the 
matrix-vector multiplication is accelerated using the DASHMM library. 



