GLY example

# run with built-in mesh generation
./afmpb --pqr-file=GLY.pqr 

# run with provided mesh file 
./afmpb --pqr-file=GLY.pqr --mesh-format=1 --mesh-file=GLY.pqr-mesh.dat-d20-r0.5

fas2 example

# run with provided mesh file, result is saved in fas2.out
./afmpb --pqr-file=fas2.pqr --mesh-format=2 --mesh-file=fas2.off 

