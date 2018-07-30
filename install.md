DAFMPB depends on two external libraries: DASHMM and HPX-5. DASHMM leverages 
the global address space of the HPX-5 runtime system to provide a unified
evaluation of the multipole methods on both shared and distributed memory computers. 
This enables the latest version of AFMPB to operate on distributed memory 
computers while at the same time maintaining backward compatibility on shared
memory computers. 

Version 4.1.0 of HPX-5 is provided in 'contrib' directory. DASHMM is 
automatically downloaded by AFMPB when the application is built. 

Users must install HPX-5 on their systems before installing the DAFMPB solver. For users
who use DAFMPB on shared memory computers only, HPX-5 can be built in the following
steps 
```
> cd /path/to/hpx
> ./configure --prefix=/path/to/install
> make
> make install
```
For users who use DAFMPB on distributed memory computers, HPX-5 currently specifies two 
network interfaces to choose from: 
1. the `ISend/IRecv` interface with the MPI transport
2. the _Put-with-completion_ (PWC) interface with the Photon transport. 
HPX-5 can be built with either transport. 

To configure HPX-5 with MPI network, one adds `--enable-mpi` to the configure line. 
The configuration will search for the appropriate way to include and link to MPI 
1. HPX-5 will try and see if `mpi.h` and `libmpi.so` are available with no additional flags. 
2. HPX-5 will test for an `mpi.h` and `-lmpi` in the current `C_INCLUDE_PATH` and `{LD}_LIBRARY_PATH`. 
3. HPX-5 will look for an `ompi pkg-config` package. 

To configure HPX-5 with the Photon network, one adds `--enable-photon` to the 
configure line. HPX-5 does not provide its own distributed job launcher, so it is necesary
to also use either the `--enable-mpi` or `--enable-pmi` option in order to build support
for `mpirun` or `aprun` bootstrapping. Note that if you are building with the Photon network, 
the libraries for the given network interconnect you are targeting need to be present
on the build system. The two supported interconnects are InfiniBand (`libverbs` and `librdmacm`)
and Cray's GEMINI and ARIES via uGNI (`libugni`). On Cray machines you need to include
`PHOTON_CARGS="--enable-ugni"` to the configure line so that Photon builds with uGNI support. 
Finally, the `--enable-hugetlbfs` option causes HPX-5 heap to be mapped with huge pages, 
which is necessary for larger heaps on some Cray Gemini machines. 

Once HPX-5 is installed, the DAFMPB package can be built in the following steps:
```
> git clone git@github.com:zhang416/dafmpb.git
> mkdir dafmpb-build
> cd dafmpb-build
> cmake ../dafmpb 
> make 
```
This put the executable `dafmpb` under `dafmpb/example` directory. 
