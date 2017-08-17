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
