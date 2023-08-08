# Alpha-transitions

This code base is accompanying our paper with title "Escape by jumps and diffusion by α-stable noise across the barrier
in a double well potential" by Ignacio del Amo and Peter Ditlevsen.
This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Alpha-transitions

It is authored by Ignacio del Amo.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.


The src folder contains the source code, which can be used to generate transition times for a stochastic process
described by a Langevin equation driven by alpha stable noise. Two parametric families of potentials are 
investigated, a quartic polinomial and a globally Lipschitz approximation of it, and a special case with an 
asymmetric potential. 

The scripts folder contains 3 files that allow to generate some of the figures of the paper. They are also ment to 
be examples on how the functions contained in the src folder are used.

For more information consult the paper or drop me a message.
