# Description

This repo contains a set of codes to benchmark various methods that
generate all pairwise separations. The work in here was inspired
(read, challenge accepted) by the blog post by Jake VanderPlas
[here](https://jakevdp.github.io/blog/2013/06/15/numba-vs-cython-take-2/).

Expectedly, C is much faster than any numba/cython/f2py combination can
get you.


## AllPairwise

The main code is ``pairwise_3d.c`` -- this contains a set of functions
that compute all pair-wise separations for a set of ``N`` random
floating point numbers. 

``pairwise_3d.ispc`` contains the [ispc](http://ispc.github.io/) code to compute
all pairwise separations. Spoiler alert: ispc seems to produce the fastest
function for a variety of N. 

``pairwise_3d.f90`` contains the fortran routine to compute the same
pairwise separations in fortran. Taken straight from the comments in the
blog link above and edited to show reasonable run times. 


``pairwise.c`` contains similar functions as ``pairwise_3d.c`` except that the
array storage is now ``xyzxyz...`` rather than ``xxx...yyy...zzz..`` . This
``pairwise`` also contains the oft-touted BLAS way of computing pairwise
separations. Notably, ``pairwise`` is about 2-3 times slower than ``pairwise_3d``,
agreeing with what we expect from array of structure vs structure of arrays 
arguments (AOS vs SOA). 

## Pairwise Histogram

In progress. Ultimately, the repo will contain a similar idea, just for the
histograms of all pairwise distances. The fastest codes will probably get
migrated into my correlation function routines [here](https://bitbucket.org/manodeep/corrfunc/).

# Author

Pairwise is written/maintained by Manodeep Sinha. Please contact the [author](mailto:manodeep@gmail.com) in
case of any issues.

# LICENSE

Pairwise is released under the MIT license. Basically, do what you want
with the code including using it in commercial application.

# Project URL
 
* version control (https://bitbucket.org/manodeep/pairwise)