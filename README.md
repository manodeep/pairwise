# Description

This repo contains a set of codes to benchmark various methods that
generate all pairwise separations. The work was inspired
(read, challenge accepted) by the blog post by Jake VanderPlas
[here](https://jakevdp.github.io/blog/2013/06/15/numba-vs-cython-take-2/).



## AllPairwise

The main code is ``pairwise_3d.c`` -- this contains a set of functions
that compute all pair-wise separations for two sets of ``N`` random
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

### Performance Table

Performing 1000 iterations over 2 sets with 1000 elements each. So, a total computation
of 1 million squared-distances, averaged over 1000 iterations. 

Function                         |  gcc timings (ms)| icc timings (ms) | clang timings (ms)
---------------------------------|------------------|------------------|--------------------
naive                            |   0.57 +- 0.00   |   0.55 +- 0.01   |   0.57 +- 0.02 
chunked                          |   0.59 +- 0.00   |   0.58 +- 0.01   | 	 0.56 +- 0.02
compiler_vectorized_chunked      |   1.56 +- 0.01   |   0.61 +- 0.01   | 	 1.56 +- 0.02
intrinsics_chunked               |   1.19 +- 0.02   |   1.18 +- 0.02   | 	 1.30 +- 0.08
intrinsics_chunked_unroll4       |   1.18 +- 0.02   |   1.18 +- 0.02   | 	 1.18 +- 0.02
pairwise_ispc                    |   0.63 +- 0.00   |   0.63 +- 0.01   | 	 0.63 +- 0.01


Now, the same table just that the distances are square-rooted. Since, square-root is a
very expensive operation, the timings are about a factor of 5 larger on average. 


Function                         |  gcc timings (ms)| icc timings (ms) | clang timings (ms)
---------------------------------|------------------|------------------|--------------------
naive                            |   3.21 +- 0.05   |   3.21 +- 0.02   |   3.19 +- 0.02 
chunked                          |   3.24 +- 0.05   |   3.22 +- 0.07   | 	 3.19 +- 0.01
compiler_vectorized_chunked      |   4.78 +- 0.06   |   3.81 +- 0.03   | 	 4.75 +- 0.01
intrinsics_chunked               |   3.20 +- 0.02   |   3.20 +- 0.02   | 	 3.20 +- 0.02
intrinsics_chunked_unroll4       |   3.29 +- 0.02   |   3.23 +- 0.01   | 	 3.20 +- 0.02
pairwise_ispc                    |   3.20 +- 0.02   |   3.20 +- 0.02   | 	 3.19 +- 0.01

Expectedly, C is much faster than any numba/cython/f2py combination can
get you. For comparison, the fastest timings with numba/cython/f2py was around 9 ms. 
See [here](https://jakevdp.github.io/blog/2013/06/15/numba-vs-cython-take-2/) for
a detailed blog post. 

Note, that you have CPU frequency scaling on *and* do not run the tests for long enough,
then you might be fooled into thinking that the tests take longer to run (simply because
the cpu is running at a lower frequency). I found out by ``diff``ing between 
linux ``perf`` output for slow and fast runs.  


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