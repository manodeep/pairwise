[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/manodeep/pairwise/master/LICENSE)
[![DOI](https://zenodo.org/badge/19184/manodeep/pairwise.svg)](https://zenodo.org/badge/latestdoi/19184/manodeep/pairwise)


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
of 1 million squared-distances, averaged over (at most) 1000 iterations. Each test is
repeated 5 times and the best average time is noted here. 

Function                         |  gcc timings (ms)| icc timings (ms) | clang timings (ms)|
---------------------------------|------------------|------------------|-------------------|
naive                            |   0.56 +-  0.01  |   0.53 +-  0.00  |   0.53 +-  0.01   | 
chunked                          |   0.57 +-  0.00  |   0.57 +-  0.01  | 	 0.55 +-  0.01   |
compiler_vectorized_chunked      |   1.45 +-  0.01  |   0.60 +-  0.01  | 	 1.54 +-  0.00   |
intrinsics_chunked               |   0.50 +-  0.01  |   0.63 +-  0.01  | 	 0.52 +-  0.01   |
intrinsics_chunked_unroll2       |   0.49 +-  0.00  |   0.59 +-  0.01  | 	 0.53 +-  0.00   |
pairwise_ispc                    |   0.62 +-  0.01  |   0.63 +-  0.01  | 	 0.63 +-  0.01   |


Now, the same table just that the distances are square-rooted. Since, square-root is a
very expensive operation, the timings are about a factor of 5 larger on average. 


Function                         |  gcc timings (ms)| icc timings (ms) | clang timings (ms)|
---------------------------------|------------------|------------------|-------------------|
naive                            |   3.19 +-  0.02  |   3.19 +-  0.03  |   3.19 +-  0.00   | 
chunked                          |   3.19 +-  0.02  |   3.19 +-  0.01  | 	 3.19 +-  0.02   |
compiler_vectorized_chunked      |   3.74 +-  0.02  |   3.81 +-  0.01  | 	 4.75 +-  0.00   |
intrinsics_chunked               |   3.19 +-  0.03  |   3.19 +-  0.03  | 	 3.19 +-  0.02   |
intrinsics_chunked_unroll4       |   3.27 +-  0.02  |   3.21 +-  0.03  | 	 3.19 +-  0.02   |
pairwise_ispc                    |   3.19 +-  0.02  |   3.19 +-  0.02  | 	 3.19 +-  0.00   |

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


Function                         |  gcc timings (ms)| icc timings (ms) | clang timings (ms)| 
---------------------------------|------------------|------------------|-------------------|
naive                            |  2.12 +-  0.02   |   2.23 +-  0.03  |   1.97 +-  0.00   |
chunked                          |  1.42 +-  0.01   |   1.50 +-  0.01  | 	 1.44 +-  0.01   |
compiler_vectorized_chunked      |  2.19 +-  0.02   |   2.21 +-  0.01  | 	 1.90 +-  0.01   |
intrinsics_chunked               |  0.60 +-  0.00   |   0.68 +-  0.01  | 	 0.71 +-  0.00   |
intrinsics_chunked_unroll2       |  0.55 +-  0.01   |   0.67 +-  0.02  |   0.63 +-  0.01   |

# Author

Pairwise is written/maintained by Manodeep Sinha. Please contact the [author](mailto:manodeep@gmail.com) in
case of any issues.

# Citing the code

If you use the code, please cite using the Zenodo DOI. The BibTeX entry is:

```
@misc{manodeep_sinha_2015_33657,
  author       = {Manodeep Sinha},
  title        = {pairwise: v1.0},
  month        = nov,
  year         = 2015,
  doi          = {10.5281/zenodo.33657},
  url          = {http://dx.doi.org/10.5281/zenodo.33657}
}
```

# LICENSE

Pairwise is released under the MIT license. Basically, do what you want
with the code including using it in commercial application.

# Project URL
 
* version control (https://bitbucket.org/manodeep/pairwise)
