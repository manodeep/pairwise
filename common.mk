## Options for all-pairwise separations 
#OPT += -DSQRT_DIST
#OPT += -DUSE_MKL

## Print the output with colors.
## Disable if you are outputting to text file or if
## the terminal does not support colors
OPT += -DCOLORED_OUTPUT

### Options that might be enabled in the future
#OPT += -DUSE_OMP

### Set the compiler -- options are icc/gcc/clang. 
CC=gcc

#### Add any compiler specific flags you want
CFLAGS=

#### Add any compiler specific link flags you want
CLINK=


DO_CHECKS := 1
ifeq (clean,$(findstring clean,$(MAKECMDGOALS)))
  DO_CHECKS := 0
endif

ifeq (distclean,$(findstring distclean,$(MAKECMDGOALS)))
  DO_CHECKS := 0
endif

## Only set everything if the command is not "make clean"
ifeq ($(DO_CHECKS), 1)
  ## Make clang the default compiler on Mac
  ## But first check for clang-omp, use that if available
  UNAME := $(shell uname)

  ## Colored text output
  ## Taken from: http://stackoverflow.com/questions/24144440/color-highlighting-of-makefile-warnings-and-errors
  ## Except, you have to use "echo -e" on linux and "echo" on Mac
  ECHO_COMMAND := echo -e
  ifeq ($(UNAME), Darwin)
    ECHO_COMMAND := echo
  endif
  ccreset :=$(shell $(ECHO_COMMAND) "\033[0;0m")
  ccred:=$(shell $(ECHO_COMMAND) "\033[0;31m")
  ccmagenta:=$(shell $(ECHO_COMMAND) "\033[0;35m")
  ccgreen:=$(shell $(ECHO_COMMAND) "\033[0;32m")
  ccblue:=$(shell $(ECHO_COMMAND) "\033[0;34m")
  ## end of colored text output

  ifeq ($(UNAME), Darwin)
    CC := clang
  endif

  ### You should NOT edit below this line
  INCLUDE= -I include -I utils -I.

  ### The POSIX_SOURCE flag is required to get the definition of strtok_r
  CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=c99 -g -m64 -D_POSIX_SOURCE -D_DARWIN_C_SOURCE -O3 -Ofast
  GSL_CFLAGS := $(shell gsl-config --cflags) 
  GSL_LIBDIR := $(shell gsl-config --prefix)/lib
  GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 

  ISPC_AVAIL :=0
  ISPC_VER := $(shell ispc --version 2>/dev/null)
  ifdef ISPC_VER
    ISPC_AVAIL := 1
  endif

  ifeq (icc,$(findstring icc,$(CC)))
    CFLAGS += -xhost #-opt-prefetch -opt-prefetch-distance=16 #-vec-report6  
  else
    ### compiler specific flags for gcc
    ifeq (gcc,$(findstring gcc,$(CC)))
		  CFLAGS += #-ftree-vectorize -funroll-loops -fprefetch-loop-arrays --param simultaneous-prefetches=4 #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction #-fprofile-generate
        ifeq ($(UNAME), Darwin)
          CFLAGS += -Wa,-q
        endif
    endif

    ### compiler specific flags for clang
    ifeq (clang,$(findstring clang,$(CC)))
		  CFLAGS += -funroll-loops
    endif

    #### common options for gcc and clang
    CFLAGS  += -march=native -mavx -mpopcnt
	  CFLAGS  += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal 
    CFLAGS  +=  -Wcast-align -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
    CLINK += -lm
  endif #compiler is not icc -> which I assume means gcc or clang

  ifeq (USE_MKL,$(findstring USE_MKL,$(OPT)))
	  BLAS_INCLUDE:=-DMKL_ILP64 -m64 -I$(MKLROOT)/include 
    ##Use the Intel MKL library. Check the compiler + openmp
	  ifneq (USE_OMP,$(findstring USE_OMP,$(OPT)))
      ##Link+include sequential libraries
		  ifeq (icc,$(findstring icc,$(CC)))
        ##icc with Intel MKL
			  BLAS_LINK:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
		  else
	      ##gcc with Intel MKL
			  BLAS_LINK:= -Wl,--no-as-needed -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
		  endif
	  else
		  ifeq (icc,$(findstring icc,$(CC)))
        ##icc with Intel MKL+OpenMP
			  BLAS_LINK:= -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
		  else
	      ##gcc with Intel MKL
			  BLAS_LINK:= -Wl,--no-as-needed -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
		  endif
    endif
  else
    ##Use some OpenMP parallel BLAS library (OpenBlas/ATLAS/GSLCBLAS, for instance)
    BLAS_INCLUDE:=$(shell gsl-config --cflags)
    BLAS_LINK:=$(shell gsl-config --libs)
  endif
endif # do checks
