#pragma once

#define ALIGNMENT   64

#define NELEMENTS 1000
#define MAXLEN    200
#define NDIM      3

#define DOUBLE_PREC
#include "function_precision.h"

#ifdef __AVX__
#include "avx_calls.h"
#include "sse_calls.h"
#endif

#ifdef __SSE_4_2__
#include "sse_calls.h"
#include 

#endif
static const int clustered_data = 1;
static const char source_galaxy_file[] = "./data/Mr19_random_subsample_100000.txt";
static const int max_galaxy_in_source_file = 100000;
static const char binfile[] = "./data/bins";
static const int max_niterations = 1000;
static const unsigned int seed = 42;
static const int repeat=5;
static const double convergence_ratio = 0.05;//sigma_runtime/avg_runtime < convergence_ratio to stop repeating the test


