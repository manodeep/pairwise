#pragma once

#define ALIGNMENT   64

#define NELEMENTS 1000
#define MAXLEN    200
#define NDIM      3

#define DOUBLE_PREC
#include "avx_calls.h"

static const int clustered_data = 1;
static const char source_galaxy_file[MAXLEN] = "./data/Mr19_random_subsample_100000.txt";
static const int max_galaxy_in_source_file = 100000;
static const char binfile[MAXLEN] = "./data/bins";



