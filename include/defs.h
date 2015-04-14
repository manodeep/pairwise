#pragma once

#define ALIGNMENT   64

#define NELEMENTS 1000
#define MAXLEN    200
#define NDIM      3

#define DOUBLE_PREC
#include "avx_calls.h"

const int clustered_data = 1;
const char source_galaxy_file[] = "./data/Mr19_random_subsample_100000.txt";
const int max_galaxy_in_source_file = 100000;
const char *binfile = "./data/bins";

