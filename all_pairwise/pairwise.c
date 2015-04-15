#include <stdio.h>
#include <stdlib.h>
#if defined (__GNUC__) && ! defined (__INTEL_COMPILER)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif

#include <sys/time.h>
#include <assert.h>
#include <math.h>

#if  __INTEL_COMPILER
#include "mkl.h"
#else
#include <gsl/gsl_cblas.h>
#endif

#include "defs.h"
#include "utils.h"

long naive(const double * restrict x1, const double * restrict y1, const int N, double * restrict d)
{
	long numcomputed=0;

	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			double *dist = d + i*N + j;
			const double *y0 = (const double *) &y1[j*NDIM];
			const double *x0 = (const double *) &x1[i*NDIM];
			const double dx = x0[0] - y0[0];
			const double dy = x0[1] - y0[1];
			const double dz = x0[2] - y0[2];
#ifndef SQRT_DIST
			*dist = dx*dx + dy*dy + dz*dz;
#else
			*dist = sqrt(dx*dx + dy*dy + dz*dz);
#endif
		}
 		numcomputed+=N;
	}
	
	return numcomputed;
}

long chunked(const double * restrict x1, const double * restrict y1, const int N, double * restrict d)
{
	long numcomputed = 0;
	const int block_size = 16;
	for(int i=0;i<N;i+=block_size) {
		const int block_size1 = (N-i) > block_size ? block_size:(N-i);
		for(int j=0;j<N;j+=block_size) {
			const int block_size2 = (N-j) > block_size ? block_size:(N-j);
			for(int ii=0;ii<block_size1;ii++) {
				const long index = (i+ii)*(long) N + j;
				double *dist = (double *) &d[index];
				for(int jj=0;jj<block_size2;jj++) {
					const double *y0 = (const double *) (&y1[(j+jj)*NDIM]);
					const double *x0 = (const double *) (&x1[(i+ii)*NDIM]);
					const double dx = x0[0] - y0[0];
					const double dy = x0[1] - y0[1];
					const double dz = x0[2] - y0[2];
					const double sqr_ds = dx*dx + dy*dy + dz*dz;
#ifndef SQRT_DIST
					dist[jj] = sqr_ds;
#else
					dist[jj] = sqrt(sqr_ds);
#endif					
				}
				numcomputed += block_size2;
			}
		}
	}

	return numcomputed;
}


long compiler_vectorized_chunked(const double * restrict x1, const double * restrict y1, const int N, double * restrict d)
{
	long numcomputed = 0;
	const int block_size = 4;
	for(int i=0;i<N;i++) {
		int j;
		for(j=0;j<=(N-block_size);j+=block_size){
			double *dist = &d[i*N + j];
			for(int jj=0;jj<block_size;jj++) {
				const double *y0 = (const double *) (&y1[(j+jj)*NDIM]);
				const double *x0 = (const double *) (&x1[i*NDIM]);
				const double dx = x0[0] - y0[0];
				const double dy = x0[1] - y0[1];
				const double dz = x0[2] - y0[2];
				dist[jj] = dx*dx + dy*dy + dz*dz;
			}
			numcomputed+=block_size;
		}			
		for(int jj=0;j<N;jj++,j++) {
			double *dist = (double *) &d[i*N + j+jj];
			const double *y0 = (const double *) (&y1[(j+jj)*NDIM]);
			const double *x0 = (const double *) (&x1[i*NDIM]);
			const double dx = x0[0] - y0[0];
			const double dy = x0[1] - y0[1];
			const double dz = x0[2] - y0[2];
			dist[jj] = dx*dx + dy*dy + dz*dz;
			numcomputed++;
		}
	}

#ifdef SQRT_DIST
	for(long index=0;index<numcomputed;index++){
		d[index] = sqrt(d[index]);
	}
#endif
	
	return numcomputed;
}


long blas_computed(const double * restrict x1, const double * restrict y1, const int N, double * restrict d)
{
	long numcomputed = (long) N * (long) N;
	const double alpha = -2.0;
	/* struct timeval t0,t1; */
	/* gettimeofday(&t0,NULL); */
	cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasTrans,
							N, N, NDIM,
							alpha, x1, NDIM , y1, NDIM, 0.0, d, N);
	/* gettimeofday(&t1,NULL); */
	/* printf(" \t\t cblas_time = %6.2lf ms\n",1e3*(t1.tv_sec-t0.tv_sec) + 1e-3*(t1.tv_usec - t0.tv_usec)); */
	double *x0 = (double *) x1;
	for(int i=0;i<N;i++){
		const double sqr_x0_norm = x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2];
		for(int j=0;j<N;j++) {
			const double *y0 = (const double *) &y1[j*NDIM];
			const double sqr_y0_norm = y0[0]*y0[0] + y0[1]*y0[1] + y0[2]*y0[2];
			double *dist = (double *) &d[i*N + j];
			const double sum = *dist + (sqr_x0_norm + sqr_y0_norm);
#ifdef SQRT_DIST
			*dist = sqrt(sum);
#else
			*dist = sum;
#endif			
		}
		x0 += NDIM;
	}
	return numcomputed;
}


long check_result(double *dist, const double *x, const double *y, const int N)
{
	int bad=0;
	int numbadprinted=0;
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			const double *tmp_x = (const double *) &(x[i*NDIM]);
			const double *tmp_y = (const double *) &(y[j*NDIM]);
			double sqr_ds = 0.0;
			for(int k=0;k<NDIM;k++) {
				const double dx = tmp_x[k] - tmp_y[k];
				sqr_ds += dx*dx;
			}
			const long index = i*(long) N + j;
			const double this_dist = dist[index];

#ifdef SQRT_DIST
			const double d = sqrt(sqr_ds);
#else
			const double d = sqr_ds;
#endif			
			
			if(fabs(d - this_dist) > 1e-8) {
				if(numbadprinted < 20) {
					fprintf(stderr,"i = %d j = %d d = %0.8lf this_dist = %0.8lf \n",i,j,d,this_dist);
					numbadprinted++;
				}
				bad++;
			}
		}
	}
	
	return bad;
}


int main(int argc, char **argv)
{
	const size_t numbytes = NDIM*NELEMENTS;
	double *x = calloc(sizeof(*x), numbytes);
	double *y = calloc(sizeof(*y), numbytes);
	const long totnpairs = (long) NELEMENTS * (long) NELEMENTS;
	double * restrict dist = calloc(sizeof(*dist), totnpairs);
	assert(x != NULL && y != NULL && dist != NULL && "memory allocation failed");
  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized_chunked","blas_computed","naive","chunked","compiler_vectorized_chunked","blas_computed"};
	const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);
	long (*allfunctions[]) (const double * restrict, const double * restrict, const int, double * restrict)      = {naive, chunked,compiler_vectorized_chunked,blas_computed,
																																																									naive, chunked,compiler_vectorized_chunked,blas_computed};
	const unsigned int seed = 42;
	srand(seed);

	int test_to_run = -1;
	if(argc > 1) {
		test_to_run = atoi(argv[1]);
		if(test_to_run >= ntests) test_to_run = -1;
	}

	if(clustered_data == 0) {
		fill_array(x, NELEMENTS);
		fill_array(y, NELEMENTS);
	} else {
		assert(NELEMENTS <= max_galaxy_in_source_file && "Clustered data does not contain enough galaxies..please reduce NELEMENTS in defs.h");
		assert(NDIM == 3 && "Clustered galaxy data contains *exactly* 3 spatial dimensions");
		read_ascii(x, NELEMENTS, source_galaxy_file);
		read_ascii(y, NELEMENTS, source_galaxy_file);
	}

	printf("#######################################################\n");
	printf("##  Function                            Time (ms)      \n");
	printf("#######################################################\n");
	for(int i=0;i<ntests;i++) {
		struct timeval t0,t1;
		double sum_x=0.0, sum_sqr_x=0.0;
		if(i == test_to_run || test_to_run == -1) {
			//warm-up
			gettimeofday(&t0,NULL);
			long ncomputed = (allfunctions[i]) (x, y, NELEMENTS, dist);
			gettimeofday(&t1,NULL);
			
			//check-result
			long numbad = check_result(dist, x, y, NELEMENTS);
			if(numbad != 0 || ncomputed != totnpairs) {
				fprintf(stderr,"ERROR: Number of incorrectly calculated distances = %ld out of a total of (%ld) possible pairs. Npairs computed = %ld\n", numbad, totnpairs, ncomputed);
				goto cleanup;
			}
#ifndef SQRT_DIST
			const long totflop = totnpairs * (8);
#else
			const long totflop = totnpairs * (8 + 3);
#endif		


			int actual_niterations=0;
			double best_time_in_ms=1e16, best_time_in_megacycles=1e16;
			uint64_t start_cycles, end_cycles;
			for(int iter=0;iter<max_niterations;iter++) {
				actual_niterations++;
				gettimeofday(&t0,NULL);
				start_cycles = rdtsc();
				(allfunctions[i]) (x, y, NELEMENTS, dist);
				end_cycles = rdtsc();
				gettimeofday(&t1,NULL);
			
				const double this_time_in_ms = 1e3*(t1.tv_sec-t0.tv_sec) +  1e-3*(t1.tv_usec - t0.tv_usec);
				if(this_time_in_ms < best_time_in_ms) best_time_in_ms = this_time_in_ms;
				const double this_time_in_megacycles = (end_cycles - start_cycles)/(1024.*1024.);
				if(this_time_in_megacycles < best_time_in_megacycles) best_time_in_megacycles = this_time_in_megacycles;

				sum_sqr_x += this_time_in_ms*this_time_in_ms;
				sum_x     += this_time_in_ms;

				if(max_niterations <= 10) {
					printf("     %-35s  %0.2lf \n",allfunction_names[i], this_time_in_ms);
				}

				const double mean_time  = sum_x/(iter+1.0);
				const double sigma_time = sqrt(sum_sqr_x/(iter+1.0) - mean_time*mean_time);
				//If the std.dev is small compared to typical runtime and
				//the code has run for more than XXX milli-seconds, then break
				if(sigma_time/mean_time < 0.05 && sum_x >= 300.0 && iter >= 10) {
					break;
				}
			}
			
			const double mean_time = sum_x/actual_niterations;
			const double sigma_time = sqrt(sum_sqr_x/actual_niterations - mean_time*mean_time);
			const double gflops = (double) totflop/(1e-3*mean_time)/1e9;
			printf(ANSI_COLOR_RED "# %-35s  %6.2lf +- %5.2lf " ANSI_COLOR_GREEN " (best -- %6.2lf ms, %6.2lf Mcycles) " ANSI_COLOR_RESET "," ANSI_COLOR_BLUE " >= %5.2lf GFlops [%04d iterations]" ANSI_COLOR_RESET "\n",
						 allfunction_names[i], mean_time, sigma_time, best_time_in_ms, best_time_in_megacycles, gflops, actual_niterations);
		}
	}
	

cleanup:
	{
		free((void *) x);free((void *) y);free((void *) dist);
		exit(EXIT_FAILURE);
	}

	exit(EXIT_SUCCESS);
}
