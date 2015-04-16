#include <stdio.h>

#define __USE_XOPEN2K
#include <stdlib.h>
#if defined (__GNUC__) && ! defined (__INTEL_COMPILER)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif

#include <sys/time.h>
#include <assert.h>
#include <math.h>

#include "defs.h"
#include "utils.h"
#include "progressbar.h"

#include "pairwise_3d_ispc.h"

long check_result(double *dist, const double *pos0, const double *pos1, const int N)
{
	long bad=0;
	int numbadprinted=0;
	const double *x0  __attribute__((aligned(ALIGNMENT))) = (const double *) pos0;
	const double *y0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[N];
	const double *z0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[2*N];

	const double *x1  __attribute__((aligned(ALIGNMENT))) = (const double *) pos1;
	const double *y1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[N];
	const double *z1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[2*N];
	
	for(int i=0;i<N;i++) {
		const double xpos  __attribute__((aligned(ALIGNMENT))) = x0[i];
		const double ypos  __attribute__((aligned(ALIGNMENT))) = y0[i];
		const double zpos  __attribute__((aligned(ALIGNMENT))) = z0[i];
		
		for(int j=0;j<N;j++) {
			const double dx = xpos - x1[j];
			const double dy = ypos - y1[j];
			const double dz = zpos - z1[j];
			const double sqr_ds = dx*dx + dy*dy + dz*dz;
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


void naive(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d)
{
	const double *x0  __attribute__((aligned(ALIGNMENT))) = (const double *) pos0;
	const double *y0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[N];
	const double *z0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[2*N];

	for(int i=0;i<N;i++) {
		const double xpos  __attribute__((aligned(ALIGNMENT))) = x0[i];
		const double ypos  __attribute__((aligned(ALIGNMENT))) = y0[i];
		const double zpos  __attribute__((aligned(ALIGNMENT))) = z0[i];
		const double *x1  __attribute__((aligned(ALIGNMENT))) = (const double *) pos1;
		const double *y1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[N];
		const double *z1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[2*N];

		double *dist = (double *) &d[i*N];	
		for(int j=0;j<N;j++) {
			const double dx = xpos - x1[j];
			const double dy = ypos - y1[j];
			const double dz = zpos - z1[j];
#ifndef SQRT_DIST
			*dist = dx*dx + dy*dy + dz*dz;
#else
			*dist = sqrt(dx*dx + dy*dy + dz*dz);
#endif
			dist++;
		}
	}
}



void chunked(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d)
{
	const int block_size = 128;
	const double *x0 = (const double *) pos0;
	const double *y0 = (const double *) &pos0[N];
	const double *z0 = (const double *) &pos0[2*N];
	
	for(int i=0;i<N;i+=block_size) {
		const int block_size1 = (N-i) > block_size ? block_size:(N-i);
		for(int j=0;j<N;j+=block_size) {
			const int block_size2 = (N-j) > block_size ? block_size:(N-j);
			const double *x1 = (const double *) &pos1[j];
			const double *y1 = (const double *) &pos1[N + j];
			const double *z1 = (const double *) &pos1[2*N + j];
			for(int ii=0;ii<block_size1;ii++) {
				const long index = (i+ii)*(long) N + j;
				double *dist = (double *) &d[index];
				const double xpos = x0[i+ii];
				const double ypos = y0[i+ii];
				const double zpos = z0[i+ii];
#ifdef __INTEL_COMPILER				
#pragma simd assert
#endif				
				for(int jj=0;jj<block_size2;jj++) {
					const double dx = xpos - x1[jj];
					const double dy = ypos - y1[jj];
					const double dz = zpos - z1[jj];
					const double sqr_ds = dx*dx + dy*dy + dz*dz;
#ifndef SQRT_DIST
					dist[jj] = sqr_ds;
#else
					dist[jj] = sqrt(sqr_ds);
#endif					
				}
			}
		}
	}
}


void compiler_vectorized_chunked(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d)
{
	const int block_size = 4;

	const double *x0 = (const double *) pos0;
	const double *y0 = (const double *) &pos0[N];
	const double *z0 = (const double *) &pos0[2*N];


	for(int i=0;i<N;i++) {
		int j;
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];

		double *x1 = (double *) pos1;
		double *y1 = (double *) &pos1[N];
		double *z1 = (double *) &pos1[2*N];


		for(j=0;j<=(N-block_size);j+=block_size){
			double *dist = &d[i*N + j];
#ifdef __INTEL_COMPILER			
#pragma simd assert
#endif
			for(int jj=0;jj<block_size;jj++) {
				const double dx = xpos - x1[jj];
				const double dy = ypos - y1[jj];
				const double dz = zpos - z1[jj];
				dist[jj] = dx*dx + dy*dy + dz*dz;
			}
			x1 += block_size;
			y1 += block_size;
			z1 += block_size;
			
		}			
		for(int jj=0;j<N;jj++,j++) {
			double *dist = (double *) &d[i*N + j+jj];
			const double dx = xpos - x1[jj];
			const double dy = ypos - y1[jj];
			const double dz = zpos - z1[jj];
			dist[jj] = dx*dx + dy*dy + dz*dz;
		}
	}

#ifdef SQRT_DIST
	long numcomputed = (long) (N) *(long) N;
	for(long index=0;index<numcomputed;index++){
		d[index] = sqrt(d[index]);
	}
#endif
	
}


void intrinsics_chunked(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d)
{
	const int block_size = 4;
	double *x0 = (double *) pos0;
	double *y0 = (double *) &pos0[N];
	double *z0 = (double *) &pos0[2*N];

	for(int i=0;i<N;i++) {

		const double xpos = *x0;
		const double ypos = *y0;
		const double zpos = *z0;
		
		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);

		x0++;y0++;z0++;

		double *x1 = (double *) pos1;
		double *y1 = (double *) &pos1[N];
		double *z1 = (double *) &pos1[2*N];

		int j;		


		for(j=0;j<=(N-block_size);j+=block_size){
			/* PREFETCH(x1 + 128); */
			/* PREFETCH(y1 + 128); */
			/* PREFETCH(z1 + 128); */
			
			double *dist  __attribute__((aligned(ALIGNMENT))) = &d[i*N + j];
			/* const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(x1); */
			/* const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(y1); */
			/* const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(z1); */

			const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_ALIGNED(x1);
			const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_ALIGNED(y1);
			const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_ALIGNED(z1);

			
			x1 += block_size;
			y1 += block_size;
			z1 += block_size;

			const AVX_FLOATS m_dx = AVX_SUBTRACT_FLOATS(m_x1, m_xpos);
			const AVX_FLOATS m_dy = AVX_SUBTRACT_FLOATS(m_y1, m_ypos);
			const AVX_FLOATS m_dz = AVX_SUBTRACT_FLOATS(m_z1, m_zpos);
			
			const AVX_FLOATS m_r2 = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dx), AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dy), AVX_SQUARE_FLOAT(m_dz)));
#ifdef SQRT_DIST
			/* AVX_STORE_FLOATS_TO_MEMORY(dist, AVX_SQRT_FLOAT(m_r2)); */
			AVX_STREAMING_STORE_FLOATS(dist, AVX_SQRT_FLOAT(m_r2));
#else
			/* AVX_STORE_FLOATS_TO_MEMORY(dist, m_r2); */
			AVX_STREAMING_STORE_FLOATS(dist, m_r2);
#endif			
		}			
		for(int jj=0;j<N;jj++,j++) {
			double *dist = (double *) &d[i*N + j+jj];
			const double dx = xpos - x1[jj];
			const double dy = ypos - y1[jj];
			const double dz = zpos - z1[jj];
#ifdef SQRT_DIST
			dist[jj] = sqrt(dx*dx + dy*dy + dz*dz);
#else
			dist[jj] = dx*dx + dy*dy + dz*dz;
#endif
		}
	}
}



void intrinsics_chunked_unroll(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d)
{
	const int block_size = 4;
#ifdef SQRT_DIST
	const int unroll_factor=4;//28 sqrt operations are supported -> can be up to 28
#else
	const int unroll_factor=2;
#endif	
	double *x0 = (double *) pos0;
	double *y0 = (double *) &pos0[N];
	double *z0 = (double *) &pos0[2*N];

#define GETXVALUE(n) AVX_FLOATS diffx##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&x1[n*block_size]), m_xpos)
#define GETYVALUE(n) AVX_FLOATS diffy##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&y1[n*block_size]), m_ypos)
#define GETZVALUE(n) AVX_FLOATS diffz##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&z1[n*block_size]), m_zpos)
#define GETVALUE(n)  GETXVALUE(n);GETYVALUE(n);GETZVALUE(n)

#ifndef SQRT_DIST	
#define GETDIST(n)  AVX_STREAMING_STORE_FLOATS(&dist[n*block_size],AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n))))
#else	
#define GETDIST(n) AVX_STREAMING_STORE_FLOATS(&dist[n*block_size],AVX_SQRT_FLOAT(AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n)))))
#endif
	
	for(int i=0;i<N;i++) {

		const double xpos = *x0;
		const double ypos = *y0;
		const double zpos = *z0;
		
		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);

		x0++;y0++;z0++;

		double *x1 = (double *) pos1;
		double *y1 = (double *) &pos1[N];
		double *z1 = (double *) &pos1[2*N];

		double *dist = (double *) &d[i*N];
		
		int offset = 0;
		const size_t bytes_offset = block_size*unroll_factor;
		for(int j=0; j <=(N/block_size-unroll_factor);j+=unroll_factor) {
			/* PREFETCH(x1 + 2*block_size*unroll_factor); */
			/* PREFETCH(y1 + 2*block_size*unroll_factor); */
			/* PREFETCH(z1 + 2*block_size*unroll_factor); */

			
#ifdef SQRT_DIST			
			GETVALUE(0);GETVALUE(1);GETVALUE(2);GETVALUE(3);
			GETDIST(0);GETDIST(1);GETDIST(2);GETDIST(3);
#else
			GETVALUE(0);GETVALUE(1);
			GETDIST(0);GETDIST(1);
#endif			
			x1 += bytes_offset;
			y1 += bytes_offset;
			z1 += bytes_offset;
			dist += bytes_offset;
			offset += bytes_offset;
		}

		dist = (double *) &d[i*N + offset];
		for(int jj=0;offset<N;jj++,offset++) {
			const double dx = xpos - *x1;
			const double dy = ypos - *y1;
			const double dz = zpos - *z1;

			x1++;y1++;z1++;
			
#ifdef SQRT_DIST
			dist[jj] = sqrt(dx*dx + dy*dy + dz*dz);
#else
			dist[jj] = dx*dx + dy*dy + dz*dz;
#endif
		}
	}

#undef GETXVALUE
#undef GETYVALUE
#undef GETZVALUE
#undef GETDIST
	
}

int main(void)
{

#if (NDIM != 3)
	#error NDIM must be set to 3
#endif
	
	const size_t numbytes = NDIM*NELEMENTS;
	double *x    __attribute__((aligned(ALIGNMENT))) = NULL;
	double *y    __attribute__((aligned(ALIGNMENT))) = NULL;
	double *dist __attribute__((aligned(ALIGNMENT))) = NULL;
	int test0 = posix_memalign((void **) &x, ALIGNMENT, sizeof(*x)*numbytes);
	int test1 = posix_memalign((void **) &y, ALIGNMENT, sizeof(*y)*numbytes); 
	
	const long totnpairs = (long) NELEMENTS * (long) NELEMENTS;
	int test2 = posix_memalign((void **) &dist, ALIGNMENT, sizeof(*dist)*totnpairs);
	assert(test0 == 0  && test1 == 0 && test2 == 0 && "memory allocation failed");

  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized_chunked","intrinsics_chunked","intrinsics_chunked_unroll","pairwise_ispc"};
	const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);
	void (*allfunctions[]) (const double * restrict x, const double * restrict y, const int, double * restrict)      = {naive,chunked,compiler_vectorized_chunked,intrinsics_chunked,intrinsics_chunked_unroll,pairwise_ispc};

	double function_best_mean_time[ntests],function_sigma_time[ntests],function_best_time_in_ms[ntests],function_best_mcycles[ntests];
	int function_niterations[ntests];
	for(int i=0;i<ntests;i++) {
		function_best_mean_time[i] = 1e16;
		function_sigma_time[i] = 0.0;
		function_best_time_in_ms[i] = 1e16;
		function_best_mcycles[i] = 1e16;
	}
	
#ifndef SQRT_DIST
	const long totflop = (long) NELEMENTS * (long) NELEMENTS * (8);
#else
	const long totflop = (long) NELEMENTS * (long) NELEMENTS * (8 + 10); //some hand-wavy thing saying that sqrt is 10 flop
#endif		

	const unsigned int seed = 42;
	srand(seed);

	if(clustered_data == 0) {
		fill_array(x, NELEMENTS);
		fill_array(y, NELEMENTS);
	} else {
		assert(NELEMENTS <= max_galaxy_in_source_file && "Clustered data does not contain enough galaxies..please reduce NELEMENTS in defs.h");
		assert(NDIM == 3 && "Clustered galaxy data contains *exactly* 3 spatial dimensions");
		read_ascii(x, NELEMENTS, source_galaxy_file);
		read_ascii(y, NELEMENTS, source_galaxy_file);
	}

	const int repeat=5;
	const int64_t totniterations = repeat*ntests*(int64_t) max_niterations;
	int64_t numdone = 0;
	int interrupted=0;

	printf("Running benchmarks...\n");
	init_my_progressbar(totniterations, &interrupted);

	for(int irep=0;irep<repeat;irep++) {
		for(int i=0;i<ntests;i++) {
			struct timeval t0,t1;
			double sum_x=0.0, sum_sqr_x=0.0;
			
			//warm-up
			(allfunctions[i]) (x, y, NELEMENTS, dist);
			
			//check-result
			long numbad = check_result(dist, x, y, NELEMENTS);
			if(numbad != 0) {
				fprintf(stderr,"ERROR: Number of incorrectly calculated distances = %ld out of a total of (%ld) possible pairs.\n", numbad, totnpairs);
				goto cleanup;
			}
			
			double best_time_in_ms=1e16, best_time_in_megacycles=1e16;
			uint64_t start_cycles, end_cycles;
			const int64_t numdone_before_iter_loop = numdone;
			for(int iter=0;iter<max_niterations;iter++) {
				gettimeofday(&t0,NULL);
				start_cycles = rdtsc();
				(allfunctions[i]) (x, y, NELEMENTS, dist);
				end_cycles = rdtsc();
				gettimeofday(&t1,NULL);
				
				numdone++;
				my_progressbar(numdone,&interrupted);
				
				const double this_time_in_ms = 1e3*(t1.tv_sec-t0.tv_sec) +  1e-3*(t1.tv_usec - t0.tv_usec);
				if(this_time_in_ms < best_time_in_ms) best_time_in_ms = this_time_in_ms;
				const double this_time_in_megacycles = (end_cycles - start_cycles)/(1024.*1024.);
				if(this_time_in_megacycles < best_time_in_megacycles) best_time_in_megacycles = this_time_in_megacycles;
				
				sum_sqr_x += this_time_in_ms*this_time_in_ms;
				sum_x     += this_time_in_ms;
				
				if(max_niterations <= 10) {
					printf("     %-35s  %0.2lf \n",allfunction_names[i], this_time_in_ms);
				}
				
				if(best_time_in_ms < function_best_time_in_ms[i]) {
					function_best_time_in_ms[i] = best_time_in_ms;
				}
				
				if(best_time_in_megacycles < function_best_mcycles[i]) {
					function_best_mcycles[i] = best_time_in_megacycles;
				}
				
				const double mean_time  = sum_x/(iter+1.0);
				const double sigma_time = sqrt(sum_sqr_x/(iter+1.0) - mean_time*mean_time);
				//If the std.dev is small compared to typical runtime and
				//the code has run for more than XXX milli-seconds, then break
				if(sigma_time/mean_time < 0.05 && sum_x >= 300.0 && iter >= 10) {
					if(mean_time < function_best_mean_time[i]) {
						function_best_mean_time[i] = mean_time;
						function_sigma_time[i] = sigma_time;
						function_niterations[i] = iter + 1;
					}
					numdone = numdone_before_iter_loop + max_niterations;
					break;
				}
			}
		}//i loop over ntests
	}//irep loop over nrepeat
	finish_myprogressbar(&interrupted);
	
	printf("##################################################\n");
	printf("##  Function                            Time (ms) \n");
	printf("##################################################\n");
	for(int i=0;i<ntests;i++) {
		const double mean_time = function_best_mean_time[i];
		const double sigma_time = function_sigma_time[i];
		const double gflops = (double) totflop/(1e-3*mean_time)/1e9;
		const double best_time_in_ms = function_best_time_in_ms[i];
		const double best_time_in_megacycles = function_best_mcycles[i];
		const int actual_niterations = function_niterations[i];
		printf(ANSI_COLOR_RED "# %-35s  %6.2lf +- %5.2lf " ANSI_COLOR_GREEN " (best -- %6.2lf ms, %6.2lf Mcycles) " ANSI_COLOR_RESET "," ANSI_COLOR_BLUE " %5.2lf GFlops [%04d iterations]" ANSI_COLOR_RESET "\n",
					 allfunction_names[i], mean_time, sigma_time, best_time_in_ms, best_time_in_megacycles, gflops, actual_niterations);

	}

	
cleanup:
	{
		free((void *) x);free((void *) y);free((void *) dist);
		exit(EXIT_FAILURE);
	}

	exit(EXIT_SUCCESS);
}
