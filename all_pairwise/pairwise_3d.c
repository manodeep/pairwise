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

#ifdef ISPC_AVAIL
#include "pairwise_3d_ispc.h"
#endif

#if (NDIM != 3)
	#error NDIM must be set to 3
#endif


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


void naive(const double * restrict x0, const double * restrict y0, const double * restrict z0,
					 const double * restrict x1, const double * restrict y1, const double * restrict z1,
					 const int N0, const int N1, 
					 double * restrict d)
{

	for(int i=0;i<N0;i++) {
		const double xpos  __attribute__((aligned(ALIGNMENT))) = x0[i];
		const double ypos  __attribute__((aligned(ALIGNMENT))) = y0[i];
		const double zpos  __attribute__((aligned(ALIGNMENT))) = z0[i];

		double *dist = (double *) &d[i*N1];	
		for(int j=0;j<N1;j++) {
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



void chunked(const double * restrict x0, const double * restrict y0, const double * restrict z0,
						 const double * restrict x1, const double * restrict y1, const double * restrict z1,
						 const int N0, const int N1, 
						 double * restrict d)
{
	const int block_size = 256;
	for(int i=0;i<N0;i+=block_size) {
		const int block_size1 = (N0-i) > block_size ? block_size:(N0-i);
		for(int j=0;j<N1;j+=block_size) {
			const int block_size2 = (N1-j) > block_size ? block_size:(N1-j);
			for(int ii=0;ii<block_size1;ii++) {
				const long index = (i+ii)*(long) N1 + j;
				double *dist = (double *) &d[index];
				const double xpos = x0[i+ii];
				const double ypos = y0[i+ii];
				const double zpos = z0[i+ii];
#ifdef __INTEL_COMPILER				
#pragma simd assert
#endif				
				for(int jj=0;jj<block_size2;jj++) {
					const double dx = xpos - x1[j+jj];
					const double dy = ypos - y1[j+jj];
					const double dz = zpos - z1[j+jj];
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


void compiler_vectorized(const double * restrict x0, const double * restrict y0, const double * restrict z0,
												 const double * restrict x1, const double * restrict y1, const double * restrict z1,
												 const int N0, const int N1, 
												 double * restrict d)
{
	for(int i=0;i<N0;i++) {
		int j;
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;


		for(j=0;j<=(N1-NVEC);j+=NVEC){
			double *dist = &d[i*N1 + j];
#ifdef __INTEL_COMPILER			
#pragma simd assert
#endif
			for(int jj=0;jj<NVEC;jj++) {
				const double dx = xpos - localx1[jj];
				const double dy = ypos - localy1[jj];
				const double dz = zpos - localz1[jj];
				dist[jj] = dx*dx + dy*dy + dz*dz;
			}
			localx1 += NVEC;
			localy1 += NVEC;
			localz1 += NVEC;
			
		}			
		double *dist = (double *) &d[i*N1 + j];
		for(int jj=0;j<N1;jj++,j++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
			dist[jj] = dx*dx + dy*dy + dz*dz;
		}
	}

#ifdef SQRT_DIST
	long numcomputed = (long) N0 *(long) N1;
	for(long index=0;index<numcomputed;index++){
		d[index] = sqrt(d[index]);
	}
#endif
	
}


#ifdef __AVX__
void avx_intrinsics(const double * restrict x0, const double * restrict y0, const double * restrict z0,
										const double * restrict x1, const double * restrict y1, const double * restrict z1,
										const int N0, const int N1, 
										double * restrict d)
{
	for(int i=0;i<N0;i++) {
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];

		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		int j;		
		for(j=0;j<=(N1-NVEC);j+=NVEC){
			double *dist  __attribute__((aligned(ALIGNMENT))) = &d[i*N1 + j];
			const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(localx1);
			const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(localy1);
			const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(localz1);
			
			localx1 += NVEC;localy1 += NVEC;localz1 += NVEC;

			const AVX_FLOATS m_dx = AVX_SUBTRACT_FLOATS(m_x1, m_xpos);
			const AVX_FLOATS m_dy = AVX_SUBTRACT_FLOATS(m_y1, m_ypos);
			const AVX_FLOATS m_dz = AVX_SUBTRACT_FLOATS(m_z1, m_zpos);
			
			const AVX_FLOATS m_r2 = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dx), AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dy), AVX_SQUARE_FLOAT(m_dz)));
#ifdef SQRT_DIST
			AVX_STORE_FLOATS_TO_MEMORY(dist, AVX_SQRT_FLOAT(m_r2));
			/* AVX_STREAMING_STORE_FLOATS(dist, AVX_SQRT_FLOAT(m_r2)); */
#else
			AVX_STORE_FLOATS_TO_MEMORY(dist, m_r2);
			/* AVX_STREAMING_STORE_FLOATS(dist, m_r2); */
#endif			
		}			

		double *dist = (double *) &d[i*N1 + j];
		for(int jj=0;j<N1;jj++,j++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
#ifdef SQRT_DIST
			dist[jj] = sqrt(dx*dx + dy*dy + dz*dz);
#else
			dist[jj] = dx*dx + dy*dy + dz*dz;
#endif
		}
	}
}



void avx_intrinsics_unroll(const double * restrict x0, const double * restrict y0, const double * restrict z0,
													 const double * restrict x1, const double * restrict y1, const double * restrict z1,
													 const int N0, const int N1, 
													 double * restrict d)
{
#ifdef SQRT_DIST
	const int unroll_factor=4;//28 sqrt operations are supported -> can be up to 28
#else
	const int unroll_factor=2;
#endif
	
#define GETXVALUE(n) AVX_FLOATS diffx##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localx1[n*NVEC]), m_xpos)
#define GETYVALUE(n) AVX_FLOATS diffy##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localy1[n*NVEC]), m_ypos)
#define GETZVALUE(n) AVX_FLOATS diffz##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localz1[n*NVEC]), m_zpos)
#define GETVALUE(n)  GETXVALUE(n);GETYVALUE(n);GETZVALUE(n)

#ifndef SQRT_DIST	
#define GETDIST(n) AVX_STORE_FLOATS_TO_MEMORY(&dist[n*NVEC],AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n))))
#else	
#define GETDIST(n) AVX_STORE_FLOATS_TO_MEMORY(&dist[n*NVEC],AVX_SQRT_FLOAT(AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n)))))
#endif

	for(int i=0;i<N0;i++) {

		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];
		
		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		double *dist = (double *) &d[i*N1];
		
		int offset = 0;
		const size_t bytes_offset = NVEC*unroll_factor;
		for(int j=0; j <=(N1/NVEC-unroll_factor);j+=unroll_factor) {
#ifdef SQRT_DIST			
			GETVALUE(0);GETVALUE(1);GETVALUE(2);GETVALUE(3);
			GETDIST(0);GETDIST(1);GETDIST(2);GETDIST(3);
#else
			GETVALUE(0);GETVALUE(1);
			GETDIST(0);GETDIST(1);
#endif			
			localx1 += bytes_offset;
			localy1 += bytes_offset;
			localz1 += bytes_offset;
			dist += bytes_offset;
			offset += bytes_offset;
		}

		dist = (double *) &d[i*N1 + offset];
		for(int jj=0;offset<N1;jj++,offset++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
			
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


static inline void avx_intrinsics_kernel(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																				 const double * restrict x1, const double * restrict y1, const double * restrict z1,
																				 const int N0, const int N1,
																				 const int totN1,
																				 double * restrict d) 
{
	for(int i=0;i<N0;i++) {
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];

		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		int j;		
		for(j=0;j<=(N1-NVEC);j+=NVEC){
			double *dist  __attribute__((aligned(ALIGNMENT))) = &d[i*totN1 + j];
			const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(localx1);
			const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(localy1);
			const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(localz1);
			
			localx1 += NVEC;localy1 += NVEC;localz1 += NVEC;

			const AVX_FLOATS m_dx = AVX_SUBTRACT_FLOATS(m_x1, m_xpos);
			const AVX_FLOATS m_dy = AVX_SUBTRACT_FLOATS(m_y1, m_ypos);
			const AVX_FLOATS m_dz = AVX_SUBTRACT_FLOATS(m_z1, m_zpos);
			
			const AVX_FLOATS m_r2 = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dx), AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dy), AVX_SQUARE_FLOAT(m_dz)));
#ifdef SQRT_DIST
			AVX_STORE_FLOATS_TO_MEMORY(dist, AVX_SQRT_FLOAT(m_r2));
			/* AVX_STREAMING_STORE_FLOATS(dist, AVX_SQRT_FLOAT(m_r2)); */
#else
			AVX_STORE_FLOATS_TO_MEMORY(dist, m_r2);
			/* AVX_STREAMING_STORE_FLOATS(dist, m_r2); */
#endif			
		}			

		double *dist = (double *) &d[i*totN1 + j];
		for(int jj=0;j<N1;jj++,j++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
#ifdef SQRT_DIST
			dist[jj] = sqrt(dx*dx + dy*dy + dz*dz);
#else
			dist[jj] = dx*dx + dy*dy + dz*dz;
#endif
		}
	}
}

void avx_intrinsics_chunked(const double * restrict x0, const double * restrict y0, const double * restrict z0,
														const double * restrict x1, const double * restrict y1, const double * restrict z1,
														const int N0, const int N1, 
														double * restrict d)
{
	const int block_size = 256;
	for(int i=0;i<N0;i+=block_size) {
		const int block_size1 = (N0-i) > block_size ? block_size:(N0-i);
		for(int j=0;j<N1;j+=block_size) {
			const int block_size2 = (N1-j) > block_size ? block_size:(N1-j);
			avx_intrinsics_kernel(&x0[i],&y0[i],&z0[i],
														&x1[j],&y1[j],&z1[j],
														block_size1,block_size2,
														N1,
														&d[i*N1 + j]);
														
		}
	}
}


#endif //AVX


int main(int argc, char **argv)
{
	int numpart = 0;
	if(argc > 1 ) {
		numpart = atoi(argv[1]);
	} else {
		numpart = NELEMENTS;
	}
	
	const size_t numbytes = NDIM*numpart;
	double *pos0    __attribute__((aligned(ALIGNMENT))) = NULL;
	double *pos1    __attribute__((aligned(ALIGNMENT))) = NULL;
	double *dist __attribute__((aligned(ALIGNMENT))) = NULL;
	int test0 = posix_memalign((void **) &pos0, ALIGNMENT, sizeof(*pos0)*numbytes);
	int test1 = posix_memalign((void **) &pos1, ALIGNMENT, sizeof(*pos1)*numbytes); 
	
	const long totnpairs = (long) numpart * (long) numpart;
	int test2 = posix_memalign((void **) &dist, ALIGNMENT, sizeof(*dist)*totnpairs);
	assert(test0 == 0  && test1 == 0 && test2 == 0 && "memory allocation failed");

  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized"
#ifdef __AVX__																						
																						,"avx_intrinsics","avx_intrinsics_unroll","avx_intrinsics_chunked"
#endif																						
#ifdef ISPC_AVAIL
																						,"pairwise_ispc"
#endif
  };
	const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);

	//Yup. This is the array of function pointers. 
	void (*allfunctions[]) (const double * restrict x0, const double * restrict y0, const double * restrict z0,
													const double * restrict x1, const double * restrict y1, const double * restrict z1,
													const int N0, const int N1, 
													double * restrict d)  = {naive,chunked,compiler_vectorized
#ifdef __AVX__																																																											
																									 ,avx_intrinsics,avx_intrinsics_unroll,avx_intrinsics_chunked
#endif																																																											
#ifdef ISPC_AVAIL
																									 ,pairwise_ispc
#endif                                                   
  };

	//end of block for function pointer array. 

	double function_best_mean_time[ntests],function_sigma_time[ntests],function_best_time_in_ms[ntests],function_best_mcycles[ntests];
	int function_niterations[ntests];
	for(int i=0;i<ntests;i++) {
		function_best_mean_time[i] = 1e16;
		function_sigma_time[i] = 0.0;
		function_best_time_in_ms[i] = 1e16;
		function_best_mcycles[i] = 1e16;
	}
	
#ifndef SQRT_DIST
	const long totflop = (long) numpart * (long) numpart * (8);
#else
	const long totflop = (long) numpart * (long) numpart * (8 + 10); //some hand-wavy thing saying that sqrt is 10 flop
#endif		

	srand(seed);

	if(clustered_data == 0) {
		fill_array(pos0, numpart);
		fill_array(pos1, numpart);
	} else {
		assert(numpart <= max_galaxy_in_source_file && "Clustered data does not contain enough galaxies..please reduce NELEMENTS in defs.h or pass a smaller number on the command-line");
		assert(NDIM == 3 && "Clustered galaxy data contains *exactly* 3 spatial dimensions");
		read_ascii(pos0, numpart, source_galaxy_file);
		read_ascii(pos1, numpart, source_galaxy_file);
	}

	const int64_t totniterations = repeat*ntests*(int64_t) max_niterations;
	int64_t numdone = 0;
	int interrupted=0;

	fprintf(stderr,"# Running benchmarks with N = %05d particles\n",numpart);
	init_my_progressbar(totniterations, &interrupted);

	const double *x0 = pos0;
	const double *y0 = &pos0[numpart];
	const double *z0 = &pos0[2*numpart];

	const double *x1 = pos1;
	const double *y1 = &pos1[numpart];
	const double *z1 = &pos1[2*numpart];
	
	
	for(int irep=0;irep<repeat;irep++) {
		for(int i=0;i<ntests;i++) {
			struct timeval t0,t1;
			double sum_x=0.0, sum_sqr_x=0.0;
			
			//warm-up
			(allfunctions[i]) (x0,y0,z0,x1,y1,z1,numpart,numpart,dist);
			
			//check-result
			long numbad = check_result(dist, pos0, pos1, numpart);
			if(numbad != 0) {
				fprintf(stderr,"ERROR: In function `%s' Number of incorrectly calculated distances = %ld out of a total of (%ld) possible pairs.\n", allfunction_names[i],numbad, totnpairs);
				interrupted=1;
				goto cleanup;
			}
			
			double best_time_in_ms=1e16, best_time_in_megacycles=1e16;
			uint64_t start_cycles, end_cycles;
			const int64_t numdone_before_iter_loop = numdone;
			for(int iter=0;iter<max_niterations;iter++) {
				gettimeofday(&t0,NULL);
				start_cycles = rdtsc();
				(allfunctions[i]) (x0,y0,z0,x1,y1,z1,numpart,numpart,dist);
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
				if(mean_time < function_best_mean_time[i]) {
					function_best_mean_time[i] = mean_time;
				}
				function_niterations[i] = iter + 1;
				function_sigma_time[i] = sigma_time;
				
        //If the std.dev is small compared to typical runtime and
				//the code has run for more than XXX milli-seconds, then break
				if(sigma_time/mean_time < convergence_ratio && sum_x >= 300 && iter >= 5) {
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
		free((void *) pos0);free((void *) pos1);free((void *) dist);
		exit(EXIT_FAILURE);
	}

	exit(EXIT_SUCCESS);
}
