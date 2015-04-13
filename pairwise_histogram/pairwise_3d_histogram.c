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
#include <stdint.h>
#include <inttypes.h>
#include "utils.h"

#define NELEMENTS 1000
#define MAXLEN    100

//Even though it seems like the code
//scales to arbitrary dimensions, the individual
//functions *only* calculate in 3-D.

#define NDIM      3

#define DOUBLE_PREC
#include "avx_calls.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define ALIGNMENT   64


void fill_array(double * restrict x, const int N)
{
  const double inv_rand_max = 1.0 / RAND_MAX;
  double *pos = (double *) x;
  for(int j=0;j<NDIM;j++) {
	for(int i=0;i<N;i++) {
	  *pos = rand() * inv_rand_max;
	  pos++;
	}
  }
}

int64_t * calculate_reference_histogram(const double * restrict pos0, const double * restrict pos1, const int N, const char *binfile)
{
  struct timeval t0, t1;
  gettimeofday(&t0, NULL);
  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");

  int64_t *npairs_reference = NULL;
  size_t numbytes = sizeof(*npairs_reference)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs_reference, ALIGNMENT, numbytes);
  memset(npairs_reference, 0, numbytes);

  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
	rupp_sqr[i] = rupp[i]*rupp[i];
  }
	
  double *d = NULL;
  const int64_t totnpairs = (int64_t) N * (int64_t) N;
  numbytes = sizeof(*d) * totnpairs;
  int test1 = posix_memalign((void **) &d, ALIGNMENT, numbytes);
  assert(test0 == 0 && test1 == 0 && "memory allocation failed");

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
	  *dist = dx*dx + dy*dy + dz*dz;
	  dist++;
	}
  }

  for(int64_t i=0;i<totnpairs;i++) {
	const DOUBLE r2 = d[i];
	for(int kbin=nrpbin-1;kbin>=1;kbin--){
	  if(r2 >= rupp_sqr[kbin-1]) {
		npairs_reference[kbin]++;
		break;
	  }
	}//searching for kbin loop
  }
  
  free(d);free(rupp);
  gettimeofday(&t1,NULL);
  const double this_time_in_ms = 1e3*(t1.tv_sec-t0.tv_sec) +  1e-3*(t1.tv_usec - t0.tv_usec);
  printf("Time taken to generate npairs reference implementation = %6.2lf ms\n",this_time_in_ms);
  return npairs_reference;
}


int check_result(const int64_t *npairs, const int64_t *npairs_reference, const int Nbins)
{
  int bad=0;
  int numbadprinted=0;
	
  for(int i=1;i<Nbins;i++) {
	if(npairs[i] != npairs_reference[i]) {
	  if(numbadprinted < 20) {
		fprintf(stderr,"i = %d npairs = %"PRId64" npairs_reference = %"PRId64" \n",i,npairs[i],npairs_reference[i]);
		numbadprinted++;
	  }
	  bad++;
	}
  }
	
  return bad;
}


int64_t * naive(const double * restrict pos0, const double * restrict pos1, const int N, const char *binfile)
{
  const double *x0  __attribute__((aligned(ALIGNMENT))) = (const double *) pos0;
  const double *y0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[N];
  const double *z0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[2*N];

  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
	rupp_sqr[i] = rupp[i]*rupp[i];
  }
	
  for(int i=0;i<N;i++) {
	const double xpos  __attribute__((aligned(ALIGNMENT))) = x0[i];
	const double ypos  __attribute__((aligned(ALIGNMENT))) = y0[i];
	const double zpos  __attribute__((aligned(ALIGNMENT))) = z0[i];
	const double *x1  __attribute__((aligned(ALIGNMENT))) = (const double *) pos1;
	const double *y1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[N];
	const double *z1  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos1[2*N];

	for(int j=0;j<N;j++) {
	  const double dx = xpos - x1[j];
	  const double dy = ypos - y1[j];
	  const double dz = zpos - z1[j];
	  const double r2 = dx*dx + dy*dy + dz*dz;
	  for(int kbin=nrpbin-1;kbin>=1;kbin--){
		if(r2 >= rupp_sqr[kbin-1]) {
		  npairs[kbin]++;
		  break;
		}
	  }//searching for kbin loop
	}
  }

  free(rupp);
  return npairs;
}



int64_t * chunked(const double * restrict pos0, const double * restrict pos1, const int N, const char *binfile)
{
  const int block_size = 64;
  const double *x0 = (const double *) pos0;
  const double *y0 = (const double *) &pos0[N];
  const double *z0 = (const double *) &pos0[2*N];
  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);

  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
	rupp_sqr[i] = rupp[i]*rupp[i];
  }

  for(int i=0;i<N;i+=block_size) {
	const int block_size1 = (N-i) > block_size ? block_size:(N-i);
	for(int j=0;j<N;j+=block_size) {
	  const int block_size2 = (N-j) > block_size ? block_size:(N-j);
	  const double *x1 = (const double *) &pos1[j];
	  const double *y1 = (const double *) &pos1[N + j];
	  const double *z1 = (const double *) &pos1[2*N + j];
	  for(int ii=0;ii<block_size1;ii++) {
		double dist[block_size2];
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
		  dist[jj] = sqr_ds;
		}

		for(int jj=0;jj<block_size2;jj++) {
		  const double r2 = dist[jj];
		  for(int kbin=nrpbin-1;kbin>=1;kbin--){
			if(r2 >= rupp_sqr[kbin-1]) {
			  npairs[kbin]++;
			  break;
			}
		  }//searching for kbin loop
		}
	  }
	}
  }

  free(rupp);
  return npairs;
}


int64_t * compiler_vectorized_chunked(const double * restrict pos0, const double * restrict pos1, const int N, const char *binfile)
{
  const int block_size = 4;
  const double *x0 = (const double *) pos0;
  const double *y0 = (const double *) &pos0[N];
  const double *z0 = (const double *) &pos0[2*N];

  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allcation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
	rupp_sqr[i] = rupp[i]*rupp[i];
  }
	
  for(int i=0;i<N;i++) {
	int j;
	const double xpos = x0[i];
	const double ypos = y0[i];
	const double zpos = z0[i];

	double *x1 = (double *) pos1;
	double *y1 = (double *) &pos1[N];
	double *z1 = (double *) &pos1[2*N];


	for(j=0;j<=(N-block_size);j+=block_size){
	  double dist[block_size];
#ifdef __INTEL_COMPILER			
#pragma simd assert
#endif
	  for(int jj=0;jj<block_size;jj++) {
		const double dx = xpos - x1[jj];
		const double dy = ypos - y1[jj];
		const double dz = zpos - z1[jj];
		dist[jj] = dx*dx + dy*dy + dz*dz;
	  }

	  for(int jj=0;jj<block_size;jj++){
		const double r2 = dist[jj];
		for(int kbin=nrpbin-1;kbin>=1;kbin--){
		  if(r2 >= rupp_sqr[kbin-1]) {
			npairs[kbin]++;
			break;
		  }
		}//searching for kbin loop
	  }
			
	  x1 += block_size;
	  y1 += block_size;
	  z1 += block_size;
			
	}			
	for(int jj=0;j<N;jj++,j++) {
	  const double dx = xpos - x1[jj];
	  const double dy = ypos - y1[jj];
	  const double dz = zpos - z1[jj];
	  const double r2 = dx*dx + dy*dy + dz*dz;
	  for(int kbin=nrpbin-1;kbin>=1;kbin--){
		if(r2 >= rupp_sqr[kbin-1]) {
		  npairs[kbin]++;
		  break;
		}
	  }//searching for kbin loop
	}
  }

  free(rupp);
  return npairs;
}


int64_t * intrinsics_chunked(const double * restrict pos0, const double * restrict pos1, const int N, const char *binfile)
{
  double *x0 = (double *) pos0;
  double *y0 = (double *) &pos0[N];
  double *z0 = (double *) &pos0[2*N];

  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(nrpbin > 0 && "Number of rp bins must be > 0");
  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  AVX_FLOATS m_rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
	rupp_sqr[i] = rupp[i]*rupp[i];
	m_rupp_sqr[i] = AVX_SET_FLOAT(rupp_sqr[i]);
  }

  const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
  const DOUBLE sqr_rpmin = rupp_sqr[0];
	
	
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
	for(j=0;j<=(N-NVEC);j+=NVEC){
	  const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_ALIGNED(x1);
	  const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_ALIGNED(y1);
	  const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_ALIGNED(z1);
	  //set constant := sqr_rpmax
	  const AVX_FLOATS m_sqr_rpmax = AVX_SET_FLOAT(sqr_rpmax);
	  //set constant := sqr_rpmin
	  const AVX_FLOATS m_sqr_rpmin = AVX_SET_FLOAT(sqr_rpmin);
			
	  x1 += NVEC;
	  y1 += NVEC;
	  z1 += NVEC;

	  const AVX_FLOATS m_dx = AVX_SUBTRACT_FLOATS(m_x1, m_xpos);
	  const AVX_FLOATS m_dy = AVX_SUBTRACT_FLOATS(m_y1, m_ypos);
	  const AVX_FLOATS m_dz = AVX_SUBTRACT_FLOATS(m_z1, m_zpos);
			
	  AVX_FLOATS r2 = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dx), AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_dy), AVX_SQUARE_FLOAT(m_dz)));
	  //Create a mask for the NVEC distances that fall within sqr_rpmin and sqr_rpmax (sqr_rpmin <= dist < sqr_rpmax)
	  const AVX_FLOATS m_rpmax_mask = AVX_COMPARE_FLOATS(r2, m_sqr_rpmax, _CMP_LT_OS);
	  const AVX_FLOATS m_rpmin_mask = AVX_COMPARE_FLOATS(r2, m_sqr_rpmin, _CMP_GE_OS);
	  AVX_FLOATS m_mask_left = AVX_BITWISE_AND(m_rpmax_mask,m_rpmin_mask);
	  if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
		continue;
	  }

	  //Update r2 such that all distances that do not satisfy sqr_rpmin <= r2 < sqr_rpmax, get set to sqr_rpmax
	  r2 = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2, m_mask_left);
	  for(int kbin=nrpbin-1;kbin>=1;kbin--) {
		//Create a mask of pairwise separations that are greater than the lower radial limit of this bin (kbin)
		const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_GE_OS);
		//Do a bitwise AND to get the mask for separations that fall into this bin
		const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left);
		//Create the mask for the remainder. This comparison should be exclusive with the comparison used for the m1 variable.
		m_mask_left = AVX_COMPARE_FLOATS(r2,m_rupp_sqr[kbin-1],_CMP_LT_OS);
		/* m_mask_left = AVX_XOR_FLOATS(m1, m_all_ones);//XOR with 0xFFFF... gives the bins that are smaller than m_rupp_sqr[kbin] (and is faster than cmp_p(s/d) in theory) */
		//Check the mask for the separations that fell into this kbin
		const int test2  = AVX_TEST_COMPARISON(m_bin_mask);

		//Do a pop-count to add the number of bits. This is somewhat wasteful, since
		//only 4 bits are set in DOUBLE_PREC mode (8 bits in regular float) but we
		//are adding up all 32 bits in the integer. However, in my massive amount of
		//testing with all sorts of faster implementations of popcount and table lookups,
		//builtin hardware popcnt always outperformed everything else. Thanks to NSA
		//for requiring a hardware popcnt I suppose.
		npairs[kbin] += AVX_BIT_COUNT_INT(test2);

		//Add the kbin variable (as float) into the m_rpbin variable.
		//This would be so much better implemented in AVX2 with support for integers
		//Check if there are any more valid points left. Break out of the kbin histogram loop if none are left
		const int test3 = AVX_TEST_COMPARISON(m_mask_left);
		if(test3 == 0) {
		  break;
		}
	  }
	}			
	for(int jj=0;j<N;jj++,j++) {
	  const double dx = xpos - x1[jj];
	  const double dy = ypos - y1[jj];
	  const double dz = zpos - z1[jj];
	  const double r2 = dx*dx + dy*dy + dz*dz;
	  for(int kbin=nrpbin-1;kbin>=1;kbin--){
		if(r2 >= rupp_sqr[kbin-1]) {
		  npairs[kbin]++;
		  break;
		}
	  }//searching for kbin loop
	}
  }

  free(rupp);
  return npairs;
}


/* int intrinsics_chunked_unroll4(const double * restrict pos0, const double * restrict pos1, const int N, double * restrict d) */
/* { */
/* 	int numcomputed = 0; */
/* 	const int block_size = 4; */
/* #ifdef SQRT_DIST */
/* 	const int unroll_factor=4;//28 sqrt operations are supported -> can be up to 28 */
/* #else */
/* 	const int unroll_factor=4; */
/* #endif	 */
/* 	double *x0 = (double *) pos0; */
/* 	double *y0 = (double *) &pos0[N]; */
/* 	double *z0 = (double *) &pos0[2*N]; */

/* #define GETXVALUE(n) AVX_FLOATS diffx##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&x1[n*block_size]), m_xpos) */
/* #define GETYVALUE(n) AVX_FLOATS diffy##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&y1[n*block_size]), m_ypos) */
/* #define GETZVALUE(n) AVX_FLOATS diffz##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_ALIGNED(&z1[n*block_size]), m_zpos) */
/* #define GETVALUE(n)  GETXVALUE(n);GETYVALUE(n);GETZVALUE(n) */

/* #ifndef SQRT_DIST	 */
/* #define GETDIST(n)  AVX_STREAMING_STORE_FLOATS(&dist[n*block_size],AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n)))) */
/* #else	 */
/* #define GETDIST(n) AVX_STREAMING_STORE_FLOATS(&dist[n*block_size],AVX_SQRT_FLOAT(AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n))))) */
/* #endif */
	
/* 	for(int i=0;i<N;i++) { */

/* 		const double xpos = *x0; */
/* 		const double ypos = *y0; */
/* 		const double zpos = *z0; */
		
/* 		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos); */
/* 		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos); */
/* 		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos); */

/* 		x0++;y0++;z0++; */

/* 		double *x1 = (double *) pos1; */
/* 		double *y1 = (double *) &pos1[N]; */
/* 		double *z1 = (double *) &pos1[2*N]; */

/* 		double *dist = (double *) &d[i*N]; */
		
/* 		int offset = 0; */
/* 		const size_t bytes_offset = block_size*unroll_factor; */
/* 		for(int j=0; j <=(N/block_size-unroll_factor);j+=unroll_factor) { */
/* 			/\* PREFETCH(x1 + 2*block_size*unroll_factor); *\/ */
/* 			/\* PREFETCH(y1 + 2*block_size*unroll_factor); *\/ */
/* 			/\* PREFETCH(z1 + 2*block_size*unroll_factor); *\/ */

			
/* #ifdef SQRT_DIST			 */
/* 			GETVALUE(0);GETVALUE(1);GETVALUE(2);GETVALUE(3);//GETVALUE(4);//GETVALUE(5);//GETVALUE(6);//GETVALUE(7);//GETVALUE(8);GETVALUE(9); */
/* 			/\* GETVALUE(10);GETVALUE(11);GETVALUE(12);GETVALUE(13);GETVALUE(14);GETVALUE(15);GETVALUE(16);GETVALUE(17);GETVALUE(18);GETVALUE(19); *\/ */
/* 			/\* GETVALUE(20);GETVALUE(21);GETVALUE(22);GETVALUE(23);GETVALUE(24);GETVALUE(25);GETVALUE(26);GETVALUE(27); *\/ */
/* 			GETDIST(0);GETDIST(1);GETDIST(2);GETDIST(3);//GETDIST(4);//GETDIST(5);//GETDIST(6);//GETDIST(7);//GETDIST(8);GETDIST(9); */
/*  			/\* GETDIST(10);GETDIST(11);GETDIST(12);GETDIST(13);GETDIST(14);GETDIST(15);GETDIST(16);GETDIST(17);GETDIST(18);GETDIST(19); *\/ */
/* 			/\* GETDIST(20);GETDIST(21);GETDIST(22);GETDIST(23);GETDIST(24);GETDIST(25);GETDIST(26);GETDIST(27); *\/ */
/* #else */
/* 			GETVALUE(0);GETVALUE(1);GETVALUE(2);GETVALUE(3); */
/* 			GETDIST(0);GETDIST(1);GETDIST(2);GETDIST(3); */
/* #endif			 */
/* 			x1 += bytes_offset; */
/* 			y1 += bytes_offset; */
/* 			z1 += bytes_offset; */
/* 			dist += bytes_offset; */
/* 			numcomputed+=bytes_offset; */
/* 			offset += bytes_offset; */
/* 		} */

/* 		dist = (double *) &d[i*N + offset]; */
/* 		for(int jj=0;offset<N;jj++,offset++) { */
/* 			const double dx = xpos - *x1; */
/* 			const double dy = ypos - *y1; */
/* 			const double dz = zpos - *z1; */

/* 			x1++;y1++;z1++; */
			
/* #ifdef SQRT_DIST */
/* 			dist[jj] = sqrt(dx*dx + dy*dy + dz*dz); */
/* #else */
/* 			dist[jj] = dx*dx + dy*dy + dz*dz; */
/* #endif */
/* 			numcomputed++; */
/* 		} */
/* 	} */

/* #undef GETXVALUE */
/* #undef GETYVALUE */
/* #undef GETZVALUE */
/* #undef GETDIST */
	
/* 	return numcomputed; */
/* } */

int main(int argc, char **argv)
{
  const size_t numbytes = NDIM*NELEMENTS;
  double *x    __attribute__((aligned(ALIGNMENT))) = NULL;
  double *y    __attribute__((aligned(ALIGNMENT))) = NULL;
  int test0 = posix_memalign((void **) &x, ALIGNMENT, sizeof(*x)*numbytes);
  int test1 = posix_memalign((void **) &y, ALIGNMENT, sizeof(*y)*numbytes); 

#if (NDIM != 3) 
#error NDIM must be set to 3
#endif
	
  assert(test0 == 0  && test1 == 0 &&  "memory allocation failed");

  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized_chunked","intrinsics_chunked"};
  const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);
  int64_t * (*allfunctions[]) (const double * restrict x, const double * restrict y, const int, const char *)      = {naive,chunked,compiler_vectorized_chunked,intrinsics_chunked};
  const long totflop = (long) NELEMENTS * (long) NELEMENTS * (8);
  const unsigned int seed = 42;
  const int niterations = 5;
  const char *binfile = "./bins";

  double *rupp;
  int Nbins ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&Nbins,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(Nbins > 0 && "Number of rp bins must be > 0");
	
  srand(seed);

  int test_to_run = -1;
  if(argc > 1) {
	test_to_run = atoi(argv[1]);
	if(test_to_run >= ntests) test_to_run = -1;
  }
	
  fill_array(x, NELEMENTS);
  fill_array(y, NELEMENTS);

  int64_t *npairs_reference = calculate_reference_histogram(x, y, NELEMENTS, binfile);

  printf("##################################################\n");
  printf("##  Function                            Time (ms) \n");
  printf("##################################################\n");
  for(int i=0;i<ntests;i++) {
	struct timeval t0,t1;
	double sum_x=0.0, sum_sqr_x=0.0;

	if(i == test_to_run || test_to_run == -1) {
	  //warm-up
	  int64_t *npairs = (allfunctions[i]) (x, y, NELEMENTS, binfile);
			
	  //check-result
	  int numbad = check_result(npairs, npairs_reference, Nbins);
	  if(numbad != 0 ) {
		fprintf(stderr,"ERROR: Number of incorrectly calculated histogram bins = %d out of total %d bins.\n", numbad, Nbins);
		free(npairs);
/* 		goto cleanup; */
	  }

/* 	  free(npairs); */
	  for(int iter=0;iter<niterations;iter++) {
		gettimeofday(&t0,NULL);
		npairs = (allfunctions[i]) (x, y, NELEMENTS, binfile);
		gettimeofday(&t1,NULL);
		free(npairs);
				
		const double this_time_in_ms = 1e3*(t1.tv_sec-t0.tv_sec) +  1e-3*(t1.tv_usec - t0.tv_usec);
		sum_sqr_x += this_time_in_ms*this_time_in_ms;
		sum_x     += this_time_in_ms;
		if(niterations < 10) {
		  printf("     %-35s  %0.2lf \n",allfunction_names[i], this_time_in_ms);
		}
	  }
	  const double mean_time = sum_x/niterations;
	  printf(ANSI_COLOR_RED "# %-35s  %0.2lf +- %0.2lf" ANSI_COLOR_RESET "," ANSI_COLOR_BLUE " %0.2lf GFlops " ANSI_COLOR_RESET "\n",allfunction_names[i], mean_time, sqrt(sum_sqr_x/niterations - mean_time*mean_time), (double) totflop*1e3/mean_time/1e9);
	}
  }

 cleanup:
  {
	free((void *) x);free((void *) y);free(rupp);free(npairs_reference);
	exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
