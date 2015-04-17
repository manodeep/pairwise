#include <stdio.h>

#define __USE_XOPEN2K
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#include "defs.h"
#include "utils.h"
#include "progressbar.h"

#if (NDIM != 3) 
#error NDIM must be set to 3
#endif
	

int64_t * calculate_reference_histogram(const double * restrict pos0, const double * restrict pos1, const int N, const char *file_with_bins)
{
  struct timeval t0, t1;
  gettimeofday(&t0, NULL);
  double *rupp;
  int nrpbin ;
  double rpmin,rpmax;
  setup_bins(file_with_bins,&rpmin,&rpmax,&nrpbin,&rupp);
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

  const double *x0  = (const double *) pos0;
  const double *y0  = (const double *) &pos0[N];
  const double *z0  = (const double *) &pos0[2*N];

  for(int i=0;i<N;i++) {
		const double xpos  = x0[i];
		const double ypos  = y0[i];
		const double zpos  = z0[i];
		const double *x1   = (const double *) pos1;
		const double *y1   = (const double *) &pos1[N];
		const double *z1   = (const double *) &pos1[2*N];

		double *dist = (double *) &d[i*N];
		for(int j=0;j<N;j++) {
			const double dx = xpos - x1[j];
			const double dy = ypos - y1[j];
			const double dz = zpos - z1[j];
			*dist = dx*dx + dy*dy + dz*dz;
			dist++;
		}
  }

	const DOUBLE sqr_rpmin = rupp_sqr[0];
	const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
	
  for(int64_t i=0;i<totnpairs;i++) {
		const DOUBLE r2 = d[i];
		if(r2 < sqr_rpmin || r2 >= sqr_rpmax) continue;
		
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


void naive(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  const double *x0  __attribute__((aligned(ALIGNMENT))) = (const double *) pos0;
  const double *y0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[N];
  const double *z0  __attribute__((aligned(ALIGNMENT))) = (const double *) &pos0[2*N];

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
  }

	const DOUBLE sqr_rpmin = rupp_sqr[0];
	const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
	
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
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
		}
  }

	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}

void chunked(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  const int block_size = 64;
  const double *x0 = (const double *) pos0;
  const double *y0 = (const double *) &pos0[N];
  const double *z0 = (const double *) &pos0[2*N];

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);

  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
  }
	const DOUBLE sqr_rpmin = rupp_sqr[0];
	const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
	
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
					if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
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


	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}


void compiler_vectorized_chunked(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  const int block_size = 4;
  const double *x0 = (const double *) pos0;
  const double *y0 = (const double *) &pos0[N];
  const double *z0 = (const double *) &pos0[2*N];

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allcation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
  }

	const DOUBLE sqr_rpmin = rupp_sqr[0];
	const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
	
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
				if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
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
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
		}
  }

	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);

}

#ifdef __AVX__
void avx_intrinsics_chunked(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  double *x0 = (double *) pos0;
  double *y0 = (double *) &pos0[N];
  double *z0 = (double *) &pos0[2*N];

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
			const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(x1);
			const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(y1);
			const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(z1);
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
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
		}
  }

	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}


void avx_intrinsics_chunked_unroll(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 4;
	const int unroll_factor=2;

  double *x0 = (double *) pos0;
  double *y0 = (double *) &pos0[N];
  double *z0 = (double *) &pos0[2*N];

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
	const AVX_FLOATS m_sqr_rpmax = AVX_SET_FLOAT(sqr_rpmax);
	const AVX_FLOATS m_sqr_rpmin = AVX_SET_FLOAT(sqr_rpmin);
			
#define GETXVALUE(n) const AVX_FLOATS diffx##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&x1[n*block_size]), m_xpos)
#define GETYVALUE(n) const AVX_FLOATS diffy##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&y1[n*block_size]), m_ypos)
#define GETZVALUE(n) const AVX_FLOATS diffz##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&z1[n*block_size]), m_zpos)
#define GETVALUE(n)  GETXVALUE(n);GETYVALUE(n);GETZVALUE(n)

#define GETDIST(n) AVX_FLOATS r2##n = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffx##n),AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(diffy##n), AVX_SQUARE_FLOAT(diffz##n)))
#define UPDATEHIST(n) const AVX_FLOATS m_rpmax_mask##n = AVX_COMPARE_FLOATS(r2##n, m_sqr_rpmax, _CMP_LT_OS); \
	const AVX_FLOATS m_rpmin_mask##n = AVX_COMPARE_FLOATS(r2##n, m_sqr_rpmin, _CMP_GE_OS); \
		AVX_FLOATS m_mask_left##n = AVX_BITWISE_AND(m_rpmax_mask##n,m_rpmin_mask##n); \
			if(AVX_TEST_COMPARISON(m_mask_left##n) != 0) { \
				r2##n = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2##n, m_mask_left##n); \
					for(int kbin##n=nrpbin-1;kbin##n>=1;kbin##n--) {							\
						const AVX_FLOATS m1 = AVX_COMPARE_FLOATS(r2##n,m_rupp_sqr[kbin##n-1],_CMP_GE_OS); \
						const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m1,m_mask_left##n);	\
						m_mask_left##n = AVX_COMPARE_FLOATS(r2##n,m_rupp_sqr[kbin##n-1],_CMP_LT_OS); \
							const int test2  = AVX_TEST_COMPARISON(m_bin_mask);				\
							npairs[kbin##n] += AVX_BIT_COUNT_INT(test2);							\
							const int test3 = AVX_TEST_COMPARISON(m_mask_left##n);		\
							if(test3 == 0) {																					\
								break;																									\
							} \
					} \
			}
#define DOCALC(n) GETVALUE(n);GETDIST(n);UPDATEHIST(n)
	
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

		int offset = 0;
		const size_t bytes_offset = block_size*unroll_factor;
		for(int j=0; j <=(N/block_size-unroll_factor);j+=unroll_factor) {

			//Doing it this way is faster than calling DOCALC(n) n-1 times
			GETVALUE(0);GETVALUE(1);
			GETDIST(0);GETDIST(1);
			UPDATEHIST(0);UPDATEHIST(1);
			
			x1 += bytes_offset;
			y1 += bytes_offset;
			z1 += bytes_offset;
			offset += bytes_offset;
		}

		for(int jj=0;offset<N;jj++,offset++) {
			const double dx = xpos - *x1;
			const double dy = ypos - *y1;
			const double dz = zpos - *z1;
			const double r2 = dx*dx + dy*dy + dz*dz;
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
			x1++;y1++;z1++;
		}
	}

#undef GETXVALUE
#undef GETYVALUE
#undef GETZVALUE
#undef GETDIST
#undef UPDATEHIST
#undef DOCALC
	
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}

#endif //AVX


#if defined (__SSE4_2__)
void sse_intrinsics_chunked(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  double *x0 = (double *) pos0;
  double *y0 = (double *) &pos0[N];
  double *z0 = (double *) &pos0[2*N];

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  SSE_FLOATS m_rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
		m_rupp_sqr[i] = SSE_SET_FLOAT(rupp_sqr[i]);
  }

  const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
  const DOUBLE sqr_rpmin = rupp_sqr[0];
	
	
  for(int i=0;i<N;i++) {
		const double xpos = *x0;
		const double ypos = *y0;
		const double zpos = *z0;
		
		const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
		const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
		const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);

		x0++;y0++;z0++;

		double *x1 = (double *) pos1;
		double *y1 = (double *) &pos1[N];
		double *z1 = (double *) &pos1[2*N];

		int j;		
		for(j=0;j<=(N-SSE_NVEC);j+=SSE_NVEC){
			const SSE_FLOATS m_x1 = SSE_LOAD_FLOATS_UNALIGNED(x1);
			const SSE_FLOATS m_y1 = SSE_LOAD_FLOATS_UNALIGNED(y1);
			const SSE_FLOATS m_z1 = SSE_LOAD_FLOATS_UNALIGNED(z1);
			//set constant := sqr_rpmax
			const SSE_FLOATS m_sqr_rpmax = SSE_SET_FLOAT(sqr_rpmax);
			//set constant := sqr_rpmin
			const SSE_FLOATS m_sqr_rpmin = SSE_SET_FLOAT(sqr_rpmin);
			
			x1 += SSE_NVEC;
			y1 += SSE_NVEC;
			z1 += SSE_NVEC;

			const SSE_FLOATS m_dx = SSE_SUBTRACT_FLOATS(m_x1, m_xpos);
			const SSE_FLOATS m_dy = SSE_SUBTRACT_FLOATS(m_y1, m_ypos);
			const SSE_FLOATS m_dz = SSE_SUBTRACT_FLOATS(m_z1, m_zpos);
			SSE_FLOATS r2 = SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(m_dx), SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(m_dy),SSE_SQUARE_FLOAT(m_dz)));
			SSE_FLOATS m_mask_left;

			{
				m_mask_left = SSE_COMPARE_FLOATS_LT(r2,m_sqr_rpmax);
				if(SSE_TEST_COMPARISON(m_mask_left) == 0) {
					continue;
				}

				SSE_FLOATS m_mask = SSE_BITWISE_AND(m_mask_left, SSE_COMPARE_FLOATS_GE(r2, m_sqr_rpmin));
				if(SSE_TEST_COMPARISON(m_mask) == 0) {
					continue;
				}
				r2 = SSE_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2, m_mask);
				m_mask_left = SSE_COMPARE_FLOATS_LT(r2,m_sqr_rpmax);
			}

			for(int kbin=nrpbin-1;kbin>=1;kbin--) {
				SSE_FLOATS m1 = SSE_COMPARE_FLOATS_GE(r2,m_rupp_sqr[kbin-1]);
				SSE_FLOATS m_bin_mask = SSE_BITWISE_AND(m1,m_mask_left);
				m_mask_left = SSE_COMPARE_FLOATS_LT(r2,m_rupp_sqr[kbin-1]);
				int test2  = SSE_TEST_COMPARISON(m_bin_mask);
				npairs[kbin] += SSE_BIT_COUNT_INT(test2);
				int test3 = SSE_TEST_COMPARISON(m_mask_left);
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
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
		}
  }

	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}


void sse_intrinsics_chunked_unroll(const double * restrict pos0, const double * restrict pos1, const int N, const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 2;
	const int unroll_factor=4;

  double *x0 = (double *) pos0;
  double *y0 = (double *) &pos0[N];
  double *z0 = (double *) &pos0[2*N];

  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  DOUBLE rupp_sqr[nrpbin];
  SSE_FLOATS m_rupp_sqr[nrpbin];
  for(int i=0; i < nrpbin;i++) {
		rupp_sqr[i] = rupp[i]*rupp[i];
		m_rupp_sqr[i] = SSE_SET_FLOAT(rupp_sqr[i]);
  }

  const DOUBLE sqr_rpmax = rupp_sqr[nrpbin-1];
  const DOUBLE sqr_rpmin = rupp_sqr[0];
	const SSE_FLOATS m_sqr_rpmax = SSE_SET_FLOAT(sqr_rpmax);
	const SSE_FLOATS m_sqr_rpmin = SSE_SET_FLOAT(sqr_rpmin);
			
#define GETXVALUE(n) const SSE_FLOATS diffx##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&x1[n*block_size]), m_xpos)
#define GETYVALUE(n) const SSE_FLOATS diffy##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&y1[n*block_size]), m_ypos)
#define GETZVALUE(n) const SSE_FLOATS diffz##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&z1[n*block_size]), m_zpos)
#define GETVALUE(n)  GETXVALUE(n);GETYVALUE(n);GETZVALUE(n)

#define GETDIST(n) SSE_FLOATS r2##n = SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(diffx##n),SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(diffy##n), SSE_SQUARE_FLOAT(diffz##n)))
#define UPDATEHIST(n) const SSE_FLOATS m_rpmax_mask##n = SSE_COMPARE_FLOATS_LT(r2##n, m_sqr_rpmax); \
	const SSE_FLOATS m_rpmin_mask##n = SSE_COMPARE_FLOATS_GE(r2##n, m_sqr_rpmin); \
		SSE_FLOATS m_mask_left##n = SSE_BITWISE_AND(m_rpmax_mask##n,m_rpmin_mask##n); \
			if(SSE_TEST_COMPARISON(m_mask_left##n) != 0) { \
				r2##n = SSE_BLEND_FLOATS_WITH_MASK(m_sqr_rpmax, r2##n, m_mask_left##n); \
					for(int kbin##n=nrpbin-1;kbin##n>=1;kbin##n--) {							\
						const SSE_FLOATS m1 = SSE_COMPARE_FLOATS_GE(r2##n,m_rupp_sqr[kbin##n-1]); \
						const SSE_FLOATS m_bin_mask = SSE_BITWISE_AND(m1,m_mask_left##n);	\
						m_mask_left##n = SSE_COMPARE_FLOATS_LT(r2##n,m_rupp_sqr[kbin##n-1]); \
							const int test2  = SSE_TEST_COMPARISON(m_bin_mask);				\
							npairs[kbin##n] += SSE_BIT_COUNT_INT(test2);							\
							const int test3 = SSE_TEST_COMPARISON(m_mask_left##n);		\
							if(test3 == 0) {																					\
								break;																									\
							} \
					} \
			}
#define DOCALC(n) GETVALUE(n);GETDIST(n);UPDATEHIST(n)

	for(int i=0;i<N;i++) {
		const double xpos = *x0;
		const double ypos = *y0;
		const double zpos = *z0;
		
		const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
		const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
		const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);

		x0++;y0++;z0++;

		double *x1 = (double *) pos1;
		double *y1 = (double *) &pos1[N];
		double *z1 = (double *) &pos1[2*N];

		int offset = 0;
		const size_t bytes_offset = block_size*unroll_factor;
		for(int j=0; j <=(N/block_size-unroll_factor);j+=unroll_factor) {

			DOCALC(0);DOCALC(1);DOCALC(2);DOCALC(3);
			
			x1 += bytes_offset;
			y1 += bytes_offset;
			z1 += bytes_offset;
			offset += bytes_offset;
		}

		for(int jj=0;offset<N;jj++,offset++) {
			const double dx = xpos - *x1;
			const double dy = ypos - *y1;
			const double dz = zpos - *z1;
			const double r2 = dx*dx + dy*dy + dz*dz;
			if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
			for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
			}//searching for kbin loop
			x1++;y1++;z1++;
		}
	}

#undef GETXVALUE
#undef GETYVALUE
#undef GETZVALUE
#undef GETDIST
#undef UPDATEHIST
#undef DOCALC
	
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] = npairs[i];
	}
	free(npairs);
}
#endif //sse 4.2


int main(int argc, char **argv)
{
	int numpart = 0;
	if(argc > 1 ) {
		numpart = atoi(argv[1]);
	} else {
		numpart = NELEMENTS;
	}
	
  const size_t numbytes = NDIM*numpart;
  double *x    __attribute__((aligned(ALIGNMENT))) = NULL;
  double *y    __attribute__((aligned(ALIGNMENT))) = NULL;
  int test0 = posix_memalign((void **) &x, ALIGNMENT, sizeof(*x)*numbytes);
  int test1 = posix_memalign((void **) &y, ALIGNMENT, sizeof(*y)*numbytes); 

  assert(test0 == 0  && test1 == 0 &&  "memory allocation failed");

  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized_chunked"
#ifdef __AVX__																						
																						,"avx_intrinsics_chunked","avx_intrinsics_chunked_unroll"
#endif																						
#ifdef __SSE4_2__
																						,"sse_intrinsics_chunked","sse_intrinsics_chunked_unroll"
#endif
	};
  const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);
  void (*allfunctions[]) (const double * restrict x, const double * restrict y, const int, const int nrpbin, const double *rupp, int64_t * npairs)
		= {naive,chunked,compiler_vectorized_chunked
#ifdef __AVX__			 
			 ,avx_intrinsics_chunked,avx_intrinsics_chunked_unroll
#endif			 
#ifdef __SSE4_2__
			 ,sse_intrinsics_chunked,sse_intrinsics_chunked_unroll
#endif			 
	};
  const double totflop = (double) numpart * (double) numpart * (8);

	double function_best_mean_time[ntests],function_sigma_time[ntests],function_best_time_in_ms[ntests],function_best_mcycles[ntests];
	int function_niterations[ntests];
	for(int i=0;i<ntests;i++) {
		function_best_mean_time[i] = 1e16;
		function_sigma_time[i] = 0.0;
		function_best_time_in_ms[i] = 1e16;
		function_best_mcycles[i] = 1e16;
	}
	
  double *rupp;
  int Nbins ;
  double rpmin,rpmax;
  setup_bins(binfile,&rpmin,&rpmax,&Nbins,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(Nbins > 0 && "Number of rp bins must be > 0");
	
  srand(seed);

	if(clustered_data == 0) {
		fill_array(x, numpart);
		fill_array(y, numpart);
	} else {
		assert(numpart <= max_galaxy_in_source_file && "Clustered data does not contain enough galaxies..please reduce NELEMENTS in defs.h or pass a smaller number on the command-line");
		assert(NDIM == 3 && "Clustered galaxy data contains *exactly* 3 spatial dimensions");
		read_ascii(x, numpart, source_galaxy_file);
		read_ascii(y, numpart, source_galaxy_file);
	}
		
  int64_t *npairs_reference = calculate_reference_histogram(x, y, numpart, binfile);
	int64_t *npairs = my_calloc(sizeof(*npairs),Nbins);
	
	const int64_t totniterations = repeat*ntests*(int64_t) max_niterations;
	int64_t numdone = 0;
	int interrupted=0;

	fprintf(stderr,"# Running benchmarks with N = %05d particles\n",numpart);
	init_my_progressbar(totniterations, &interrupted);

	for(int irep=0;irep<repeat;irep++) {
		for(int i=0;i<ntests;i++) {
			struct timeval t0,t1;
			double sum_x=0.0, sum_sqr_x=0.0;
			
			//warm-up
			(allfunctions[i]) (x, y, numpart, Nbins, rupp, npairs);
			
			//check-result
			int numbad = check_result(npairs, npairs_reference, Nbins);
			if(numbad != 0 ) {
				fprintf(stderr,"ERROR: In function `%s' Number of incorrectly calculated histogram bins = %d out of total %d bins.\n", allfunction_names[i],numbad, Nbins);
				goto cleanup;
			}
			
			double best_time_in_ms=1e16, best_time_in_megacycles=1e16;
			uint64_t start_cycles, end_cycles;
			const int64_t numdone_before_iter_loop = numdone;
			for(int iter=0;iter<max_niterations;iter++) {
				gettimeofday(&t0,NULL);
				start_cycles = rdtsc();
				(allfunctions[i]) (x, y, numpart, Nbins, rupp, npairs);
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
				
				
				if(best_time_in_ms < function_best_time_in_ms[i]) {
					function_best_time_in_ms[i] = best_time_in_ms;
				}
				
				if(best_time_in_megacycles < function_best_mcycles[i]) {
					function_best_mcycles[i] = best_time_in_megacycles;
				}

				if(max_niterations <= 10) {
					printf("     %-35s  %0.2lf max_niterations = %d\n",allfunction_names[i], this_time_in_ms,max_niterations);
					interrupted=1;
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
		printf(ANSI_COLOR_RED "# %-35s  %6.2lf +- %5.2lf " ANSI_COLOR_GREEN " (best -- %6.2lf ms, %6.2lf Mcycles) " ANSI_COLOR_RESET "," ANSI_COLOR_BLUE " >= %5.2lf GFlops [%04d iterations]" ANSI_COLOR_RESET "\n",
					 allfunction_names[i], mean_time, sigma_time, best_time_in_ms, best_time_in_megacycles, gflops, actual_niterations);

	}

cleanup:
  {
		free((void *) x);free((void *) y);free(rupp);free(npairs_reference);free(npairs);
		exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
