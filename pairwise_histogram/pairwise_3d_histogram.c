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
  int num_zero = 0;
  
  for(int i=1;i<Nbins;i++) {
		if(npairs[i] != npairs_reference[i]) {
			if(numbadprinted < 20) {
				fprintf(stderr,"i = %d npairs = %"PRId64" npairs_reference = %"PRId64" \n",i,npairs[i],npairs_reference[i]);
				numbadprinted++;
			}
			bad++;
		}

    if(npairs[i] == 0)
      num_zero++;
  }
	if(num_zero > 0) {
    fprintf(stderr,"Strange: found %d zeros in the pair counts\n", num_zero);
  }
  return bad;
}


static inline void naive(const double * restrict x0, const double * restrict y0, const double * restrict z0,
												 const double * restrict x1, const double * restrict y1, const double * restrict z1,
												 const int N0, const int N1,
												 const int nrpbin, const double *rupp, int64_t *results_npairs)
{
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
	
  for(int i=0;i<N0;i++) {
		const double xpos  __attribute__((aligned(ALIGNMENT))) = x0[i];
		const double ypos  __attribute__((aligned(ALIGNMENT))) = y0[i];
		const double zpos  __attribute__((aligned(ALIGNMENT))) = z0[i];

		for(int j=0;j<N1;j++) {
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
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}

static inline void chunked(const double * restrict x0, const double * restrict y0, const double * restrict z0,
													 const double * restrict x1, const double * restrict y1, const double * restrict z1,
													 const int N0, const int N1,
													 const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  const int block_size = 64;

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
	
  for(int i=0;i<N0;i+=block_size) {
		const int block_size1 = (N0-i) > block_size ? block_size:(N0-i);
		for(int j=0;j<N1;j+=block_size) {
			const int block_size2 = (N1-j) > block_size ? block_size:(N1-j);
			for(int ii=0;ii<block_size1;ii++) {
				double dist[block_size2];
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
					dist[jj] = dx*dx + dy*dy + dz*dz;
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
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}


static inline void compiler_vectorized(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																			 const double * restrict x1, const double * restrict y1, const double * restrict z1,
																			 const int N0, const int N1,
																			 const int nrpbin, const double *rupp, int64_t *results_npairs)
{
  const int block_size = 4;
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
	
  for(int i=0;i<N0;i++) {
		int j;
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;
		
		for(j=0;j<=(N1-block_size);j+=block_size){
			double dist[block_size];
#ifdef __INTEL_COMPILER			
#pragma simd assert
#endif
			for(int jj=0;jj<block_size;jj++) {
				const double dx = xpos - localx1[jj];
				const double dy = ypos - localy1[jj];
				const double dz = zpos - localz1[jj];
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
			
			localx1 += block_size;
			localy1 += block_size;
			localz1 += block_size;
		}			
		for(int jj=0;j<N1;jj++,j++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
			
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
		results_npairs[i] += npairs[i];
	}
	free(npairs);

}

#ifdef __AVX__
static inline void avx_intrinsics(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																	const double * restrict x1, const double * restrict y1, const double * restrict z1,
																	const int N0, const int N1,
																	const int nrpbin, const double *rupp, int64_t *results_npairs)
{
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
	
	double *localx0 = (double *) x0;
	double *localy0 = (double *) y0;
	double *localz0 = (double *) z0;
	
  for(int i=0;i<N0;i++) {
		/* const double xpos = x0[i]; */
		/* const double ypos = y0[i]; */
		/* const double zpos = z0[i]; */
		
		const AVX_FLOATS m_xpos = AVX_SET_FLOAT(*localx0);
		const AVX_FLOATS m_ypos = AVX_SET_FLOAT(*localy0);
		const AVX_FLOATS m_zpos = AVX_SET_FLOAT(*localz0);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		int j;		
		for(j=0;j<=(N1-NVEC);j+=NVEC){
			const AVX_FLOATS m_x1 = AVX_LOAD_FLOATS_UNALIGNED(localx1);
			const AVX_FLOATS m_y1 = AVX_LOAD_FLOATS_UNALIGNED(localy1);
			const AVX_FLOATS m_z1 = AVX_LOAD_FLOATS_UNALIGNED(localz1);
			//set constant := sqr_rpmax
			const AVX_FLOATS m_sqr_rpmax = AVX_SET_FLOAT(sqr_rpmax);
			//set constant := sqr_rpmin
			const AVX_FLOATS m_sqr_rpmin = AVX_SET_FLOAT(sqr_rpmin);
			
			localx1 += NVEC;
			localy1 += NVEC;
			localz1 += NVEC;

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

		for(int jj=0;j<N1;jj++,j++) {
		  const double dx = *localx0 - localx1[jj];
		  const double dy = *localy0 - localy1[jj];
		  const double dz = *localz0 - localz1[jj];
		  const double r2 = dx*dx + dy*dy + dz*dz;
		  if(r2 >= sqr_rpmax || r2 < sqr_rpmin) continue;
		  for(int kbin=nrpbin-1;kbin>=1;kbin--){
				if(r2 >= rupp_sqr[kbin-1]) {
					npairs[kbin]++;
					break;
				}
		  }//searching for kbin loop
		}
		localx0++;localy0++;localz0++;
	}

	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}


static inline void avx_intrinsics_unroll(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																				 const double * restrict x1, const double * restrict y1, const double * restrict z1,
																				 const int N0, const int N1,
																				 const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 4;
	const int unroll_factor=2;
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
			
#define GETXVALUE(n) const AVX_FLOATS diffx##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localx1[n*block_size]), m_xpos)
#define GETYVALUE(n) const AVX_FLOATS diffy##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localy1[n*block_size]), m_ypos)
#define GETZVALUE(n) const AVX_FLOATS diffz##n = AVX_SUBTRACT_FLOATS(AVX_LOAD_FLOATS_UNALIGNED(&localz1[n*block_size]), m_zpos)
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

		int offset = 0;
		const size_t bytes_offset = block_size*unroll_factor;
		for(int j=0; j <=(N1/block_size-unroll_factor);j+=unroll_factor) {

			//Doing it this way is faster than calling DOCALC(n) n-1 times
			GETVALUE(0);GETVALUE(1);
			GETDIST(0);GETDIST(1);
			UPDATEHIST(0);UPDATEHIST(1);
			
			localx1 += bytes_offset;
			localy1 += bytes_offset;
			localz1 += bytes_offset;
			offset += bytes_offset;
		}

		for(int jj=0;offset<N1;jj++,offset++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
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

#undef GETXVALUE
#undef GETYVALUE
#undef GETZVALUE
#undef GETDIST
#undef UPDATEHIST
#undef DOCALC
	
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}

static inline void avx_intrinsics_chunked(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																					const double * restrict x1, const double * restrict y1, const double * restrict z1,
																					const int N0, const int N1,
																					const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 1024;
  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  for(int i=0;i<N0;i+=block_size) {
		const int block_size1 = (N0-i) > block_size ? block_size:(N0-i);
		for(int j=0;j<N1;j+=block_size) {
			const int block_size2 = (N1-j) > block_size ? block_size:(N1-j);
			avx_intrinsics_unroll(&x0[i],&y0[i],&z0[i],
														&x1[j],&y1[j],&z1[j],
														block_size1,block_size2,
														nrpbin, rupp,npairs);
		}
	}
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}


#endif //AVX


#if defined (__SSE4_2__)
void sse_intrinsics(const double * restrict x0, const double * restrict y0, const double * restrict z0,
										const double * restrict x1, const double * restrict y1, const double * restrict z1,
										const int N0, const int N1,
										const int nrpbin, const double *rupp, int64_t *results_npairs)
{
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
	
	
  for(int i=0;i<N0;i++) {
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];
		
		const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
		const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
		const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		int j;		
		for(j=0;j<=(N1-SSE_NVEC);j+=SSE_NVEC){
			const SSE_FLOATS m_x1 = SSE_LOAD_FLOATS_UNALIGNED(localx1);
			const SSE_FLOATS m_y1 = SSE_LOAD_FLOATS_UNALIGNED(localy1);
			const SSE_FLOATS m_z1 = SSE_LOAD_FLOATS_UNALIGNED(localz1);
			//set constant := sqr_rpmax
			const SSE_FLOATS m_sqr_rpmax = SSE_SET_FLOAT(sqr_rpmax);
			//set constant := sqr_rpmin
			const SSE_FLOATS m_sqr_rpmin = SSE_SET_FLOAT(sqr_rpmin);
			
			localx1 += SSE_NVEC;
			localy1 += SSE_NVEC;
			localz1 += SSE_NVEC;

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

		for(int jj=0;j<N1;jj++,j++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
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
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}


void sse_intrinsics_unroll(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																	 const double * restrict x1, const double * restrict y1, const double * restrict z1,
																	 const int N0, const int N1,
																	 const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 2;
	const int unroll_factor=4;
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
			
#define GETXVALUE(n) const SSE_FLOATS diffx##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&localx1[n*block_size]), m_xpos)
#define GETYVALUE(n) const SSE_FLOATS diffy##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&localy1[n*block_size]), m_ypos)
#define GETZVALUE(n) const SSE_FLOATS diffz##n = SSE_SUBTRACT_FLOATS(SSE_LOAD_FLOATS_UNALIGNED(&localz1[n*block_size]), m_zpos)
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

	for(int i=0;i<N0;i++) {
		const double xpos = x0[i];
		const double ypos = y0[i];
		const double zpos = z0[i];
		
		const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
		const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
		const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);

		double *localx1 = (double *) x1;
		double *localy1 = (double *) y1;
		double *localz1 = (double *) z1;

		int offset = 0;
		const size_t bytes_offset = block_size*unroll_factor;
		for(int j=0; j <=(N1/block_size-unroll_factor);j+=unroll_factor) {

			DOCALC(0);DOCALC(1);DOCALC(2);DOCALC(3);
			
			localx1 += bytes_offset;
			localy1 += bytes_offset;
			localz1 += bytes_offset;
			offset += bytes_offset;
		}

		for(int jj=0;offset<N1;jj++,offset++) {
			const double dx = xpos - localx1[jj];
			const double dy = ypos - localy1[jj];
			const double dz = zpos - localz1[jj];
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

#undef GETXVALUE
#undef GETYVALUE
#undef GETZVALUE
#undef GETDIST
#undef UPDATEHIST
#undef DOCALC
	
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] += npairs[i];
	}
	free(npairs);
}

static inline void sse_intrinsics_chunked(const double * restrict x0, const double * restrict y0, const double * restrict z0,
																					const double * restrict x1, const double * restrict y1, const double * restrict z1,
																					const int N0, const int N1,
																					const int nrpbin, const double *rupp, int64_t *results_npairs)
{
	const int block_size = 1024;
  int64_t *npairs = NULL;
  size_t numbytes = sizeof(*npairs)*nrpbin;
  const int test0 = posix_memalign((void **) &npairs, ALIGNMENT, numbytes);
  assert(test0 == 0 && "memory allocation failed");
  memset(npairs, 0, numbytes);
	
  for(int i=0;i<N0;i+=block_size) {
		const int block_size1 = (N0-i) > block_size ? block_size:(N0-i);
		for(int j=0;j<N1;j+=block_size) {
			const int block_size2 = (N1-j) > block_size ? block_size:(N1-j);
			sse_intrinsics(&x0[i],&y0[i],&z0[i],
										 &x1[j],&y1[j],&z1[j],
										 block_size1,block_size2,
										 nrpbin, rupp,npairs);
		}
	}
	for(int i=0;i<nrpbin;i++) {
		results_npairs[i] += npairs[i];
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
  double *pos0    __attribute__((aligned(ALIGNMENT))) = NULL;
  double *pos1    __attribute__((aligned(ALIGNMENT))) = NULL;
  int test0 = posix_memalign((void **) &pos0, ALIGNMENT, sizeof(*pos0)*numbytes);
  int test1 = posix_memalign((void **) &pos1, ALIGNMENT, sizeof(*pos1)*numbytes); 

  assert(test0 == 0  && test1 == 0 &&  "memory allocation failed");

  const char allfunction_names[][MAXLEN] = {"naive","chunked","compiler_vectorized"
#ifdef __AVX__																						
																						,"avx_intrinsics","avx_intrinsics_unroll","avx_intrinsics_chunked"
#endif																						
#ifdef __SSE4_2__
																						,"sse_intrinsics","sse_intrinsics_unroll","sse_intrinsics_chunked"
#endif
	};

	const int ntests = sizeof(allfunction_names)/(sizeof(char)*MAXLEN);



	//Seriously this is the declaration for the function pointers..
	void (*allfunctions[]) ( const double * restrict x0, const double * restrict y0, const double * restrict z0,
													 const double * restrict x1, const double * restrict y1, const double * restrict z1,
													 const int N0, const int N1,
													 const int nrpbin, const double *rupp, int64_t *results_npairs)
		= {naive,chunked,compiler_vectorized
#ifdef __AVX__			 
			 ,avx_intrinsics,avx_intrinsics_unroll,avx_intrinsics_chunked
#endif			 
#ifdef __SSE4_2__
			 ,sse_intrinsics,sse_intrinsics_unroll,sse_intrinsics_chunked
#endif			 
	};
	//end of function pointers declaration.
	
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
  char file_with_bins[MAXLEN];
  my_snprintf(file_with_bins, MAXLEN,"%s",binfile);
  FILE *fp = fopen(file_with_bins,"r");
  if(fp == NULL) {
	my_snprintf(file_with_bins, MAXLEN,"../%s",binfile);
	fp = fopen(file_with_bins,"r");
	if(fp == NULL) {
	  fprintf(stderr,"ERROR: Could not find file containing radial bins either as `%s' or as `%s'..exiting\n",
			  binfile, file_with_bins);
		exit(EXIT_FAILURE);
	}
  } 

  if(fp != NULL) {
	fclose(fp);
  }


  setup_bins(file_with_bins,&rpmin,&rpmax,&Nbins,&rupp);
  assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
  assert(Nbins > 0 && "Number of rp bins must be > 0");
	
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
		
  int64_t *npairs_reference = calculate_reference_histogram(pos0,pos1, numpart, file_with_bins);
	int64_t *npairs = my_calloc(sizeof(*npairs),Nbins);
	
	const int64_t totniterations = repeat*ntests*(int64_t) max_niterations;
	int64_t numdone = 0;
	int interrupted=0;

	fprintf(stderr,"# Running benchmarks with N = %06d particles\n",numpart);
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
			
			memset(npairs,0,sizeof(*npairs)*Nbins);
			//warm-up
			(allfunctions[i]) (x0,y0,z0,x1,y1,z1, numpart,numpart, Nbins, rupp, npairs);
			
			//check-result
			int numbad = check_result(npairs, npairs_reference, Nbins);
			if(numbad != 0 ) {
				fprintf(stderr,"ERROR: In function `%s' Number of incorrectly calculated histogram bins = %d out of total %d bins.\n", allfunction_names[i],numbad, Nbins);
				interrupted=1;
				goto cleanup;
			}

      fflush(stdout);
      fprintf(stdout,"%s\n",allfunction_names[i]);
      for(int j=1;j<Nbins;j++) {
        fprintf(stdout,"%12.6lf   %12.6lf  %12"PRId64"\n", rupp[j-1], rupp[j], npairs[j]);
      }
      fprintf(stdout,"--------------------------------------------\n");
      fflush(stdout);
      
			double best_time_in_ms=1e16, best_time_in_megacycles=1e16;
			uint64_t start_cycles, end_cycles;
			const int64_t numdone_before_iter_loop = numdone;
			for(int iter=0;iter<max_niterations;iter++) {
				gettimeofday(&t0,NULL);
				start_cycles = rdtsc();
				(allfunctions[i]) (x0,y0,z0,x1,y1,z1,numpart, numpart, Nbins, rupp, npairs);
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
		free((void *) pos0);free((void *) pos1);free(rupp);free(npairs_reference);free(npairs);
		exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
