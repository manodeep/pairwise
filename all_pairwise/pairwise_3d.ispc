export void pairwise_ispc(const uniform double uniform pos0[], const uniform double uniform pos1[],
																 const uniform int N,
																 uniform double uniform d[])
{
/* 	uniform int Ncomputed = 0; */
	const uniform double * uniform x0 = pos0;
	const uniform double * uniform y0 = pos0 + N;
	const uniform double * uniform z0 = pos0 + 2*N;

	const uniform double * uniform x1 = pos1;
	const uniform double * uniform y1 = pos1 + N;
	const uniform double * uniform z1 = pos1 + 2*N;

	foreach(i=0...N, j=0...N) {
			int index = i*N + j;
			double x1pos = x0[i];
			double y1pos = y0[i];
			double z1pos = z0[i];
			double x2pos = x1[j];
			double y2pos = y1[j];
			double z2pos = z1[j];
			double sqr_dist = (x1pos-x2pos)*(x1pos-x2pos) + (y1pos-y2pos)*(y1pos-y2pos) + (z1pos-z2pos)*(z1pos-z2pos);
#ifdef SQRT_DIST
			d[index] = sqrt(sqr_dist);
#else
			d[index] = sqr_dist;
#endif
		
	}
}


/* inline void mxm(uniform float * uniform A, */
/* 								uniform float * uniform B, */
/* 								uniform float * uniform C, */
/* 								uniform int M, */
/* 								uniform int N, */
/* 								uniform int K, */
/* 								uniform int nmat, */
/* 								int idx) */
/* { */
/* 	for(uniform int i = 0; i < M; i++){ */
/* 		for(uniform int j = 0; j < N; j++){ */
/* 			float sum = 0; */
/* 			for(uniform int k = 0; k < K; k++){ */
/* 				sum += A[i*K*nmat + k*nmat + idx] * B[k*N*nmat + j*nmat + idx]; */
/* 			} */
/* 			C[i*N*nmat + j*nmat + idx] = sum; */
/* 		} */
/* 	} } */
/* export int64 pairwise(uniform double * uniform A, */
/* 											uniform float * uniform B, */
/* 											uniform float * uniform C, */
/* 											uniform int M, */
/* 											uniform int N, */
/* 											uniform int K, */
/* 											uniform int nmat) */
/* { */
/* 	foreach(i = 0 ... nmat){ */
/* 		mxm(A, B, C, M, N, K, nmat, i); */
/* 	} */
/* } */
