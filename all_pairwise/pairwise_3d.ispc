
export void pairwise_ispc(const uniform double uniform x0[], const uniform double uniform y0[], const uniform double uniform z0[],
													const uniform double uniform x1[], const uniform double uniform y1[], const uniform double uniform z1[], 
													const uniform int N0, const uniform int N1, 
													uniform double uniform d[])
{
	foreach(i=0...N0, j=0...N1) {
		int64 index = i*(int64) N1 + j;
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
