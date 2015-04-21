//
// all_pairwise/pairwise_3d_ispc.h
// (Header automatically generated by the ispc compiler.)
// DO NOT EDIT THIS FILE.
//

#ifndef ISPC_ALL_PAIRWISE_PAIRWISE_3D_ISPC_H
#define ISPC_ALL_PAIRWISE_PAIRWISE_3D_ISPC_H

#include <stdint.h>



#ifdef __cplusplus
namespace ispc { /* namespace */
#endif // __cplusplus

///////////////////////////////////////////////////////////////////////////
// Functions exported from ispc code
///////////////////////////////////////////////////////////////////////////
#if defined(__cplusplus) && !defined(__ISPC_NO_EXTERN_C)
extern "C" {
#endif // __cplusplus
    extern void pairwise_ispc(const double * x0, const double * y0, const double * z0, const double * x1, const double * y1, const double * z1, const int32_t N0, const int32_t N1, double * d);
#if defined(__cplusplus) && !defined(__ISPC_NO_EXTERN_C)
} /* end extern C */
#endif // __cplusplus


#ifdef __cplusplus
} /* namespace */
#endif // __cplusplus

#endif // ISPC_ALL_PAIRWISE_PAIRWISE_3D_ISPC_H
