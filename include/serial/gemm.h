#ifndef GEMM_H
#define GEMM_H

void mxv(scalar_t * A, scalar_t * B, scalar_t * C, int dim, int N);
void mxm(scalar_t * A, scalar_t * B, scalar_t * C, int dim, int N, int idx);
void mxv_mul(scalar_t A[], scalar_t B[], scalar_t C[], int size, int dim);
void mxm_mul(scalar_t A[], scalar_t B[], scalar_t C[], int size, int dim);

#endif
