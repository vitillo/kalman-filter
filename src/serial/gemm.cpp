#include "common.h"

inline void mxv(scalar_t * A, scalar_t * B, scalar_t * C, int dim, int N){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      C[i] += A[i*dim + j] * B[j];
    }
  }
}

inline void mxm(scalar_t * A, scalar_t * B, scalar_t * C, int dim, int N, int idx){
  scalar_t sum = 0;

  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
	sum += A[i*dim + k] * B[k*dim + j];
      }

      C[i*dim + j] = sum;
    }
  }
}

void mxv_mul(scalar_t A[], scalar_t B[], scalar_t C[], int size, int dim){
  int s = dim * dim;
  for(int i = 0; i < size; i++)
    mxv(A + s, B + dim, C + dim, dim, size);
}

void mxm_mul(scalar_t A[], scalar_t B[], scalar_t C[], int size, int dim){
  int s = dim * dim;
  for(int i = 0; i < size; i++){
    mxm(A + i*s, B + i*s, C + i*dim, dim, size, i);
  }
}


