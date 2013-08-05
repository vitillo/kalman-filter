#include <iostream>
#include <stdlib.h>
#include <immintrin.h>
#include <assert.h>
#include <cstring>
#include <cmath>

#include "ROOT_Ttree_io/event_writer.h"
#include "serial/KalmanFilterSerial.h"
#include "ispc/KalmanFilter.h"

#include "common.h"
#include "clock.h"
#include "test.h"

using namespace std;

void print_soa(const scalar_t *A, int M, int N, int nmat){
  for(int i = 0; i < M*N; i++){
    for(int k = 0; k < nmat; k++){
      cout << A[i*nmat + k] << " ";
    }
    cout << endl;
  }
}

void print_aos(const scalar_t *A, int M, int N, int nmat){
  for(int k = 0; k < nmat; k++){
    for(int i = 0; i < M; i++){
      for(int j = 0; j < N; j++){
	cout << A[k*M*N + i*N + j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}

void aos_to_soa(scalar_t *A, scalar_t *B, int M, int N, int nmat){
  scalar_t *tmp = new scalar_t[M*N*nmat];

  for(int k = 0; k < nmat; k++){
    for(int i = 0; i < M; i++){
      for(int j = 0; j < N; j++){
	tmp[i*N*nmat + j*nmat + k] = A[k*M*N + i*N + j];
      }
    }
  }

  memcpy(B, tmp, sizeof(scalar_t)*M*N*nmat);
  delete tmp;
}

void soa_to_aos(scalar_t *A, scalar_t *B, int M, int N, int nmat){
  scalar_t *tmp = new scalar_t[M*N*nmat];

  for(int k = 0; k < nmat; k++){
    for(int i = 0; i < M; i++){
      for(int j = 0; j < N; j++){
	tmp[k*M*N + i*N + j] = A[i*N*nmat + j*nmat + k];
      }
    }
  }

  memcpy(B, tmp, sizeof(scalar_t)*M*N*nmat);
  delete tmp;
}

void initialize_aos(scalar_t *A, int M, int N, int nmat){
  for(int k = 0; k < nmat; k++){
    for(int i = 0; i < M; i++){
      for(int j = 0; j < N; j++){
	A[k*N*M + i*N + j] = i*N + j;
      }
    }
  }
}

void initialize_scalar(scalar_t *A, int nmat){
  for(int i = 0; i < nmat; i++)
    A[i] = 3;
}

void compare(scalar_t *A, scalar_t *B, int M, int N, int nmat){
  for(int k = 0; k < nmat; k++){
    for(int i = 0; i < M; i++){
      for(int j = 0; j < N; j++){
	assert(abs(A[k*M*N + i*N + j] - B[i*N + j]) < 0.00001);
      }
    }
  }
}

void testLinearAlgebra(){
  const int nmat = 16, M = 3, N = 2;

  scalar_t *A = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat*M*N, ALIGNMENT);
  scalar_t *AA = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat*M*N, ALIGNMENT);
  scalar_t *B = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat*M*N, ALIGNMENT);
  scalar_t *BB = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat*M*N, ALIGNMENT);
  scalar_t *C = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat*M*M, ALIGNMENT);
  scalar_t *s = (scalar_t *) _mm_malloc(sizeof(scalar_t)*nmat, ALIGNMENT);

  initialize_aos(A, M, N, nmat);
  initialize_aos(B, N, M, nmat);
  initialize_scalar(s, nmat);

  aos_to_soa(A, AA, M, N, nmat);
  aos_to_soa(B, BB, N, M, nmat);
  ispc::gemm_test(0, 0, AA, BB, C, M, M, N, nmat);
  soa_to_aos(C, C, M, M, nmat);
  scalar_t r[] = {3, 4, 5, 9, 14, 19, 15, 24, 33};
  compare(C, r, M, M, nmat);

  ispc::gemm_test(0, 1, AA, BB, C, M, M, N, nmat);
  soa_to_aos(C, C, M, M, nmat);
  scalar_t r1[] = {1, 3, 5, 3, 13, 23, 5, 23, 41};
  compare(C, r1, M, M, nmat);

  ispc::gemm_test(1, 0, AA, BB, C, N, N, M, nmat);
  soa_to_aos(C, C, N, N, nmat);
  scalar_t r11[] = {20, 26, 26, 35};
  compare(C, r11, N, N, nmat);

  initialize_aos(AA, 2, 2, nmat);
  aos_to_soa(AA, C, 2, 2, nmat);
  ispc::inverse2x2_test(C, C, nmat);
  soa_to_aos(C, C, 2, 2, nmat);
  scalar_t r12[] = {-1.5, 0.5, 1, 0};
  compare(C, r12, 2, 2, nmat);

  aos_to_soa(A, C, N, M, nmat);
  ispc::add_test(C, C, C, M, N, nmat);
  soa_to_aos(C, C, M, N, nmat);
  scalar_t r2[] = {0, 2, 4, 6, 8, 10};
  compare(C, r2, M, N, nmat);

  aos_to_soa(A, C, N, M, nmat);
  ispc::subtract_test(C, C, C, M, N, nmat);
  soa_to_aos(C, C, M, N, nmat);
  scalar_t r3[] = {0, 0, 0, 0, 0, 0};
  compare(C, r3, M, N, nmat);

  aos_to_soa(A, C, N, M, nmat);
  ispc::multiply_test(C, s, C, M, N, nmat);
  soa_to_aos(C, C, M, N, nmat);
  scalar_t r4[] = {0, 3, 6, 9, 12, 15};
  compare(C, r4, M, N, nmat);

  aos_to_soa(A, C, N, M, nmat);
  ispc::divide_test(C, s, C, M, N, nmat);
  soa_to_aos(C, C, M, N, nmat);
  scalar_t r5[] = {0, 0.333333, 0.666667, 1, 1.33333, 1.66667};
  compare(C, r5, M, N, nmat);

  ispc::identity_test(C, M, nmat);
  soa_to_aos(C, C, M, M, nmat);
  scalar_t r6[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  compare(C, r6, M, M, nmat);
}

void testKalmanFilterNTracksMHitsMixed(int ntracks, int hits){
  const double eps = 0.1;

  KalmanFilterSerial kalmanFilter1D, kalmanFilter2D;
  Track_t track1D;
  initializeTrack(track1D, hits, false);
  kalmanFilter1D.setTrackData(track1D);
  kalmanFilter1D.startFilter();

  Track_t track2D;
  initializeTrack(track2D, hits, true);
  kalmanFilter2D.setTrackData(track2D);
  kalmanFilter2D.startFilter();

  ispc::KalmanFilter filter;
  ispc::KalmanFilterParameter param;
  ispc::KalmanFilter_initialize(&filter, ntracks);
  KalmanFilterParallelNTracksMHitsMixed(&param, ntracks, hits);
  ispc::startFilter(&filter, &param);

  for(int t = 0; t < ntracks; t++){
    KalmanFilterSerial &serial = (t % 2 == 0) ? kalmanFilter1D : kalmanFilter2D;

    for(int i = 0; i < 5; i++){                                      
      assert(abs(GSL_VECTOR_GET(serial.mP_k1, i) - filter.P_k1[i*ntracks + t]) < eps);
      for(int j = 0; j < 5; j++){
	assert(abs(GSL_MATRIX_GET(serial.mC_k1, i, j) - filter.C_k1[i*5*ntracks + j*ntracks + t]) < eps);
      }
    }
  }

  ispc::KalmanFilter_deallocate(&filter);
  ispc::KalmanFilterParameter_deallocate(&param);
}

int main(){
  testLinearAlgebra();
  testKalmanFilterNTracksMHitsMixed(1, 1);
  testKalmanFilterNTracksMHitsMixed(1, 20);
  testKalmanFilterNTracksMHitsMixed(1000, 1);
  testKalmanFilterNTracksMHitsMixed(1000, 20);
  return 0;
}
