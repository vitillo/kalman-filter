#include <iostream>
#include <stdlib.h>
#include <immintrin.h>
#include <assert.h>

#include "common.h"
#include "clock.h"
#include "test.h"
#include "ispc/KalmanFilter.h"
#include "serial/gemm.h"

using namespace std;

#define DIM 5
#define NMAT 1000
#define ITER 10000

void mm(){
  scalar_t *A = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM*DIM, ALIGNMENT);
  scalar_t *B = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM*DIM, ALIGNMENT);
  scalar_t *C = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM*DIM, ALIGNMENT);
  double ts, tp;

  cout << "Matrix" << DIM << "/Matrix" << DIM << " Multiplication" << endl;

  start_clock();
  for(int i = 0; i < ITER; i++)
    ispc::gemm_test(0, 0, A, B, C, DIM, DIM, DIM, NMAT);
  tp = stop_clock();

  start_clock();
  for(int i = 0; i < ITER; i++)
    mxm_mul(A, B, C, NMAT, DIM);
  ts = stop_clock();
  cout << "speedup: " << ts/tp << endl;
}

void mv() {
  scalar_t *A = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM*DIM, ALIGNMENT);
  scalar_t *B = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM, ALIGNMENT);
  scalar_t *C = (scalar_t *) _mm_malloc(sizeof(scalar_t)*NMAT*DIM, ALIGNMENT);
  double ts, tp;

  cout << "Matrix" << DIM << "/Vector" << DIM << " Multiplication" << endl;
  
  for(int i = 0; i < ITER; i++)
    ispc::gemm_test(0, 0, A, B, C, DIM, 1, DIM, NMAT);

  start_clock();
  for(int i = 0; i < ITER; i++)
    ispc::gemm_test(0, 0, A, B, C, DIM, 1, DIM, NMAT);
  tp = stop_clock();

  start_clock();
  for(int i = 0; i < ITER; i++)
    mxv_mul(A, B, C, NMAT, DIM);
  ts = stop_clock();
  cout << "speedup: " << ts/tp << endl;
}

void kalman(){
  cout << "Kalman Filter" << endl;

  double t1, tn, tn1;
  int hits = 15;
  int ntracks = 1000;

  Track_t track1D, track2D;
  initializeTrack(track1D, hits, false);
  initializeTrack(track2D, hits, true);
  KalmanFilterSerial kalmanFilter;
  kalmanFilter.setTrackData(track1D);

  start_clock();
  for(int j = 0; j < ntracks; j++){
    kalmanFilter.setTrackData((j & 1) ? track1D : track2D);
    kalmanFilter.startFilter();
  }
  t1 = stop_clock();

  ispc::KalmanFilter filter;
  ispc::KalmanFilterParameter param;
  ispc::KalmanFilter_initialize(&filter, ntracks);
  KalmanFilterParallelNTracksMHitsMixed(&param, ntracks, hits);

  start_clock();
  ispc::startFilter(&filter, &param);
  tn = stop_clock();

  ispc::KalmanFilter filter_serial;
  ispc::KalmanFilterParameter param_serial;
  ispc::KalmanFilter_initialize(&filter_serial, 1);
  KalmanFilterParallelNTracksMHitsMixed(&param_serial, 1, hits);

  start_clock();
  for(int j = 0; j < ntracks; j++)
    ispc::startFilter(&filter_serial, &param_serial);
  tn1 = stop_clock();

  ispc::KalmanFilter_deallocate(&filter);
  ispc::KalmanFilterParameter_deallocate(&param);

  cout << "absolute speedup: " << t1/tn << endl;
  cout << "relative speedup: " << tn1/tn << endl;
}

int main(){
  mv();
  mm();
  kalman();
}
