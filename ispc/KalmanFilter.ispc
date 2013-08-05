#include "common.h"

#define M5(matrix, i, j, lane_id) filter->matrix[i*5*filter->ntracks + j*filter->ntracks + lane_id]
#define M2(matrix, i, j, lane_id) filter->matrix[i*2*filter->ntracks + j*filter->ntracks + lane_id]
#define M1(matrix, i, lane_id) filter->matrix[i*filter->ntracks + lane_id]
#define V(matrix, i, lane_id) filter->matrix[i*filter->ntracks + lane_id]

struct TrackHits {
  scalar_t *normal;
  scalar_t *ref;
  scalar_t *err_locX;
  scalar_t *err_locY;
  scalar_t *cov_locXY;
  scalar_t *Jacobi;
  int *is2Dim;
};

struct KalmanFilterParameter {
  int max_hit_storage;
  int max_track_storage;

  int max_hit_count;
  int ntracks;
  int *hit_count;

  TrackHits *hits;
  /*trackHitData *fitsForward;
  trackHitData *fitsBackward;
  trackHitData *fitsResult;
  scalar_t *C_Forward;
  scalar_t *C_Inverse; // also stores result in the end (after smoothing)
  scalar_t *C_Backward;
  TrackInfo_t* fittedPerigee;
  scalar_t *C_fittedPerigee;*/
};

struct KalmanFilter{
  scalar_t *C_k1;      // 5x5 Initial covariance estimate used for the prediction.
  scalar_t *C_k;       // 5x5 The predicted covariance.

  scalar_t *Dim1_H;    // 1x5 measurment matrix to map from track parameter space to measurement space
  scalar_t *Dim2_H;    // 2x5 measurement matrix

  scalar_t *Dim1_K;    // 5x1 Kalman gain matrix.
  scalar_t *Dim2_K;    // 5x2 Kalman gain matrix.
  scalar_t *Dim5_K;

  scalar_t *Dim1_V;    // 1x1 variance matrix.
  scalar_t *Dim2_V;    // 2x2 variance matrix.

  scalar_t *Dim1_m;    // 1D measurement vector
  scalar_t *Dim2_m;    // 2D measurement vector.

  scalar_t *P_k1;      // 5x1 Initial state estimate used for the prediction.
  scalar_t *P_k;       // 5x1 The predicted state.

  scalar_t *M55_tmp;
  scalar_t *M15_tmp;
  scalar_t *M11_tmp;
  scalar_t *M55_tmp1;
  scalar_t *M25_tmp;
  scalar_t *M22_tmp;

  scalar_t *v2_tmp;
  scalar_t *v5_tmp;

  scalar_t *Inverse;

  int ntracks;
  int max_track_storage;
};

inline void mxm(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int K, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < N; j++){
      scalar_t sum = 0;

      for(uniform int k = 0; k < K; k++){
	sum += A[i*K*num + k*num + idx] * B[k*N*num + j*num + idx];
      }

      C[i*N*num + j*num + idx] = sum;
    }
  }
}

inline void mtxm(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int K, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < N; j++){
      scalar_t sum = 0;

      for(uniform int k = 0; k < K; k++){
	sum += A[k*M*num + i*num + idx] * B[k*N*num + j*num + idx];
      }

      C[i*N*num + j*num + idx] = sum;
    }
  }
}

inline void mxmt(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int K, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < N; j++){
      scalar_t sum = 0;

      for(uniform int k = 0; k < K; k++){
	sum += A[i*K*num + k*num + idx] * B[j*K*num + k*num + idx];
      }

      C[i*N*num + j*num + idx] = sum;
    }
  }
}

inline void gemm(uniform int transpose_A, uniform int transpose_B, uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int K, uniform int num, int idx){
  if(!transpose_A && !transpose_B)
    mxm(A, B, C, M, N, K, num, idx);
  else if(transpose_B && !transpose_A)
    mxmt(A, B, C, M, N, K, num, idx);
  else if(transpose_A && !transpose_B)
    mtxm(A, B, C, M, N, K, num, idx);
  else
    assert(0); //not implemented
}

inline void add(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    uniform int row = i*N*num;

    for(uniform int j = 0; j < N; j++){
      uniform int col = j*num;

      C[row + col + idx] = A[row + col + idx] + B[row + col + idx];
    }
  }
}

inline void subtract(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    uniform int row = i*N*num;

    for(uniform int j = 0; j < N; j++){
      uniform int col = j*num;

      C[row + col + idx] = A[row + col + idx] - B[row + col + idx];
    }
  }
}

inline void multiply(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < N; j++){
      C[i*N*num + j*num + idx] = A[i*N*num + j*num + idx] * B[idx];
    }
  }
}

inline void divide(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < N; j++){
      C[i*N*num + j*num + idx] = A[i*N*num + j*num + idx] / B[idx];
    }
  }
}

inline void identity(uniform scalar_t * uniform A, uniform int M, uniform int num, int idx){
  for(uniform int i = 0; i < M; i++){
    for(uniform int j = 0; j < M; j++){
      A[i*M*num + j*num + idx] = (i == j) ? 1. : 0.;
    }
  }
}

inline void inverse2x2(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform int num, int idx){
  scalar_t a = A[idx];
  scalar_t b = A[num + idx];
  scalar_t c = A[2*num + idx];
  scalar_t d = A[3*num + idx];

  scalar_t det_inv = 1./(a*d - b*c);
  B[idx] = d*det_inv;
  B[num + idx] = -b*det_inv;
  B[2*num + idx] = -c*det_inv;
  B[3*num + idx] = a*det_inv;
}

inline void predictState(uniform KalmanFilter * uniform filter, uniform KalmanFilterParameter * uniform param, uniform int hit, int idx){
  gemm(0, 0, param->hits[hit].Jacobi, filter->P_k1, filter->P_k, 5, 1, 5, filter->ntracks, idx);
}

inline void predictCovariance(uniform KalmanFilter * uniform filter, uniform KalmanFilterParameter * uniform param, uniform int hit, int idx){
  gemm(0, 0, param->hits[hit].Jacobi, filter->C_k1, filter->M55_tmp, 5, 5, 5, filter->ntracks, idx);
  gemm(0, 1, filter->M55_tmp, param->hits[hit].Jacobi, filter->C_k, 5, 5, 5, filter->ntracks, idx);
}

inline void correctGain1D(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim1_H, filter->C_k, filter->M15_tmp, 1, 5, 5, filter->ntracks, idx);
  gemm(0, 1, filter->M15_tmp, filter->Dim1_H, filter->M11_tmp, 1, 1, 5, filter->ntracks, idx);
  add(filter->M11_tmp, filter->Dim1_V, filter->M11_tmp, 1, 1, filter->ntracks, idx);

  gemm(0, 1, filter->C_k, filter->Dim1_H, filter->Dim1_K, 5, 1, 5, filter->ntracks, idx);
  divide(filter->Dim1_K, filter->M11_tmp, filter->Dim1_K, 5, 1, filter->ntracks, idx);
}

inline void correctGain2D(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim2_H, filter->C_k, filter->M25_tmp, 2, 5, 5, filter->ntracks, idx);
  gemm(0, 1, filter->M25_tmp, filter->Dim2_H, filter->M22_tmp, 2, 2, 5, filter->ntracks, idx);
  add(filter->M22_tmp, filter->Dim2_V, filter->M22_tmp, 2, 2, filter->ntracks, idx);
  inverse2x2(filter->M22_tmp, filter->M22_tmp, filter->ntracks, idx);
  gemm(1, 0, filter->Dim2_H, filter->M22_tmp, filter->M25_tmp, 5, 2, 2, filter->ntracks, idx); //check better...
  gemm(0, 0, filter->C_k, filter->M25_tmp, filter->Dim2_K, 5, 2, 5, filter->ntracks, idx);
}

inline void correctState(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim1_H, filter->P_k, filter->M11_tmp, 1, 1, 5, filter->ntracks, idx);
  subtract(filter->Dim1_m, filter->M11_tmp, filter->M11_tmp, 1, 1, filter->ntracks, idx);
  multiply(filter->Dim1_K, filter->M11_tmp, filter->M15_tmp, 5, 1, filter->ntracks, idx);
  add(filter->P_k, filter->M15_tmp, filter->P_k1, 5, 1, filter->ntracks, idx);
}

inline void correctState2D(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim2_H, filter->P_k, filter->v2_tmp, 2, 1, 5, filter->ntracks, idx);
  subtract(filter->Dim2_m, filter->v2_tmp, filter->v2_tmp, 2, 1, filter->ntracks, idx);
  gemm(0, 0, filter->Dim2_K, filter->v2_tmp, filter->v5_tmp, 5, 1, 2, filter->ntracks, idx);
  add(filter->v5_tmp, filter->P_k, filter->P_k1, 5, 1, filter->ntracks, idx);
}

inline void correctCovariance(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim1_K, filter->Dim1_H, filter->M55_tmp, 5, 5, 1, filter->ntracks, idx);
  identity(filter->M55_tmp1, 5, filter->ntracks, idx);
  subtract(filter->M55_tmp1, filter->M55_tmp, filter->M55_tmp, 5, 5, filter->ntracks, idx);
  gemm(0, 0, filter->M55_tmp, filter->C_k, filter->C_k1, 5, 5, 5, filter->ntracks, idx);
}

inline void correctCovariance2D(uniform KalmanFilter * uniform filter, int idx){
  gemm(0, 0, filter->Dim2_K, filter->Dim2_H, filter->M55_tmp, 5, 5, 2, filter->ntracks, idx);
  identity(filter->M55_tmp1, 5, filter->ntracks, idx);
  subtract(filter->M55_tmp1, filter->M55_tmp, filter->M55_tmp, 5, 5, filter->ntracks, idx);
  gemm(0, 0, filter->M55_tmp, filter->C_k, filter->C_k1, 5, 5, 5, filter->ntracks, idx);
}

inline void kalman_filter(uniform KalmanFilter * uniform filter, uniform KalmanFilterParameter * uniform param, int i){
  M5(P_k1, 0, 0, i) = 0.;
  M5(P_k1, 0, 1, i) = 0.;
  M5(P_k1, 0, 2, i) = 0.;
  M5(P_k1, 0, 3, i) = 0.;
  M5(P_k1, 0, 4, i) = 0.;

  // Initial error covariance matrix estimate
  M5(C_k1, 0, 0, i) = 250.;
  M5(C_k1, 1, 1, i) = 250.;
  M5(C_k1, 2, 2, i) = 0.25;
  M5(C_k1, 3, 3, i) = 0.25;
  M5(C_k1, 4, 4, i) = 0.000001;

  // Set the transport matrix
  M5(Dim1_H, 0, 0, i) = 1.0;
  M5(Dim2_H, 0, 0, i) = 1.0;
  M5(Dim2_H, 1, 1, i) = 1.0;

  for(uniform int h = 0; h < param->max_hit_count; h++){
    if(h >= param->hit_count[i])
      continue;

    predictState(filter, param, h, i);
    predictCovariance(filter, param, h, i);

    if(param->hits[h].is2Dim[i]){
      V(Dim2_m, 0, i) = param->hits[h].normal[i] - param->hits[h].ref[i];
      V(Dim2_m, 1, i) = param->hits[h].normal[filter->ntracks + i] - param->hits[h].ref[filter->ntracks + i];

      M2(Dim2_V, 0, 0, i) = param->hits[h].err_locX[i]*param->hits[h].err_locX[i];
      M2(Dim2_V, 0, 1, i) = param->hits[h].cov_locXY[i];
      M2(Dim2_V, 1, 0, i) = param->hits[h].cov_locXY[h];
      M2(Dim2_V, 1, 1, i) = param->hits[h].err_locY[i]*param->hits[h].err_locY[i];

      correctGain2D(filter, i);
      correctState2D(filter, i);
      correctCovariance2D(filter, i);
    }else{
      V(Dim1_m, 0, i) = param->hits[h].normal[i] - param->hits[h].ref[i];
      M1(Dim1_V, 0, i) = param->hits[h].err_locX[i]*param->hits[h].err_locX[i];

      correctGain1D(filter, i);
      correctState(filter, i);
      correctCovariance(filter, i);
    }
  }
}

export void startFilter(uniform KalmanFilter * uniform filter, uniform KalmanFilterParameter * uniform param){
  foreach(i = 0 ... filter->ntracks){
    kalman_filter(filter, param, i);
  }
}

export void foo(uniform TrackHits * uniform, uniform KalmanFilterParameter * uniform){} // Needed to generate the C++ version of TrackHit


// Test Routines
export void identity_test(uniform scalar_t * uniform A,  uniform int M, uniform int num){
  foreach(i = 0 ... num)
    identity(A, M, num, i);
}

export void divide_test(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num){
  foreach(i = 0 ... num)
    divide(A, B, C, M, N, num, i);
}

export void multiply_test(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num){
  foreach(i = 0 ... num)
    multiply(A, B, C, M, N, num, i);
}

export void subtract_test(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num){
  foreach(i = 0 ... num)
    subtract(A, B, C, M, N, num, i);
}

export void add_test(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int num){
  foreach(i = 0 ... num)
    add(A, B, C, M, N, num, i);
}

export void gemm_test(uniform int transpose_A, uniform int transpose_B, uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform scalar_t * uniform C, uniform int M, uniform int N, uniform int K, uniform int num){
  if(N == 1){ // Ugly but generates more efficient code for matrix/vector case...
    foreach(i = 0 ... num){
      gemm(transpose_A, transpose_B, A, B, C, M, 1, K, num, i);
    }
  }else{
    foreach(i = 0 ... num){
      gemm(transpose_A, transpose_B, A, B, C, M, N, K, num, i);
    }
  }
}

export void inverse2x2_test(uniform scalar_t * uniform A, uniform scalar_t * uniform B, uniform int num){
  foreach(i = 0 ... num){
    inverse2x2(A, B, num, i);
  }
}
