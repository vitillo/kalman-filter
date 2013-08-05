#ifndef COMMON_H_
#define COMMON_H_

#include "../common.h"

// 0: none 1: Results only 2:reserved 3:everything
#define DBG_LVL 0
//#define DEBUG

#ifndef HAS_COMPETING_ROTS
#define HAS_COMPETING_ROTS 0
#endif
#ifndef HAS_MATERIAL_EFFECTS
#define HAS_MATERIAL_EFFECTS 1
#endif

#ifndef CUDA_PROJECT
#ifndef USE_DOUBLE_CALCULATIONS
#define USE_DOUBLE_CALCULATIONS 1
#endif
#ifndef BENCHMARK
#define BENCHMARK 1
#endif
#if USE_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#if USE_DOUBLE_CALCULATIONS
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double gpu_scalar_t;
#else
typedef float gpu_scalar_t;
#endif
#ifndef BUILDSERVER
#define BUILDSERVER 0
#endif

#if BUILDSERVER || defined(__APPLE__)
#define SELECTED_DEVICETYPE CL_DEVICE_TYPE_CPU
#define TARGET_DEVICE_NAME "Intel(R)"
#else
#define SELECTED_DEVICETYPE CL_DEVICE_TYPE_GPU
#define TARGET_DEVICE_NAME "NVIDIA"
#endif
#else
typedef double scalar_t;
#endif

#define ORDER 5
#define MAXROTS 30

#define PionMassSquared 19479. //139.57*139.57

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1
#endif

struct trackHitDataStruct {
  scalar_t locX;
  scalar_t locY;
  scalar_t Phi;
  scalar_t Theta;
  scalar_t QoverP;
};

struct trackHitDataStructArray {
  //struct trackHitDataStruct member;
  scalar_t data[5];
};

typedef struct trackHitDataStructArray trackHitData;

struct TrackInfoStruct {
  scalar_t d0;
  scalar_t z0;
  scalar_t phi;
  scalar_t theta;
  scalar_t qoverp;
};

struct TrackHitStruct {
  scalar_t normal[ORDER]; // needs to be only 2-dim
  scalar_t ref[ORDER];    // does not need to be transfered
  scalar_t err_locX;
  scalar_t err_locY;
  scalar_t cov_locXY;
#if HAS_MATERIAL_EFFECTS
  scalar_t sigmaDeltaThetaSquared;
  scalar_t sigmaDeltaPhiSquared;
  scalar_t sigmaDeltaQoverPSquared;
#endif
  scalar_t jacobi[ORDER * ORDER];  // need to transfer only upper triangle (15 out of 25)
  scalar_t jacobiInverse[ORDER * ORDER];  // does not need to be transfered
  char is2Dim;
  int detType;// does not need to be transfered
  int bec;// does not need to be transfered
#if HAS_COMPETING_ROTS
  unsigned int competingRotIds;// does not need to be transfered
  scalar_t rotLocX[MAXROTS ];
  scalar_t rotLocY[MAXROTS ];
  int rotIds[MAXROTS];

  scalar_t rotProbs[MAXROTS];// does not need to be transfered
  scalar_t rotCov[MAXROTS][4]; // [00;01;10;11]
#endif
};

typedef struct TrackHitStruct TrackHit_t;
typedef struct TrackInfoStruct TrackInfo_t;

//Rene FilterParams
struct kalmanFilterParamStruct {
  TrackHit_t *hits;
  scalar_t *C_k_1k_1_global;
  int *hitCount;
  trackHitData *fitsForward;
  trackHitData *fitsBackward;
  trackHitData *fitsResult;
  scalar_t *C_Values;
  scalar_t *C_Inverse;
  scalar_t *C_Result;
  TrackInfo_t* fittedPerigee;
  scalar_t *C_fittedPerigee;
};

typedef struct kalmanFilterParamStruct KalmanFilterParameter_t;

#endif /* COMMON_H_ */
