#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include <cassert>
#include <vector>
#include <immintrin.h>
#include <string.h>

#include "KalmanFilter.ispc.h"
#include "serial/KalmanFilterSerial.h"

namespace ispc{

inline void KalmanFilterParameter_initialize(ispc::KalmanFilterParameter *param, int ntracks, int nhits){
  param->max_track_storage = ntracks;
  param->max_hit_storage = nhits;
  param->ntracks = ntracks;
  param->max_hit_count = nhits;
  param->hit_count = (int32_t *) _mm_malloc(sizeof(int32_t)*ntracks, ALIGNMENT);

  param->hits = (ispc::TrackHits *) _mm_malloc(sizeof(ispc::TrackHits)*nhits, ALIGNMENT);
  for(int i = 0; i < param->max_hit_count; i++){
    param->hits[i].Jacobi = (scalar_t *) _mm_malloc(sizeof(scalar_t)*5*5*ntracks, ALIGNMENT);
    param->hits[i].is2Dim = (int32_t *) _mm_malloc(sizeof(int32_t)*ntracks, ALIGNMENT);
    param->hits[i].normal = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);
    param->hits[i].ref = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);
    param->hits[i].err_locX = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
    param->hits[i].err_locY = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
    param->hits[i].cov_locXY = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
  }
}


inline void KalmanFilterParameter_initialize(ispc::KalmanFilterParameter *param, std::vector<Track_t> &tracks, int nhits){
  int ntracks = tracks.size();

  param->max_track_storage = ntracks;
  param->max_hit_storage = nhits;
  param->ntracks = ntracks;
  param->max_hit_count = nhits;
  param->hit_count = (int32_t *) _mm_malloc(sizeof(int32_t)*ntracks, ALIGNMENT);

  param->hits = (ispc::TrackHits *) _mm_malloc(sizeof(ispc::TrackHits)*nhits, ALIGNMENT);
  for(int i = 0; i < param->max_hit_count; i++){
    param->hits[i].Jacobi = (scalar_t *) _mm_malloc(sizeof(scalar_t)*5*5*ntracks, ALIGNMENT);
    param->hits[i].is2Dim = (int32_t *) _mm_malloc(sizeof(int32_t)*ntracks, ALIGNMENT);
    param->hits[i].normal = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);
    param->hits[i].ref = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);
    param->hits[i].err_locX = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
    param->hits[i].err_locY = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
    param->hits[i].cov_locXY = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);

    for(int t = 0; t < ntracks; t++){
      Track_t &track = tracks[t];

      if(i >= track.track.size())
	continue;

      TrackHit_t &hit = track.track.at(i);

      param->hit_count[t] = track.track.size();
      param->hits[i].is2Dim[0*ntracks + t] = hit.is2Dim;
      param->hits[i].normal[0*ntracks + t] = hit.normal[0];
      param->hits[i].normal[1*ntracks + t] = hit.normal[1];
      param->hits[i].ref[0*ntracks + t] = hit.ref[0];
      param->hits[i].ref[1*ntracks + t] = hit.ref[1];
      param->hits[i].err_locX[t] = hit.err_locX;
      param->hits[i].err_locY[t] = hit.err_locY;
      param->hits[i].cov_locXY[t] = hit.cov_locXY;
      param->hits[i].Jacobi[0*ntracks + t] = hit.jacobi[0];
      param->hits[i].Jacobi[1*ntracks + t] = hit.jacobi[1];
      param->hits[i].Jacobi[2*ntracks + t] = hit.jacobi[2];
      param->hits[i].Jacobi[3*ntracks + t] = hit.jacobi[3];
      param->hits[i].Jacobi[4*ntracks + t] = hit.jacobi[4];
      param->hits[i].Jacobi[5*ntracks + t] = hit.jacobi[5];
      param->hits[i].Jacobi[6*ntracks + t] = hit.jacobi[6];
      param->hits[i].Jacobi[7*ntracks + t] = hit.jacobi[7];
      param->hits[i].Jacobi[8*ntracks + t] = hit.jacobi[8];
      param->hits[i].Jacobi[9*ntracks + t] = hit.jacobi[9];
      param->hits[i].Jacobi[10*ntracks + t] = hit.jacobi[10];
      param->hits[i].Jacobi[11*ntracks + t] = hit.jacobi[11];
      param->hits[i].Jacobi[12*ntracks + t] = hit.jacobi[12];
      param->hits[i].Jacobi[13*ntracks + t] = hit.jacobi[13];
      param->hits[i].Jacobi[14*ntracks + t] = hit.jacobi[14];
      param->hits[i].Jacobi[15*ntracks + t] = hit.jacobi[15];
      param->hits[i].Jacobi[16*ntracks + t] = hit.jacobi[16];
      param->hits[i].Jacobi[17*ntracks + t] = hit.jacobi[17];
      param->hits[i].Jacobi[18*ntracks + t] = hit.jacobi[18];
      param->hits[i].Jacobi[19*ntracks + t] = hit.jacobi[19];
      param->hits[i].Jacobi[20*ntracks + t] = hit.jacobi[20];
      param->hits[i].Jacobi[21*ntracks + t] = hit.jacobi[21];
      param->hits[i].Jacobi[22*ntracks + t] = hit.jacobi[22];
      param->hits[i].Jacobi[23*ntracks + t] = hit.jacobi[23];
      param->hits[i].Jacobi[24*ntracks + t] = hit.jacobi[24];
    }
  }
}

inline void KalmanFilterParameter_deallocate(KalmanFilterParameter *param){
  _mm_free(param->hit_count);

  for(int i = 0; i < param->max_hit_count; i++){
    _mm_free(param->hits[i].Jacobi);
    _mm_free(param->hits[i].is2Dim);
    _mm_free(param->hits[i].normal);
    _mm_free(param->hits[i].ref);
    _mm_free(param->hits[i].err_locX);
    _mm_free(param->hits[i].err_locY);
    _mm_free(param->hits[i].cov_locXY);
  }
}

inline void KalmanFilterParameter_reinitialize(ispc::KalmanFilterParameter *param, std::vector<Track_t> &tracks, int nhits){
  if(param->max_hit_storage < nhits || param->max_track_storage < tracks.size()){
    KalmanFilterParameter_deallocate(param);
    KalmanFilterParameter_initialize(param, tracks, nhits);
  }else{
    int ntracks = tracks.size();

    param->ntracks = ntracks;
    param->max_hit_count = nhits;

    for(int i = 0; i < param->max_hit_count; i++){
      for(int t = 0; t < ntracks; t++){
	Track_t &track = tracks[t];

	if(i >= track.track.size())
	  continue;

	TrackHit_t &hit = track.track.at(i);

	param->hit_count[t] = track.track.size();
	param->hits[i].is2Dim[0*ntracks + t] = hit.is2Dim;
	param->hits[i].normal[0*ntracks + t] = hit.normal[0];
	param->hits[i].normal[1*ntracks + t] = hit.normal[1];
	param->hits[i].ref[0*ntracks + t] = hit.ref[0];
	param->hits[i].ref[1*ntracks + t] = hit.ref[1];
	param->hits[i].err_locX[t] = hit.err_locX;
	param->hits[i].err_locY[t] = hit.err_locY;
	param->hits[i].cov_locXY[t] = hit.cov_locXY;
	param->hits[i].Jacobi[0*ntracks + t] = hit.jacobi[0];
	param->hits[i].Jacobi[1*ntracks + t] = hit.jacobi[1];
	param->hits[i].Jacobi[2*ntracks + t] = hit.jacobi[2];
	param->hits[i].Jacobi[3*ntracks + t] = hit.jacobi[3];
	param->hits[i].Jacobi[4*ntracks + t] = hit.jacobi[4];
	param->hits[i].Jacobi[5*ntracks + t] = hit.jacobi[5];
	param->hits[i].Jacobi[6*ntracks + t] = hit.jacobi[6];
	param->hits[i].Jacobi[7*ntracks + t] = hit.jacobi[7];
	param->hits[i].Jacobi[8*ntracks + t] = hit.jacobi[8];
	param->hits[i].Jacobi[9*ntracks + t] = hit.jacobi[9];
	param->hits[i].Jacobi[10*ntracks + t] = hit.jacobi[10];
	param->hits[i].Jacobi[11*ntracks + t] = hit.jacobi[11];
	param->hits[i].Jacobi[12*ntracks + t] = hit.jacobi[12];
	param->hits[i].Jacobi[13*ntracks + t] = hit.jacobi[13];
	param->hits[i].Jacobi[14*ntracks + t] = hit.jacobi[14];
	param->hits[i].Jacobi[15*ntracks + t] = hit.jacobi[15];
	param->hits[i].Jacobi[16*ntracks + t] = hit.jacobi[16];
	param->hits[i].Jacobi[17*ntracks + t] = hit.jacobi[17];
	param->hits[i].Jacobi[18*ntracks + t] = hit.jacobi[18];
	param->hits[i].Jacobi[19*ntracks + t] = hit.jacobi[19];
	param->hits[i].Jacobi[20*ntracks + t] = hit.jacobi[20];
	param->hits[i].Jacobi[21*ntracks + t] = hit.jacobi[21];
	param->hits[i].Jacobi[22*ntracks + t] = hit.jacobi[22];
	param->hits[i].Jacobi[23*ntracks + t] = hit.jacobi[23];
	param->hits[i].Jacobi[24*ntracks + t] = hit.jacobi[24];
      }
    }
  }
}

inline void KalmanFilter_initialize(KalmanFilter *filter, int ntracks){
  filter->C_k1 = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);
  filter->C_k = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);

  filter->Dim1_H = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);
  filter->Dim2_H = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*2, ALIGNMENT);

  filter->Dim1_K = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);
  filter->Dim2_K = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*2, ALIGNMENT);
  filter->Dim5_K = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);

  filter->Dim1_V = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
  filter->Dim2_V = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*4, ALIGNMENT);

  filter->Dim1_m = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
  filter->Dim2_m = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);

  filter->P_k1 = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);
  filter->P_k = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);

  filter->Inverse = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);

  filter->M55_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);
  filter->M15_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);
  filter->M11_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks, ALIGNMENT);
  filter->M25_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2*5, ALIGNMENT);
  filter->M22_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2*2, ALIGNMENT);
  filter->M55_tmp1 = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5*5, ALIGNMENT);

  filter->v2_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*2, ALIGNMENT);
  filter->v5_tmp = (scalar_t *) _mm_malloc(sizeof(scalar_t)*ntracks*5, ALIGNMENT);

  memset(filter->C_k1, 0, sizeof(scalar_t)*ntracks*5*5);
  memset(filter->C_k, 0, sizeof(scalar_t)*ntracks*5*5);
  memset(filter->Dim1_H, 0, sizeof(scalar_t)*ntracks*5);
  memset(filter->Dim2_H, 0, sizeof(scalar_t)*ntracks*5*2);
  memset(filter->Dim1_K, 0, sizeof(scalar_t)*ntracks*5);
  memset(filter->Dim2_K, 0, sizeof(scalar_t)*ntracks*5*2);
  memset(filter->Dim5_K, 0, sizeof(scalar_t)*ntracks*5*5);
  memset(filter->Dim1_V, 0, sizeof(scalar_t)*ntracks);
  memset(filter->Dim2_V, 0, sizeof(scalar_t)*ntracks*4);
  memset(filter->Dim1_m, 0, sizeof(scalar_t)*ntracks);
  memset(filter->Dim2_m, 0, sizeof(scalar_t)*ntracks*2);
  memset(filter->P_k1, 0, sizeof(scalar_t)*ntracks*5);
  memset(filter->P_k, 0, sizeof(scalar_t)*ntracks*5);

  filter->ntracks = ntracks;
  filter->max_track_storage = ntracks;
}


inline void KalmanFilter_deallocate(ispc::KalmanFilter *filter){
  _mm_free(filter->C_k1);
  _mm_free(filter->C_k);
  _mm_free(filter->Dim1_H);
  _mm_free(filter->Dim2_H);
  _mm_free(filter->Dim1_K);
  _mm_free(filter->Dim2_K);
  _mm_free(filter->Dim5_K);
  _mm_free(filter->Dim1_V);
  _mm_free(filter->Dim2_V);
  _mm_free(filter->Dim1_m);
  _mm_free(filter->Dim2_m);
  _mm_free(filter->P_k1);
  _mm_free(filter->P_k);
  _mm_free(filter->Inverse);
  _mm_free(filter->M55_tmp);
  _mm_free(filter->M15_tmp);
  _mm_free(filter->M11_tmp);
  _mm_free(filter->M55_tmp1);
  _mm_free(filter->M25_tmp);
  _mm_free(filter->v2_tmp);
  _mm_free(filter->v5_tmp);
}

inline void KalmanFilter_reinitialize(KalmanFilter *filter, int ntracks){
  if(filter->max_track_storage < ntracks){
    KalmanFilter_deallocate(filter);
    KalmanFilter_initialize(filter, ntracks);
  }else{
    filter->ntracks = ntracks;
  }
}

}

#endif
