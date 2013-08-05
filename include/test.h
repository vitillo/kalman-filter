#ifndef TEST_H
#define TEST_H

#include <vector>

#include "ROOT_Ttree_io/event_reader.h"
#include "ispc/KalmanFilter.h"

using namespace std;

void KalmanFilterParallelNTracksMHitsMixed(ispc::KalmanFilterParameter *param, int ntracks, int nhits){
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
      param->hit_count[t] = nhits;
      param->hits[i].is2Dim[0*ntracks + t] = (t & 1);
      param->hits[i].normal[0*ntracks + t] = -23.9713;
      param->hits[i].normal[1*ntracks + t] = -13.9713;
      param->hits[i].ref[0*ntracks + t] = -23.9716;
      param->hits[i].ref[1*ntracks + t] = -13.9716;
      param->hits[i].err_locX[t] = 0.0140155;
      param->hits[i].err_locY[t] = 0.0240155;
      param->hits[i].cov_locXY[t] = 0.0040155;
      param->hits[i].Jacobi[0*ntracks + t] = 1.09201;
      param->hits[i].Jacobi[1*ntracks + t] = 3.88707e-18;
      param->hits[i].Jacobi[2*ntracks + t] = 330.602;
      param->hits[i].Jacobi[3*ntracks + t] = 33.952;
      param->hits[i].Jacobi[4*ntracks + t] = -45077.2;
      param->hits[i].Jacobi[5*ntracks + t] = 0.834501;
      param->hits[i].Jacobi[6*ntracks + t] = 1;
      param->hits[i].Jacobi[7*ntracks + t] = 118.912;
      param->hits[i].Jacobi[8*ntracks + t] = -772;
      param->hits[i].Jacobi[9*ntracks + t] = -28567.4;
      param->hits[i].Jacobi[10*ntracks + t] = -0.000351449;
      param->hits[i].Jacobi[11*ntracks + t] = -5.44478e-22;
      param->hits[i].Jacobi[12*ntracks + t] = 0.949725;
      param->hits[i].Jacobi[13*ntracks + t] = 0.202329;
      param->hits[i].Jacobi[14*ntracks + t] = -270.811;
      param->hits[i].Jacobi[15*ntracks + t] = 3.11037e-07;
      param->hits[i].Jacobi[16*ntracks + t] = 4.8187e-25;
      param->hits[i].Jacobi[17*ntracks + t] =  0.000631694;
      param->hits[i].Jacobi[18*ntracks + t] =  0.999819;
      param->hits[i].Jacobi[19*ntracks + t] = 0.102849;
      param->hits[i].Jacobi[20*ntracks + t] = 0;
      param->hits[i].Jacobi[21*ntracks + t] = 0;
      param->hits[i].Jacobi[22*ntracks + t] = 0;
      param->hits[i].Jacobi[23*ntracks + t] = 0;
      param->hits[i].Jacobi[24*ntracks + t] = 1;
    }
  }
}

void initializeTrack(Track_t &track, int nhits, bool is2D){
  TrackHit_t hit;

  hit.is2Dim = is2D;
  hit.normal[0] = -23.9713;
  hit.normal[1] = -13.9713;
  hit.ref[0] = -23.9716;
  hit.ref[1] = -13.9716;
  hit.err_locX = 0.0140155;
  hit.err_locY = 0.0240155;
  hit.cov_locXY = 0.040155;
  hit.jacobi[0] = 1.09201;
  hit.jacobi[1] = 3.88707e-18;
  hit.jacobi[2] = 330.602;
  hit.jacobi[3] = 33.952;
  hit.jacobi[4] = -45077.2;
  hit.jacobi[5] = 0.834501;
  hit.jacobi[6] = 1;
  hit.jacobi[7] = 118.912;
  hit.jacobi[8] = -772;
  hit.jacobi[9] = -28567.4;
  hit.jacobi[10] = -0.000351449;
  hit.jacobi[11] = -5.44478e-22;
  hit.jacobi[12] = 0.949725;
  hit.jacobi[13] = 0.202329;
  hit.jacobi[14] = -270.811;
  hit.jacobi[15] = 3.11037e-07;
  hit.jacobi[16] = 4.8187e-25;
  hit.jacobi[17] =  0.000631694;
  hit.jacobi[18] =  0.999819;
  hit.jacobi[19] = 0.102849;
  hit.jacobi[20] = 0;
  hit.jacobi[21] = 0;
  hit.jacobi[22] = 0;
  hit.jacobi[23] = 0;
  hit.jacobi[24] = 1;

  for(int i = 0; i < nhits; i++)
    track.track.push_back(hit);
}

#endif
