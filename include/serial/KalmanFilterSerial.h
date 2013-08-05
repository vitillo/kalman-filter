//==================================================================================================
//
// Name         : KalmanFilterSerial.h
// Description  : Simple serial Kalman Filter class.
//
//==================================================================================================

#ifndef H_GUARD_KALMAN_FILTER_SERIAL
#define H_GUARD_KALMAN_FILTER_SERIAL

#include <unistd.h>
#include <cmath>

#include "ROOT_Ttree_io/event_reader.h"
#include "serial/common_gsl.h"

#define PRINT_ERR(x) \
  std::cerr << "ERROR: " << __FILE__ << ", line: " << __LINE__ << " " << x << std::endl;

#define PRINT_WARN(x) \
  std::cerr << "WARNING: " << __FILE__ << ", line: " << __LINE__ << " " << x << std::endl;

class KalmanFilterSerial {
  public:
    KalmanFilterSerial  (void);
    ~KalmanFilterSerial (void);

    Track_t filteredTrack (void);

    int setTrackData (Track_t track);
    int startFilter ();

  private:
    int predictCovariance (GSL_MATRIX *dest, GSL_MATRIX *F_k, GSL_MATRIX *C_k1) const;
    int predictState (GSL_VECTOR *dest, GSL_MATRIX *F_k, GSL_VECTOR *P_k1) const;
    int correctCovariance (GSL_MATRIX *dest, GSL_MATRIX *K_k, GSL_MATRIX *H_k, GSL_MATRIX *C_k1) const;
    int correctGain (GSL_MATRIX *dest, GSL_MATRIX *C_k1, GSL_MATRIX *H_k, GSL_MATRIX *V_k) const;
    int correctState (GSL_VECTOR *dest, GSL_VECTOR *P_k1, GSL_MATRIX *K_k, GSL_VECTOR *m_k, GSL_MATRIX *H_k) const;

    Track_t mTrack;
    Track_t mFilteredTrack;
    Track_t mFilteredTrackBackwards;

    GSL_MATRIX * mC_k1;
    GSL_MATRIX * mC_k;

    GSL_MATRIX * mDim1_H;
    GSL_MATRIX * mDim2_H;

    GSL_MATRIX * mDim1_K;
    GSL_MATRIX * mDim2_K;
    GSL_MATRIX * mDim5_K;

    GSL_MATRIX * mDim1_V;
    GSL_MATRIX * mDim2_V;

    GSL_VECTOR * mDim1_m;
    GSL_VECTOR * mDim2_m;

    GSL_VECTOR * mP_k1;
    GSL_VECTOR * mP_k;

    GSL_MATRIX * mInverse;

    friend void testKalmanFilterNTracksMHitsMixed(int ntracks, int hits);
};

#endif
