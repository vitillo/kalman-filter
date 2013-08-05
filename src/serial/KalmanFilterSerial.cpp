//==================================================================================================
//
// Name         : KalmanFilterSerial.h
// Description  : Simple serial Kalman Filter class.
//
//==================================================================================================

#include "serial/KalmanFilterSerial.h"

// Constructor
KalmanFilterSerial::KalmanFilterSerial (void)
{
  // Allocate vectors and matrices. Using fixed values, since they won't change
  // anyway and are easier to read.
  mC_k1   = GSL_MATRIX_CALLOC (5, 5);
  mC_k    = GSL_MATRIX_CALLOC (5, 5);

  mDim1_H = GSL_MATRIX_CALLOC (1, 5);
  mDim2_H = GSL_MATRIX_CALLOC (2, 5);

  mDim1_K = GSL_MATRIX_CALLOC (5, 1);
  mDim2_K = GSL_MATRIX_CALLOC (5, 2);
  mDim5_K = GSL_MATRIX_CALLOC (5, 5);

  mDim1_V = GSL_MATRIX_CALLOC (1, 1);
  mDim2_V = GSL_MATRIX_CALLOC (2, 2);

  mDim1_m = GSL_VECTOR_CALLOC (1);
  mDim2_m = GSL_VECTOR_CALLOC (2);

  mP_k1   = GSL_VECTOR_CALLOC (5);
  mP_k    = GSL_VECTOR_CALLOC (5);

  mInverse = GSL_MATRIX_CALLOC (5, 5);
}

// Destructor
KalmanFilterSerial::~KalmanFilterSerial (void)
{
  // Free the resources
  GSL_MATRIX_FREE (mC_k1);
  GSL_MATRIX_FREE (mC_k);
  GSL_MATRIX_FREE (mDim1_H);
  GSL_MATRIX_FREE (mDim2_H);
  GSL_MATRIX_FREE (mDim1_K);
  GSL_MATRIX_FREE (mDim2_K);
  GSL_MATRIX_FREE (mDim1_V);
  GSL_MATRIX_FREE (mDim2_V);
  GSL_MATRIX_FREE (mInverse);

  GSL_VECTOR_FREE (mDim1_m);
  GSL_VECTOR_FREE (mDim2_m);
  GSL_VECTOR_FREE (mP_k1);
  GSL_VECTOR_FREE (mP_k);
}

// Set the track data
  int
KalmanFilterSerial::setTrackData (Track_t track)
{
  if (track.track.empty ( ))
    return -1;

  mTrack = track;

  return 0;
}

// Return the filtered track
  Track_t
KalmanFilterSerial::filteredTrack (void)
{
  return mFilteredTrack;
}

// Predict the error covariance ahead using the formula:
// C_(k|k-1) = F_(k) C_(k|k-1) F_(k)^T [+ P_(k) Q_(k) P_(k)^T]
// The term in [ ] brackets is ommited right now..
int
KalmanFilterSerial::predictCovariance (GSL_MATRIX *dest, GSL_MATRIX *F_k, GSL_MATRIX *C_k1) const
{
  if (!dest || !F_k || !C_k1) {
    PRINT_WARN ("Parameter is NULL!");
    return -1;
  }

#if DBG_LVL > 2
  std::cout << "predictCovariance:" << std::endl;
  std::cout << "F_k: ";
  printMatrix (F_k);
  std::cout << "C_k: ";
  printMatrix (C_k1);
#endif

  GSL_MATRIX *tmp;
  tmp = GSL_MATRIX_CALLOC (F_k->size1, F_k->size2);

  // tmp = F_(k) * C_(k|k-1)
  GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0f, F_k, C_k1, 0.0f, tmp);

  // dest = tmp * F_(k)^T
  GSL_BLAS_XGEMM (CblasNoTrans, CblasTrans, 1.0f, tmp, F_k, 0.0f, dest);

  GSL_MATRIX_FREE (tmp);

#if DBG_LVL > 2
  std::cout << "maxima predict Covariance:" << std::endl;
  std::cout << "F_k . C_k . transpose(F_k);" << std::endl;
  std::cout << "res: ";
  printMatrix (dest);
#endif

  return 0;
}

// Predict the state ahead using the formula:
// P_(k|k-1) = F_(k) P_(k|k-1)
int
KalmanFilterSerial::predictState (GSL_VECTOR *dest, GSL_MATRIX *F_k, GSL_VECTOR *P_k1) const
{
  if (!dest || !F_k || !P_k1) {
    PRINT_WARN ("Parameter is NULL!");
    return -1;
  }

#if DBG_LVL > 2
  std::cout << "predictState:" << std::endl;
  std::cout << "F_k: ";
  printMatrix (F_k);
  std::cout << "P_k: ";
  printVector (P_k1);
#endif

  GSL_BLAS_XGEMV (CblasNoTrans, 1.0f, F_k, P_k1, 0.0f, dest);

#if DBG_LVL > 2
  std::cout << "maxima predict State:" << std::endl;
  std::cout << "F_k . P_k;" << std::endl;
  std::cout << "res: ";
  printVector (dest);
#endif

  return 0;
}

// Update the error covariance using the formula:
// C_(k|k) = (I - K_(k) H_(k)) C_(k|k-1)
int
KalmanFilterSerial::correctCovariance (GSL_MATRIX *dest, GSL_MATRIX *K_k,
    GSL_MATRIX *H_k, GSL_MATRIX *C_k1) const
{
  if (!dest || !K_k || !H_k || !C_k1) {
    PRINT_WARN ("Parameter is NULL!");
    return -1;
  }

#if DBG_LVL > 2
  std::cout << "correctCovariance:" << std::endl;
  std::cout << "H_k: ";
  printMatrix (H_k);
  std::cout << "K_k: ";
  printMatrix (K_k);
  std::cout << "C_k: ";
  printMatrix (C_k1);
#endif

  GSL_MATRIX *tmp, *identity;

  tmp      = GSL_MATRIX_CALLOC (K_k->size1, K_k->size1);
  identity = GSL_MATRIX_CALLOC (K_k->size1, K_k->size1);

  GSL_MATRIX_SET_IDENTITY (identity);

  // tmp = K_(k) * H_(k)
  GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0, K_k, H_k, 0.0, tmp);

#if DBG_LVL > 2
  std::cout << "tmp: ";
  printMatrix (tmp);
  std::cout << "identity: ";
  printMatrix (identity);
#endif

  // identity = identity - tmp
  GSL_MATRIX_SUB (identity, tmp);

  // dest = identity * C_(k|k-1)
  GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0, identity, C_k1, 0.0, dest);

  GSL_MATRIX_FREE (identity);
  GSL_MATRIX_FREE (tmp);

#if DBG_LVL > 2
  std::cout << "maxima correct Covariance:" << std::endl;
  std::cout << "(ident (5) - K_k . H_k) . C_k;" << std::endl;
  std::cout << "res: ";
  printMatrix (dest);
#endif

  return 0;
}

// Compute Kalman gain using the formula:
// K_(k) = C_(k|k-1) H_(k)^T (V_(k) + H_(k) C_(k|k-1) H_(k)^T)^-1
int
KalmanFilterSerial::correctGain (GSL_MATRIX *dest, GSL_MATRIX *C_k1,
    GSL_MATRIX *H_k,  GSL_MATRIX *V_k) const
{
  if (!dest || !C_k1 || !H_k || !V_k) {
    PRINT_WARN ("Parameter is NULL!");
    return -1;
  }

#if DBG_LVL > 2
  std::cout << "correctGain:" << std::endl;
  std::cout << "H_k: ";
  printMatrix (H_k);
  std::cout << "C_k: ";
  printMatrix (C_k1);
  std::cout << "V_k: ";
  printMatrix (V_k);
#endif

  GSL_MATRIX *tmp, *tmp2, *tmp3, *inverse;

  tmp     = GSL_MATRIX_CALLOC (H_k->size1, H_k->size2);
  tmp2    = GSL_MATRIX_CALLOC (H_k->size1, H_k->size1);
  tmp3    = GSL_MATRIX_CALLOC (H_k->size2, H_k->size1);
  inverse = GSL_MATRIX_CALLOC (H_k->size1, H_k->size1);

  GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0, H_k, C_k1, 0.0, tmp);
  GSL_BLAS_XGEMM (CblasNoTrans, CblasTrans,   1.0, tmp, H_k,  0.0, tmp2);

  GSL_MATRIX_ADD (tmp2, V_k);

  GSL_MATRIX_MEMCPY (inverse, tmp2);

  matrixInverseSymmetric (inverse);

  GSL_BLAS_XGEMM (CblasTrans,   CblasNoTrans, 1.0, H_k,  inverse, 0.0, tmp3);
  GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0, C_k1, tmp3,    0.0, dest);

  GSL_MATRIX_FREE (inverse);
  GSL_MATRIX_FREE (tmp);
  GSL_MATRIX_FREE (tmp2);
  GSL_MATRIX_FREE (tmp3);

#if DBG_LVL > 2
  std::cout << "maxima correct Gain:" << std::endl;
  std::cout << "C_k . transpose (H_k) . invert ((V_k + H_k . C_k . transpose (H_k)));" << std::endl;
  std::cout << "res: ";
  printMatrix (dest);
#endif

  return 0;
}

// Update the state prediction using the formula:
// p_(k|k) = p_(k|k-1) + K_k * (m_k - H_k * p_(k|k-1))
int
KalmanFilterSerial::correctState (GSL_VECTOR *dest, GSL_VECTOR *P_k1, GSL_MATRIX *K_k, 
    GSL_VECTOR *m_k,  GSL_MATRIX *H_k) const
{
  if (!dest || !P_k1 || !K_k || !m_k || !H_k) {
    PRINT_WARN ("Parameter is NULL!");
    return -1;
  }

#if DBG_LVL > 2
  std::cout << "correctState:" << std::endl;
  std::cout << "H_k: ";
  printMatrix (H_k);
  std::cout << "K_k: ";
  printMatrix (K_k);
  std::cout << "P_k: ";
  printVector (P_k1);
  std::cout << "m_k: ";
  printVector (m_k);
#endif

  GSL_VECTOR *tmp, *tmp2, *tmp3;

  tmp     = GSL_VECTOR_CALLOC (H_k->size1);
  tmp2    = GSL_VECTOR_CALLOC (H_k->size1);
  tmp3    = GSL_VECTOR_CALLOC (H_k->size2);

  GSL_BLAS_XGEMV (CblasNoTrans, 1.0, H_k, P_k1, 0.0, tmp);

  GSL_VECTOR_MEMCPY (tmp2, m_k);
  GSL_VECTOR_SUB    (tmp2, tmp);

  GSL_BLAS_XGEMV (CblasNoTrans, 1.0, K_k, tmp2, 0.0, tmp3);

  GSL_VECTOR_MEMCPY (dest, P_k1);
  GSL_VECTOR_ADD    (dest, tmp3);

#if DBG_LVL > 2
  std::cout << "maxima correct State:" << std::endl;
  std::cout << "res: ";
  printVector (dest);
#endif

  GSL_VECTOR_FREE (tmp);
  GSL_VECTOR_FREE (tmp2);
  GSL_VECTOR_FREE (tmp3);

  return 0;
}

// Start the Kalman Filter process.
//
// Variable naming scheme (DimX prefix referring to Xdimensional m_k):
// -------------------------------------------------------------------
// mC_k1 = Initial covariance estimate used for the prediction. Stores
//         the result of the measurement update to become the next
//         iterations prediction.
// mC_k  = The predicted covariance.
//
// mDim1_H = 1x5 measurement matrix to map from track parameter space
//           to measurement space.
// mDim2_H = 2x5 measurement matrix to map from track parameter space
//           to measurement space.
//
// mDim1_K = 5x1 Kalman gain matrix.
// mDim2_K = 5x2 Kalman gain matrix.
//
// mDim1_V = 1x1 variance matrix.
// mDim2_V = 2x2 variance matrix.
//
// mDim1_m = 1 dimensional measurement vector.
// mDim2_m = 2 dimensional measurement vector.
//
// mP_k1 = Initial state estimate used for the prediction. Stores
//         the result of the measurement update to become the next
//         iterations prediction.
// mP_k  = The predicted state.

  int
KalmanFilterSerial::startFilter()
{
  // Nothing to do...
  if (mTrack.track.empty ( )) {
    PRINT_WARN ("Empty Track1!");
    return -1;
  }

  KalmanFilterParameter_t kalmanParams;
  size_t hitCount = mTrack.track.size ( );

  kalmanParams.hits = (TrackHit_t*) calloc (hitCount, sizeof (TrackHit_t));
  kalmanParams.fitsForward  = (trackHitData*) calloc (hitCount, sizeof (trackHitData));
  kalmanParams.fitsBackward = (trackHitData*) calloc (hitCount, sizeof (trackHitData));
  kalmanParams.C_Values     = (scalar_t*) calloc (hitCount, ORDER * ORDER * sizeof (scalar_t));
  kalmanParams.C_Inverse    = (scalar_t*) calloc (hitCount, ORDER * ORDER * sizeof (scalar_t));

  for (size_t i = 0; i < hitCount; ++i) {
    kalmanParams.hits[i] = mTrack.track.at(i);
  }

  GSL_VECTOR_SET_ZERO (mP_k1);
  GSL_MATRIX_SET_ZERO (mC_k1);
  GSL_MATRIX_SET_ZERO (mDim2_H);

  // Initial error covariance matrix estimate
  GSL_MATRIX_SET (mC_k1, 0, 0, 250.);
  GSL_MATRIX_SET (mC_k1, 1, 1, 250.);
  GSL_MATRIX_SET (mC_k1, 2, 2, 0.25);
  GSL_MATRIX_SET (mC_k1, 3, 3, 0.25);
  GSL_MATRIX_SET (mC_k1, 4, 4, 0.000001);

  // Set the transport matrix H_(k)
  GSL_MATRIX_SET (mDim1_H, 0, 0, 1.0);

  GSL_MATRIX_SET (mDim2_H, 0, 0, 1.0);
  GSL_MATRIX_SET (mDim2_H, 1, 1, 1.0);

  for (size_t i = 0; i < hitCount; ++i) {
    GSL_MATRIX_VIEW jacobi;

    // Get the Transport-Jacobians
    jacobi = GSL_MATRIX_VIEW_ARRAY (kalmanParams.hits[i].jacobi, 5, 5);

    predictState (mP_k, &jacobi.matrix, mP_k1);
    predictCovariance (mC_k, &jacobi.matrix, mC_k1);

#if HAS_MATERIAL_EFFECTS
    scalar_t covVal = GSL_MATRIX_GET(mC_k, 2, 2);
    GSL_MATRIX_SET(mC_k, 2, 2, covVal+kalmanParams.hits[i].sigmaDeltaPhiSquared);
    covVal = GSL_MATRIX_GET(mC_k, 3, 3);
    GSL_MATRIX_SET(mC_k, 3, 3, covVal+kalmanParams.hits[i].sigmaDeltaThetaSquared);
    covVal = GSL_MATRIX_GET(mC_k, 4, 4);
    GSL_MATRIX_SET(mC_k, 4, 4, covVal+kalmanParams.hits[i].sigmaDeltaQoverPSquared);
#endif

    if (kalmanParams.hits[i].is2Dim) { 
      scalar_t value1, value2;
      value1 = kalmanParams.hits[i].normal[0] - kalmanParams.hits[i].ref[0];
      value2 = kalmanParams.hits[i].normal[1] - kalmanParams.hits[i].ref[1];

      GSL_VECTOR_SET (mDim2_m, 0, value1);
      GSL_VECTOR_SET (mDim2_m, 1, value2);

      GSL_MATRIX_SET (mDim2_V, 0, 0, kalmanParams.hits[i].err_locX * kalmanParams.hits[i].err_locX);
      GSL_MATRIX_SET (mDim2_V, 0, 1, kalmanParams.hits[i].cov_locXY);
      GSL_MATRIX_SET (mDim2_V, 1, 0, kalmanParams.hits[i].cov_locXY);
      GSL_MATRIX_SET (mDim2_V, 1, 1, kalmanParams.hits[i].err_locY * kalmanParams.hits[i].err_locY);

      correctGain  (mDim2_K, mC_k, mDim2_H, mDim2_V);
      correctState (mP_k1, mP_k, mDim2_K, mDim2_m, mDim2_H);
      correctCovariance (mC_k1, mDim2_K, mDim2_H, mC_k);
    } else {
      scalar_t value;
      value = kalmanParams.hits[i].normal[0] - kalmanParams.hits[i].ref[0];

      GSL_VECTOR_SET (mDim1_m, 0, value);
      GSL_MATRIX_SET (mDim1_V, 0, 0, kalmanParams.hits[i].err_locX * kalmanParams.hits[i].err_locX);

      correctGain  (mDim1_K, mC_k, mDim1_H, mDim1_V);
      correctState (mP_k1, mP_k, mDim1_K, mDim1_m, mDim1_H);
      correctCovariance (mC_k1, mDim1_K, mDim1_H, mC_k);
    }

    /* TODO: implemenent
    for (int c = 0; c < 5; ++c) {
      for (int k = 0; k < 5; ++k)
	kalmanParams.C_Values[c*5+k+i*25] = GSL_MATRIX_GET (mC_k1, c, k);
    }

    for (int c = 0; c < 5; ++c)
      kalmanParams.fitsForward[i].data[c] = GSL_VECTOR_GET (mP_k1, c);

#if DBG_LVL > 0
    std::cout << "P_k1: ";
    printVector (mP_k1);
#endif
*/

  }

  /* TODO: implement
  // Initial error covariance matrix estimate
  GSL_MATRIX_SET (mC_k1, 0, 0, 250.);
  GSL_MATRIX_SET (mC_k1, 1, 1, 250.);
  GSL_MATRIX_SET (mC_k1, 2, 2, 0.25);
  GSL_MATRIX_SET (mC_k1, 3, 3, 0.25);
  GSL_MATRIX_SET (mC_k1, 4, 4, 0.000001);

  for (int i = hitCount-1; i >= 0; --i) {
    GSL_MATRIX_VIEW jacobi;

    // Get the Transport-Jacobians
    jacobi = GSL_MATRIX_VIEW_ARRAY (kalmanParams.hits[i].jacobiInverse, 5, 5);

    predictState (mP_k, &jacobi.matrix, mP_k1);
    predictCovariance (mC_k, &jacobi.matrix, mC_k1);

#if HAS_MATERIAL_EFFECTS
    // add uncertainties from material effects:
    scalar_t covVal = GSL_MATRIX_GET(mC_k, 2, 2);
    GSL_MATRIX_SET(mC_k, 2, 2, covVal+kalmanParams.hits[i].sigmaDeltaPhiSquared);
    covVal = GSL_MATRIX_GET(mC_k, 3, 3);
    GSL_MATRIX_SET(mC_k, 3, 3, covVal+kalmanParams.hits[i].sigmaDeltaThetaSquared);
    covVal = GSL_MATRIX_GET(mC_k, 4, 4);
    GSL_MATRIX_SET(mC_k, 4, 4, covVal+kalmanParams.hits[i].sigmaDeltaQoverPSquared);
#endif

    if (kalmanParams.hits[i].is2Dim) { 
      scalar_t value1, value2;
      value1 = kalmanParams.hits[i].normal[0] - kalmanParams.hits[i].ref[0];
      value2 = kalmanParams.hits[i].normal[1] - kalmanParams.hits[i].ref[1];

      GSL_VECTOR_SET (mDim2_m, 0, value1);
      GSL_VECTOR_SET (mDim2_m, 1, value2);

      GSL_MATRIX_SET (mDim2_V, 0, 0, kalmanParams.hits[i].err_locX *
	  kalmanParams.hits[i].err_locX);
      GSL_MATRIX_SET (mDim2_V, 0, 1, kalmanParams.hits[i].cov_locXY);
      GSL_MATRIX_SET (mDim2_V, 1, 0, kalmanParams.hits[i].cov_locXY);
      GSL_MATRIX_SET (mDim2_V, 1, 1, kalmanParams.hits[i].err_locY *
	  kalmanParams.hits[i].err_locY);

      correctGain  (mDim2_K, mC_k, mDim2_H, mDim2_V);
      correctState (mP_k1, mP_k, mDim2_K, mDim2_m, mDim2_H);
      correctCovariance (mC_k1, mDim2_K, mDim2_H, mC_k);
    } else {
      scalar_t value;
      value = kalmanParams.hits[i].normal[0] - kalmanParams.hits[i].ref[0];

      GSL_VECTOR_SET (mDim1_m, 0, value);

      GSL_MATRIX_SET (mDim1_V, 0, 0, kalmanParams.hits[i].err_locX *
	  kalmanParams.hits[i].err_locX);

      correctGain  (mDim1_K, mC_k, mDim1_H, mDim1_V);
      correctState (mP_k1, mP_k, mDim1_K, mDim1_m, mDim1_H);
      correctCovariance (mC_k1, mDim1_K, mDim1_H, mC_k);
    }

    // FIXME: C_Inverse = C values im backward filter? (mit Rene klären)
    for (int c = 0; c < 5; ++c) {
      for (int k = 0; k < 5; ++k)
	kalmanParams.C_Inverse[c*5+k+i*25] = GSL_MATRIX_GET (mC_k1, c, k);
    }

    for (int c = 0; c < 5; ++c)
      kalmanParams.fitsBackward[i].data[c] = GSL_VECTOR_GET (mP_k1, c);

#if DBG_LVL > 0
    std::cout << "P_k1: ";
    printVector (mP_k1);
#endif

  }
  GSL_VECTOR *tmpVec;
  tmpVec = GSL_VECTOR_CALLOC(5);

  // Smoothing. FIXME: Save results to a useful struct..
  for (size_t i = 0; i < hitCount; ++i) {
    // Dim5_K = C_f * (C_f + C_b)^-1
    GSL_MATRIX_VIEW C_f, C_b;

    C_f = GSL_MATRIX_VIEW_ARRAY (&kalmanParams.C_Values[i*25], 5, 5);
    C_b = GSL_MATRIX_VIEW_ARRAY (&kalmanParams.C_Inverse[i*25], 5, 5);

    // mC_k = C_f + C_b
    GSL_MATRIX_MEMCPY (mC_k, &C_f.matrix);
    GSL_MATRIX_ADD (mC_k, &C_b.matrix);

    GSL_MATRIX_MEMCPY (mInverse, mC_k);

    // mInverse = (C_f + C_b)^-1
    matrixInverseSymmetric (mInverse);

    // mDim5_K = C_f * (C_f + C_b)^-1
    GSL_BLAS_XGEMM (CblasNoTrans, CblasNoTrans, 1.0f, &C_f.matrix, mInverse, 0.0f, mDim5_K);

    // P = p_f + K * p_b // SF: BUG!
    // P = p_f + K * (p_b - p_f)
    GSL_VECTOR_VIEW p_f, p_b;

    p_f = GSL_VECTOR_VIEW_ARRAY (kalmanParams.fitsForward[i].data, 5);
    p_b = GSL_VECTOR_VIEW_ARRAY (kalmanParams.fitsBackward[i].data, 5);
    //GSL_BLAS_XGEMV (CblasNoTrans, 1.0, mDim5_K, &p_b.vector, 0.0, mP_k);

    //tmpVec=p_b - p_f
    GSL_VECTOR_MEMCPY(tmpVec, &p_b.vector);
    GSL_VECTOR_SUB (tmpVec, &p_f.vector);

    // mP_k = K * (p_b - p_f)
    GSL_BLAS_XGEMV (CblasNoTrans, 1.0, mDim5_K, tmpVec, 0.0, mP_k);
    // mP_k = p_f + K * (p_b - p_f)
    GSL_VECTOR_ADD (mP_k, &p_f.vector);

#if BENCHMARK == 0
    x->addWert (GSL_VECTOR_GET (mP_k, 0));
    if (kalmanParams.hits[i].is2Dim)
      y->addWert (GSL_VECTOR_GET (mP_k, 1));
#endif

    // FIXME: calculation of smoothed covariance missing!
#if DBG_LVL > 2
    std::cout << "C_f: ";
    printMatrix (&C_f.matrix);
    std::cout << "C_b: ";
    printMatrix (&C_b.matrix);
    std::cout << "p_f: ";
    printVector (&p_f.vector);
    std::cout << "p_b: ";
    printVector (&p_b.vector);
    std::cout << "SmoothedP: p_f + (C_f . invert(C_f + C_b)) . p_b;";
    std::cout << "res: ";
    printVector (mP_k);
#endif

#if DBG_LVL > 0
    std::cout << "P_k: ";
    printVector (mP_k);
#endif
  }
*/

  free (kalmanParams.hits);
  free (kalmanParams.fitsForward);
  free (kalmanParams.fitsBackward);
  free (kalmanParams.C_Values);
  free (kalmanParams.C_Inverse);

  return 0;
}
