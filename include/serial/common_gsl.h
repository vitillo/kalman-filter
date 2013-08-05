/*
 * common_gsl.h
 *
 *  Created on: 31.05.2013
 *      Author: phil
 */

#ifndef COMMON_GSL_H_
#define COMMON_GSL_H_

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_vector_float.h>

#include <math.h>

#include "serial/common.h"

#if USE_DOUBLE

#define GSL_VECTOR                          gsl_vector
#define GSL_VECTOR_VIEW                     gsl_vector_view
#define GSL_MATRIX                          gsl_matrix
#define GSL_MATRIX_VIEW                     gsl_matrix_view

#define GSL_VECTOR_CALLOC(size1)            gsl_vector_calloc (size1)
#define GSL_VECTOR_FREE(v)                  gsl_vector_free (v)
#define GSL_VECTOR_GET(v, i)                gsl_vector_get (v, i)
#define GSL_VECTOR_SET(v, i, x)             gsl_vector_set (v, i, x)
#define GSL_VECTOR_MEMCPY(dest, src)        gsl_vector_memcpy (dest, src)
#define GSL_VECTOR_SUB(a, b)                gsl_vector_sub (a, b)
#define GSL_VECTOR_ADD(a, b)                gsl_vector_add (a, b)
#define GSL_VECTOR_VIEW_ARRAY(v, n)         gsl_vector_view_array (v, n)
#define GSL_VECTOR_SET_ZERO(v)              gsl_vector_set_zero (v)

#define GSL_MATRIX_CALLOC(size1, size2)     gsl_matrix_calloc (size1, size2)
#define GSL_MATRIX_SET(m, i, j, x)          gsl_matrix_set (m, i, j, x)
#define GSL_MATRIX_GET(m, i, j)             gsl_matrix_get (m, i, j)
#define GSL_MATRIX_FREE(m)                  gsl_matrix_free (m)
#define GSL_MATRIX_SET_IDENTITY(m)          gsl_matrix_set_identity (m)
#define GSL_MATRIX_SUB(a, b)                gsl_matrix_sub (a, b)
#define GSL_MATRIX_ADD(a, b)                gsl_matrix_add (a, b)
#define GSL_MATRIX_VIEW_ARRAY(b, n1, n2)    gsl_matrix_view_array (b, n1, n2)
#define GSL_MATRIX_MEMCPY(dest, src)		gsl_matrix_memcpy (dest, src);
#define GSL_MATRIX_SET_ZERO(m)              gsl_matrix_set_zero (m)

#define GSL_BLAS_XGEMM(TransA, TransB, alpha, A, B, beta, C) \
  gsl_blas_dgemm (TransA, TransB, alpha, A, B, beta, C)

#define GSL_BLAS_XGEMV(TransA, alpha, A, X, beta, Y) \
  gsl_blas_dgemv (TransA, alpha, A, X, beta, Y)

#else

#define GSL_VECTOR                          gsl_vector_float
#define GSL_VECTOR_VIEW                     gsl_vector_float_view
#define GSL_MATRIX                          gsl_matrix_float
#define GSL_MATRIX_VIEW                     gsl_matrix_float_view

#define GSL_VECTOR_CALLOC(size1)            gsl_vector_float_calloc (size1)
#define GSL_VECTOR_FREE(v)                  gsl_vector_float_free (v)
#define GSL_VECTOR_GET(v, i)                gsl_vector_float_get (v, i)
#define GSL_VECTOR_SET(v, i, x)             gsl_vector_float_set (v, i, x)
#define GSL_VECTOR_MEMCPY(dest, src)        gsl_vector_float_memcpy (dest, src)
#define GSL_VECTOR_SUB(a, b)                gsl_vector_float_sub (a, b)
#define GSL_VECTOR_ADD(a, b)                gsl_vector_float_add (a, b)
#define GSL_VECTOR_VIEW_ARRAY(v, n)         gsl_vector_float_view_array (v, n)
#define GSL_VECTOR_SET_ZERO(v)              gsl_vector_float_set_zero (v)

#define GSL_MATRIX_CALLOC(size1, size2)     gsl_matrix_float_calloc (size1, size2)
#define GSL_MATRIX_SET(m, i, j, x)          gsl_matrix_float_set (m, i, j, x)
#define GSL_MATRIX_GET(m, i, j)             gsl_matrix_float_get (m, i, j)
#define GSL_MATRIX_FREE(m)                  gsl_matrix_float_free (m)
#define GSL_MATRIX_SET_IDENTITY(m)          gsl_matrix_float_set_identity (m)
#define GSL_MATRIX_SUB(a, b)                gsl_matrix_float_sub (a, b)
#define GSL_MATRIX_ADD(a, b)                gsl_matrix_float_add (a, b)
#define GSL_MATRIX_VIEW_ARRAY(b, n1, n2)    gsl_matrix_float_view_array (b, n1, n2)
#define GSL_MATRIX_MEMCPY(dest, src)		gsl_matrix_float_memcpy (dest, src);
#define GSL_MATRIX_SET_ZERO(m)              gsl_matrix_float_set_zero (m)

#define GSL_BLAS_XGEMM(TransA, TransB, alpha, A, B, beta, C) \
  gsl_blas_sgemm (TransA, TransB, alpha, A, B, beta, C)

#define GSL_BLAS_XGEMV(TransA, alpha, A, X, beta, Y) \
  gsl_blas_sgemv (TransA, alpha, A, X, beta, Y)

#endif

inline void printMatrix (GSL_MATRIX *m) {
  std::cout << "matrix (\n";
  for (size_t i = 0; i < m->size1; ++i) {
    std::cout << "[";

    for (size_t j = 0; j < m->size2; ++j) {
      std::cout << GSL_MATRIX_GET (m, i, j);
      if (j < m->size2 - 1)
	std::cout << ", ";
    }

    if (i < m->size1 - 1)
      std::cout << "],\n";
    else
      std::cout << "]";
  }
  std::cout << ");";
  std::cout << std::endl;
}

inline void printVector (GSL_VECTOR *v) {
  std::cout << "[";
  for (size_t i = 0; i < v->size; ++i) {
    std::cout << GSL_VECTOR_GET (v, i);
    if (i < v->size - 1)
      std::cout << ", ";
  }
  std::cout << "];" << std::endl;
}

//Custom Error Handler so Program does not completely shutdown, instead print error and return nan-values for now
static void privGSLHandler (const char * reason, const char * file, int line, int gsl_errno) {
  printf ("error: %s. reason: %s. In file %s line %d\n", gsl_strerror (gsl_errno), reason, file, line);
}

// Calculate the pseudoinverse matrix of src
inline void
pseudoInverse (GSL_MATRIX *src)
{
  size_t m = src->size1, n = src->size2;

  gsl_matrix * U   = gsl_matrix_calloc (m, n);
  gsl_matrix * tmp = gsl_matrix_calloc (m, n);

  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j)
      gsl_matrix_set (U, i, j, (double) GSL_MATRIX_GET (src, i, j));
  }

  gsl_vector * Sdiag = gsl_vector_calloc (n);
  gsl_vector * work  = gsl_vector_calloc (n);
  gsl_matrix * V     = gsl_matrix_calloc (n, n);

  gsl_linalg_SV_decomp (U, V, Sdiag, work);
  gsl_vector_free (work);

  for (size_t i = 0; i < n; ++i)
    if (fabs (gsl_vector_get (Sdiag, i)) < 0.0000000001)
      gsl_vector_set (Sdiag, i, 0.0);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double Spinv_c_c = 0.0;

      if (gsl_vector_get (Sdiag, j) != 0.0)
	Spinv_c_c = 1.0 / gsl_vector_get (Sdiag, j);

      gsl_matrix_set (V, i, j, Spinv_c_c * gsl_matrix_get (V, i, j));
    }
  }

  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, V, U, 0.0, tmp);

  for (size_t i = 0; i < src->size1; ++i) {
    for (size_t j = 0; j < src->size2; ++j)
      GSL_MATRIX_SET (src, i, j, (scalar_t) gsl_matrix_get (tmp, i, j));
  }

  gsl_vector_free (Sdiag);
  gsl_matrix_free (U);
  gsl_matrix_free (V);
  gsl_matrix_free (tmp);

  return;
}

// Calculate the inverse of a matrix src.
inline void
matrixInverse (GSL_MATRIX *src)
{
  gsl_matrix * LU  = gsl_matrix_calloc (src->size1, src->size2);
  gsl_matrix * tmp = gsl_matrix_calloc (src->size1, src->size2);

  for (size_t i = 0; i < src->size1; ++i) {
    for (size_t j = 0; j < src->size2; ++j)
      gsl_matrix_set (LU, i, j, (double) GSL_MATRIX_GET (src, i, j));
  }

  int s;
  gsl_permutation * perm = gsl_permutation_calloc (LU->size1);
  gsl_linalg_LU_decomp (LU, perm, &s);

  // LU decomp singular?
  for (size_t i = 0; i < LU->size1; i++) {
    if (gsl_matrix_get (LU, i, i) == 0) {

#if DBG_LVL > 2
      std::cout << "Matrix is singular! Using pseudoinverse matrix." << std::endl;
#endif
      gsl_matrix_free      (LU);
      gsl_permutation_free (perm);
      gsl_matrix_free      (tmp);

      pseudoInverse (src);
      return;
    }
  }

  gsl_linalg_LU_invert (LU, perm, tmp);

  for (size_t i = 0; i < src->size1; ++i) {
    for (size_t j = 0; j < src->size2; ++j)
      GSL_MATRIX_SET (src, i, j, (scalar_t) gsl_matrix_get (tmp, i, j));
  }

  gsl_matrix_free (LU);
  gsl_matrix_free (tmp);
  gsl_permutation_free (perm);

  return;
}

inline void 
matrixInverseSymmetric(GSL_MATRIX *src)
{
  // Speedup matrix inversion for our specific case...
  if(src->size1 == 1){
    GSL_MATRIX_SET(src, 0, 0, 1/GSL_MATRIX_GET(src, 0, 0));

    return;
  }else if(src->size1 == 2){
    scalar_t a = GSL_MATRIX_GET(src, 0, 0); 
    scalar_t b = GSL_MATRIX_GET(src, 0, 1); 
    scalar_t c = GSL_MATRIX_GET(src, 1, 0); 
    scalar_t d = GSL_MATRIX_GET(src, 1, 1); 

    scalar_t inv_det = 1/(a*d - b*c);
    GSL_MATRIX_SET(src, 0, 0, d*inv_det);
    GSL_MATRIX_SET(src, 0, 1, -b*inv_det);
    GSL_MATRIX_SET(src, 1, 0, -c*inv_det);
    GSL_MATRIX_SET(src, 1, 1, a*inv_det);

    return;
  }

  gsl_matrix *tmp = gsl_matrix_calloc (src->size1, src->size2);
  gsl_matrix *dst = gsl_matrix_calloc (src->size1, src->size2);
  gsl_permutation *perm = gsl_permutation_calloc(src->size1);
  int s;

  for (size_t i = 0; i < src->size1; ++i) {
    for (size_t j = 0; j < src->size2; ++j)
      gsl_matrix_set (tmp, i, j, (double) GSL_MATRIX_GET (src, i, j));
  }

  gsl_linalg_LU_decomp(tmp, perm, &s);
  gsl_linalg_LU_invert(tmp, perm, dst);

  for (size_t i = 0; i < tmp->size1; ++i) {
    for (size_t j = 0; j < tmp->size2; ++j)
      GSL_MATRIX_SET (src, i, j, (scalar_t) gsl_matrix_get (dst, i, j));
  }

  gsl_matrix_free (tmp); 
  gsl_matrix_free (dst); 
  gsl_permutation_free(perm);

  return;

  /*
  // Cholesky decomp only works with gsl_matrix (no float support)
  gsl_matrix *tmp = gsl_matrix_calloc (src->size1, src->size2);

  for (size_t i = 0; i < src->size1; ++i) {
  for (size_t j = 0; j < src->size2; ++j)
  gsl_matrix_set (tmp, i, j, (double) GSL_MATRIX_GET (src, i, j));
  }

  gsl_linalg_cholesky_decomp (tmp);
  gsl_linalg_cholesky_invert (tmp);

  for (size_t i = 0; i < tmp->size1; ++i) {
  for (size_t j = 0; j < tmp->size2; ++j)
  GSL_MATRIX_SET (src, i, j, (scalar_t) gsl_matrix_get (tmp, i, j));
  }

  gsl_matrix_free (tmp); 

  return;*/
}

inline void
matrixInverse (scalar_t* matrix, scalar_t* mInverse, unsigned int columns, unsigned int rows)
{
  gsl_error_handler_t *old_handler = gsl_set_error_handler (&privGSLHandler);

  GSL_MATRIX * gA = GSL_MATRIX_CALLOC (rows, columns);

  for (unsigned int i = 0; i < rows; i++)
    for (unsigned int j = 0; j < columns; j++)
      GSL_MATRIX_SET (gA, i, j, matrix[i * rows + j]);

  matrixInverse (gA);

  for (unsigned int i = 0; i < columns; i++)
    for (unsigned int j = 0; j < rows; j++)
      mInverse[i * rows + j] = GSL_MATRIX_GET (gA, i, j);

  gsl_set_error_handler (old_handler);

  GSL_MATRIX_FREE (gA);
}

inline void
matrixInverse (scalar_t* matrixStart, scalar_t* matrixInverseStart, unsigned int columns,
    unsigned int rows, unsigned int batch)
{
  gsl_error_handler_t *old_handler = gsl_set_error_handler (&privGSLHandler);

  GSL_MATRIX * gA   = GSL_MATRIX_CALLOC (rows, columns);

  for (unsigned int k = 0; k < batch; k++) {
    scalar_t *matrix = matrixStart + (rows*columns*k);
    scalar_t *mInverse = matrixInverseStart + (rows*columns*k);

    for (unsigned int i = 0; i < rows; i++)
      for (unsigned int j = 0; j < columns; j++)
	GSL_MATRIX_SET (gA, i, j, matrix[i * rows + j]);

    matrixInverse (gA);

    for (unsigned int i = 0; i < columns; i++)
      for (unsigned int j = 0; j < rows; j++)
	mInverse[i * rows + j] = GSL_MATRIX_GET (gA, i, j);
  }

  gsl_set_error_handler (old_handler);

  GSL_MATRIX_FREE (gA);
}

inline void
matrixInverseSymmetric (scalar_t* matrixStart, scalar_t* matrixInverseStart, unsigned int columns,
    unsigned int rows, unsigned int batch)
{
  gsl_error_handler_t *old_handler = gsl_set_error_handler (&privGSLHandler);

  GSL_MATRIX * gA   = GSL_MATRIX_CALLOC (rows, columns);

  for (unsigned int k = 0; k < batch; k++) {
    scalar_t *matrix = matrixStart + (rows*columns*k);
    scalar_t *mInverse = matrixInverseStart + (rows*columns*k);

    for (unsigned int i = 0; i < rows; i++)
      for (unsigned int j = 0; j < columns; j++)
	GSL_MATRIX_SET (gA, i, j, matrix[i * rows + j]);

    matrixInverseSymmetric (gA);
    for (unsigned int i = 0; i < columns; i++)
      for (unsigned int j = 0; j < rows; j++)
	mInverse[i * rows + j] = (gpu_scalar_t)GSL_MATRIX_GET (gA, i, j);
  }

  gsl_set_error_handler(old_handler);
  GSL_MATRIX_FREE (gA);
}

#endif /* COMMON_GSL_H_ */
