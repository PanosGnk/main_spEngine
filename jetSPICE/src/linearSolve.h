/*
 * linearSolve.h
 *
 *  Created on: Nov 4, 2011
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 *      Version: 0.1
 *  Description: Implementation of Direct Methods (LU Decomposition, Cholesky Factorization, upper and lower triangular solvers) and
 *  			 Iterative Methods (Conjugate Gradient, Bi-Conjugate Gradient and GMRES).
 */

#ifndef LINEARSOLVE_H_
#define LINEARSOLVE_H_

#include <stdlib.h>
#include <stdio.h>
#include "csparse.h"

#define EPS 1.0e-14

/********************************************************************************
 *                                                                              *
 *                            FUNCTION DECLARATIONS                             *
 *                                                                              *
 ********************************************************************************/

/**
 * Function that computes Cholesky factorization of a matrix.
 * @param rowVector rectangular matrix to factorize and factorized matrix on output.
 * @param size size of matrix.
 * @param diag vector which contains the diagonal elements of factorized matrix.
 */
void choleskyDecomp ( double **rowVector, int size, double *diag );


/**
 * Function that solves x=A\b. using Cholesky factorization.
 * @param rowVector rectangular factorized matrix.
 * @param size size of matrix.
 * @param diag vector which contains the diagonal elements of factorized matrix.
 * @param bVector right-hand vector.
 * @return solution vector.
 */
double *choleskySolve ( double **rowVector, int size, double *diag, double *bVector );


/**
 * Function for performing LU decomposition with partial (row) pivoting.
 * @param rowVector rectangular factorized matrix.
 * @param size size of matrix.
 * @param pivotTable pivot vector.
 * @result a matrix with L in the lower triangular and U in the upper triangular.
 */
double **luDecomp ( double **rowVector, int size, int *pivotTable );


/**
 * Function for performing Backward Substitution.
 * @param rowVector LU matrix.
 * @param size size of matrix and vectors.
 * @param pivotTable pivot vector.
 * @param bVector right-hand vector.
 * @result solution vector.
 */
double *luBackSubst ( double **rowVector, int size, int *pivotTable, double *bVector );


/**
 * Function that computes Jacobi Preconditioner, M = diag(A).
 * @param M vector that contains the diagonal elements of A on output.
 * @param A input matrix.
 * @param size size of rectangular matrix.
 */
void jacobiPrecond (double *M, double **A, int size);


/**
 * Function that performs Preconditioner Solve, M * z^[i-1]=r^[i-1].
 * @param M Jacobi preconditioner.
 * @param r residual vector.
 * @param z result vector.
 * @param size size of vectors.
 */
void precSolve (double *M, double *r, double *z, int size);


/**
 * Conjugate Gradient Method
 * @param rowVector input rectangular matrix.
 * @param bVector right-hand vector.
 * @param size size of matrix and vectors.
 * @param iter iteration number on output.
 * @param itol tolerance threshold.
 * @result solution vector.
 */
double *cg (double **rowVector, double *bVector, int size, int *iter, double itol );


/**
 * Function that implements Bi-Conjugate Gradient Method.
 * @param rowVector input rectangular matrix.
 * @param bVector right-hand vector.
 * @param size size of matrix and vectors.
 * @param iter iteration number on output.
 * @param itol tolerance threshold.
 * @return solution vector.
 */
double *bicg (double **rowVector, double *bVector, int size, int *iter, double itol );

/**
 * Update Function for GMRES(m) Method.
 */
void update(double *x, double *u, double **A, double *b,double *h, double *s, int resd,int hsize, int size,int iter);


/**
 * Function that implements GMRES Method.
 * @param rowVector input rectangular matrix.
 * @param bVector right-hand vector.
 * @param size size of matrix and vectors.
 * @param iter iteration number on output.
 * @param res
 * @param itol tolerance threshold.
 * @return solution vector.
 */
double *gmres(double **rowVector, double *bVector, int size, int *iter,int res, double itol);


/**
 * Function that computes Jacobi Preconditioner of a Sparse matrix, M = diag(A).
 * @param M vector that contains the diagonal elements of A on output.
 * @param A sparse matrix.
 * @param size size of rectangular matrix.
 */
void cs_jacobiPrecond (double *M, cs *A, int size);


/**
 * Function that implements Conjugate Gradient Method using Sparse Matrices.
 * @param A sparse matrix.
 * @param bVector right-hand vector.
 * @param size size of matrix and vectors.
 * @param iter iteration number on output.
 * @param itol tolerance threshold.
 * @return solution vector.
 */
double *cs_cg (cs *A, double *bVector, int size, int *iter, double itol );


/**
 * Function that implements Bi-Conjugate Gradient Method using Sparse Matrices.
 * @param A sparse matrix.
 * @param bVector right-hand vector.
 * @param size size of matrix and vectors.
 * @param iter iteration number on output.
 * @param itol tolerance threshold.
 * @return solution vector.
 */
double *cs_bicg (cs *A, double *bVector, int size, int *iter, double itol );


/**
 * tester function-no actual use
 */
void solveSystem ();


/********************************************************************************
 *                                                                              *
 *                            UTILITY FUNCTIONS                                 *
 *                                                                              *
 ********************************************************************************/

/**
 * Function that pivots vectorX (string vector).
 * @param pivots pivot vector.
 * @param vec string vector.
 * @param size size of n both vectors.
 * @return transposed vector.
 */
char **transposeVectorX ( int *pivots, char **vec, int size);

/**
 *
 */
int findVarPos ( char *str, char **vector , int size);


/**
 * Function that clears a vector, x[i]=0.
 * @param x input vector and solution on output.
 * @param size size of vector.
 */
void clearVec (double *x, int size);


/**
 * Function that computes the Inner Product of two vectors.
 * @param x input vector.
 * @param x input vector.
 * @param size size of vectors.
 * @return the inner product of the two vectors.
 */
double innerProd(double *x, double *y, int size);


/**
 * Function that computes the Norm-2 of a vector.
 * @param x input vector.
 * @param size size of vector.
 * @return the norm-2 of the vector x.
 */
double norm2(double *x, int size);

/**
 * Function that performs a Vector Copy, x[i] = y[i].
 * @param x solution on output.
 * @param y input vector.
 * @param size size of vectors.
 */
void vecCpy (double *x, double *y, int size);

/**
 * Function that performs a Scalar Vector Addition, x[i] += s * y[i].
 * @param x solution on output.
 * @param y input vector.
 * @param s multiplication factor for vector y.
 * @param size size of vectors.
 */
void vecCpyScal(double *x, double *y, double s, int size);


/**
 * Function that performs a Scalar Vector Sum, x[i] = y[i] + s * z[i].
 * @param x solution on output.
 * @param y input vector.
 * @param z input vector.
 * @param s multiplication factor for vector z.
 * @param size size of vectors.
 */
void vecSumScal(double *x, double *y, double *z, double s, int size);


/**
 * Function that performs a Scalar Vector Division, x[i] = y[i] / s.
 * @param x solution on output.
 * @param y input vector.
 * @param s division factor for vector y.
 * @param size size of vector.
 */
void vecDivScal(double *x, double *y, double s, int size);


/**
 * Function that performs a Vector Difference, x[i] = y[i] - z[i].
 * @param x solution on output.
 * @param y input vector.
 * @param z input vector.
 * @param size size of vectors.
 */
void vecDiff(double *x, double *y, double *z, int size);


/**
 * Function that performs a Scalar Vector Difference, x[i] = y[i] - ( s * z[i] ).
 * @param x solution on output.
 * @param y input vector.
 * @param z input vector.
 * @param s multiplication factor for vector z.
 * @param size size of vectors.
 */
void vecDiffScal(double *x, double *y, double *z, double s, int size);


/**
 * Funtion that performs a Vector-Number Multiplication x[i]= s * y[i].
 * @param x solution on output.
 * @param y input vector.
 * @param s multiplication factor for vector y.
 * @param size size of vector.
 */
void vecMult(double *x, double *y, double s, int size);

/**
 * Function that inverses a vector, .
 * @param x input vector and solution on output.
 * @param size size of vector.
 */
void invVec (double *x, int size);


/**
 * Function that implements a rectangular Matrix-Vector Product of size n, Ax=A*x.
 *
 * @param Ax solution on output.
 * @param A multiplicant Matrix.
 * @param x Multiplier Vector.
 * @param size size n of both vector and matrix.
 * @param method regular multiplication (method=0) or transpose matrix multiplication (method=1).
 */
void mvProd (double *Ax, double **A, double *x, int size, int method);


/**
 * Function that prints vector.
 *
 * @param x vector.
 * @param size size of vector.
 */
void printVec(double *x, int size);


/**
 * Function that add two Matrices.
 * @param M result matrix.
 * @param alpha multiplicant of A.
 * @param A matrix A
 * @param beta multiplicant of A
 * @param B matirx B.
 * @param size size of the system.
 *
 */
void addMMScal(double *M, double alpha, double *A, double beta, double *B, int size);

#endif /* LINEARSOLVE_H_ */
