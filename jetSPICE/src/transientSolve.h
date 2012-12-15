/*
 * transientSolve.h
 *
 *  Created on: Feb 19, 2012
 *      Author: gepinita
 */

#ifndef TRANSIENTSOLVE_H_
#define TRANSIENTSOLVE_H_


#include <stdio.h>
#include <stdlib.h>
#include "csparse.h"
#include "linearSolve.h"
#include "mnaSparse.h"
#include "math.h"

void backwardEuler(mnaSystem *system, zeroCircuit *circuit, hashTable *hashTab);
void trapezoidal(mnaSystem *system, zeroCircuit *circuit, hashTable *hashTab);
void cs_backwardEuler(mnaSpSystem *system, zeroCircuit *circuit, hashTable *hashTab);
void cs_trapezoidal(mnaSpSystem *system, zeroCircuit *circuit, hashTable *hashTab);
void constructE(double *e, double t, int *s, int *periods, hashTable *hashTab, zeroCircuit *circuit);
void formatRowVecM ( double **rV, double *A, int size );

#endif /* TRANSIENTSOLVE_H_ */
