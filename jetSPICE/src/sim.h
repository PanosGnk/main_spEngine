/*
 * sim.h
 *
 *  Created on: Nov 7, 2011
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 * Description:	Configure and run simulation routines.
 *
 */

#ifndef SIM_H_
#define SIM_H_

//commandList *cmdList;
#include "types.h"
#include "mnaSparse.h"
#include "mna.h"


void initializeCommandList ();
void setMemoryPolicy ( hashTable *namTab, int lnum, int vnum );
void configOptions ( int *options );
void configDcSweep ( double *values, char **srcs, int num );
void configPlot( char **plots, int num );
void configItol ( double tolerance );
int configMethod ( char *str );
traceFileList *initializeTraceFiles ();
void initializeFiles ( traceFileList *flist );
int *findVarPositions ( char **vec, int size );
int *findPositions ( hashTable *namTab );
int checkDCPos (  zeroCircuit *circuit );
int checkDCPosSp (  zeroCircuit *circuit );

void dcPointSim ( mnaSystem *sys );
void dcSweepSim ( mnaSystem *sys, traceFileList *flist, hashTable *names, zeroCircuit *cir );
void dcPointSimSparse ( mnaSpSystem *sys );
void dcSweepSimSparce ( mnaSystem *sys, traceFileList *flist, hashTable *names, zeroCircuit *cir );

#endif /* SIM_H_ */
