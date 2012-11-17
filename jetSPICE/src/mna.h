/*
 * mna.h
 *
 *  Created on: Oct 23, 2011
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 *      Version	: 0.1.1
 *  Description : Functions to form the system of Modified Nodal Analysis method
 */

#ifndef MNA_H_
#define MNA_H_

#include "types.h"
#include "hashtable.h"



struct nodeTableZero
{
	int *table;
	int size;

};

typedef struct nodeTableZero nodeTable;

struct MNA_system
{
	double 	*array_A;
	double 	*vector_B;
	char 	**vector_X;

	int 	dim;
	double 	**rowVector;
};

typedef struct MNA_system mnaSystem;



nodeTable initializeTable ();
//int findPlace(int connector,nodeTable *table );
nodeTable  analyseNodes( zeroCircuit *node, nodeTable nodes );
void printNodes( nodeTable nodes );
void initVectorX ( zeroCircuit *node, mnaSystem *system, hashTable *table  );
void fillVectorX ( zeroCircuit *node, mnaSystem *system, int k, int origsize );
void addElementStamp( zeroCircuit *node, mnaSystem *system, hashTable *namTab, int row, int kPos );
void formatRowVec ( mnaSystem *finalSystem );
mnaSystem *formatMnaTable ( zeroCircuit *circuit, hashTable *nameTab, int v_num, int l_num );
void printMNA (hashTable *names, mnaSystem *systemMna);

#endif /* MNA_H_ */
