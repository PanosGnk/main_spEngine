/*
 * mnaSparse.h
 *
 *  Created on: Oct 23, 2011
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 *      Version	: 0.1.1
 *  Description : Functions to form the system of Modified Nodal Analysis method
 */

#ifndef MNASPARSE_H_
#define MNASPARSE_H_



#include "mna.h"
#include "csparse.h"



struct MNA_SpSystem
{
	double 	*vector_B;
	char 	**vector_X;
	cs *array_A;
	cs *array_C;

	int 	dim;
	double 	**rowVector;
	double 	*dcPointRes;
	double 	*trRes;
};

typedef struct MNA_SpSystem mnaSpSystem;

/**
 * Function that leaves an element stamp on G matrix.
 * @param node the element from the list.
 * @param systemMna final Mna system.
 * @param arrayA the G matrix.
 * @param nodes table with the nodes of the system.
 * @param row row number of the G table.
 * @param kPos counter that extends the C matrix
 */
void addElementStampSparse( zeroCircuit *node, mnaSpSystem *system, cs *arrayA, hashTable *namTab, int row, int kPos);

/**
 * Function that leaves an element stamp on C matrix.
 * @param node the element from the list.
 * @param arrayA the C matrix.
 * @param nodes table with the nodes of the system.
 * @param kPos counter that extends the C matrix
 */
void buildCSp ( zeroCircuit *node, cs *arrayC, hashTable *namTab, int kPos );

/**
 * Function that builds the MNA matrices (G in DC, both G and C in transient analysis).
 * @param circuit list with the circuit elements
 * @param table table with the nodes of the system.
 * @param v_num number of voltage sources.
 * @param l_num number of inductors.
 * @return the MNA system.
 */
mnaSpSystem *formatMnaTableSparse ( zeroCircuit *circuit, hashTable *names, int v_num, int l_num ,int r_num, int i_num );

/**
 * Function that prints the MNA System.
 * @param tableN table with the nodes of the system.
 * @param systemMna final MNA system.
 */
void printMNASparse (hashTable *names, mnaSpSystem *systemMna);

#endif /* MNASPARSE_H_ */
