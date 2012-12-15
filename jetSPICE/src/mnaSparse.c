/*
 * mnaSparse.c
 *
 *  Created on: Jan 18, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */


#include "mnaSparse.h"
#include <stdio.h>

//#include "mna.h"
//#include "csparse.h"


extern int vnumber;
extern int lnumber;
extern int rnumber;
extern int inumber;

extern commandList *cmdList;

void addElementStampSparse( zeroCircuit *node, mnaSpSystem *system, cs *arrayA, hashTable *namTab, int row, int kPos)
{
	int posn1;
	int posn2;
	double tmp;

	struct linear_element *element;

	element = node->linElement;


	switch ( node->type )
	{
	case 0:
		switch( node->linElement->element )
		{
		case 1:

			if ( ( strcmp( node->linElement->connectors[0], "0\0" ) ) && ( strcmp( node->linElement->connectors[1], "0\0" ) ) )
			{
				//posn1 = findPlace ( node->linElement->connectors[0], &nodes );
				//posn2 = findPlace ( node->linElement->connectors[1], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );

				tmp = (-1.0) * ( 1.0 / node->linElement->value);
				cs_entry(arrayA,posn1,posn2,tmp);
				cs_entry(arrayA,posn2,posn1,tmp);

			}

			if ( strcmp( node->linElement->connectors[0], "0\0" ) )
			{
				//posn1 = findPlace ( node->linElement->connectors[0], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				tmp = ( 1.0 / node->linElement->value );
				cs_entry(arrayA,posn1,posn1,tmp);
			}

			if ( strcmp( node->linElement->connectors[1], "0\0" ) )
			{
				//posn2 = findPlace ( node->linElement->connectors[1], &nodes );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );
				tmp = ( 1.0 / node->linElement->value );
				cs_entry(arrayA,posn2,posn2,tmp);
			}

			break;
		case 3:

			if ( strcmp( node->linElement->connectors[0], "0\0" ) )
			{
				//posn1 = findPlace ( node->linElement->connectors[0], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				system->vector_B[posn1] += (-1.0) * element->value;
			}

			if ( strcmp( node->linElement->connectors[1], "0\0" ) )
			{
				//posn2 = findPlace ( node->linElement->connectors[1], &nodes );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );
				system->vector_B[posn2] += element->value;
			}

			break;
		case 2: //case 2same as 5
		case 5:
			element->k = kPos;
			if (element->element == 2)
			{
				system->vector_B[ namTab->currPos -1 + kPos ] += element->value;

				//printf ( "bug ib B: %f\n", element->value );
			}


			if ( ( strcmp( node->linElement->connectors[0], "0\0" ) ) && ( strcmp( node->linElement->connectors[1], "0\0" ) ) )
			{
				//posn1 = findPlace ( element->connectors[0], &nodes );
				//posn2 = findPlace ( element->connectors[1], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );

				cs_entry(arrayA,kPos + namTab->currPos -1,posn1,1.0);
				cs_entry(arrayA,kPos + namTab->currPos -1,posn2,-1.0);

				cs_entry(arrayA,posn1,kPos + namTab->currPos -1,1.0);
				cs_entry(arrayA,posn2,kPos + namTab->currPos -1,-1.0);
				break;

			}

			if ( strcmp( node->linElement->connectors[0], "0\0" ) )
			{
				//posn2 = findPlace ( element->connectors[1], &nodes );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );

				cs_entry(arrayA,kPos + namTab->currPos -1,posn2,-1.0);
				cs_entry(arrayA,posn2,kPos + namTab->currPos -1,-1.0);

				break;
			}

			if ( strcmp( node->linElement->connectors[1], "0\0" ) )
			{
				//posn1 = findPlace ( element->connectors[0], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );

				cs_entry(arrayA,kPos + namTab->currPos -1,posn1,1.0);
				cs_entry(arrayA,posn1,kPos + namTab->currPos -1,1.0);

			}

			break;
		}

		break;
	case 1:
		//irrelevant
		break;
	}

}

void buildCSp ( zeroCircuit *node, cs *arrayC, hashTable *namTab, int kPos )
{
	int posn1;
	int posn2;
	double tmp;

	struct linear_element *element;

	element = node->linElement;


	switch ( node->type )
	{
	case 0:
		switch( node->linElement->element )
		{
		case 4:
			if ( ( strcmp( node->linElement->connectors[0], "0\0" ) ) && ( strcmp( node->linElement->connectors[1], "0\0" ) ) )
			{
				//posn1 = findPlace ( node->linElement->connectors[0], &nodes );
				//posn2 = findPlace ( node->linElement->connectors[1], &nodes );

				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );

				tmp = (-1.0) * node->linElement->value;
				cs_entry(arrayC,posn1,posn2,tmp);
				cs_entry(arrayC,posn2,posn1,tmp);

			}

			if ( strcmp( node->linElement->connectors[0], "0\0" ) )
			{
				//posn1 = findPlace ( node->linElement->connectors[0], &nodes );
				posn1 = findHash ( namTab, node->linElement->connectors[0] );
				tmp = node->linElement->value;
				cs_entry(arrayC,posn1,posn1,tmp);
			}

			if ( strcmp( node->linElement->connectors[1], "0\0" ) )
			{
				//posn2 = findPlace ( node->linElement->connectors[1], &nodes );
				posn2 = findHash ( namTab, node->linElement->connectors[1] );
				tmp = node->linElement->value;
				cs_entry(arrayC,posn2,posn2,tmp);
			}

			break;

		case 5:

			tmp = -node->linElement->value;
			cs_entry(arrayC, kPos + namTab->currPos - 1, kPos + namTab->currPos - 1, tmp);
			break;
		}
		break;

	case 1:
		//irrelevant
		break;
	}

}

mnaSpSystem *formatMnaTableSparse ( zeroCircuit *circuit, hashTable *names, int v_num, int l_num, int r_num, int i_num )
{
    int rows;
    int k_pos;
    int maxEstimated;

    zeroCircuit *cur;
    mnaSpSystem *finalSystem;
    cs *A, *Cm;

    //initialize k position
    k_pos = 0;

    finalSystem = ( mnaSpSystem * )malloc ( sizeof( mnaSpSystem ) );

    rows = ( names->currPos ) + v_num + l_num ;		//dimension of original matrix
    maxEstimated = 2 * ( 4 * r_num + 1 * v_num + 2 * i_num );

    A = cs_spalloc(rows, rows, maxEstimated, 1, 1);

    if ( A == NULL )
    {
    	printf( "!!!WARNING!!! Table A was NOT allocated successfully\n" );
    	printf( "\t-Dimension of system to be allocated : %d\n\n", rows );
    }
    //Transient Analysis : Build C
    if (cmdList->_transient == 1)
    {
    	Cm = cs_spalloc(rows,rows,maxEstimated,1,1);
    }

    finalSystem->vector_B = (double *)calloc(rows, sizeof(double));
    finalSystem->vector_X = (char **)malloc( rows * sizeof (char *) );

    finalSystem->dim = rows;

    //for ( cur = circuit->next; (cur->linElement == NULL && cur->nonlinElement == NULL); cur = cur->next )
    for( cur = circuit->next; ( (cur->linElement != NULL) || (cur->nonlinElement != NULL) ); cur = cur->next )
    {
    	if (( cur->linElement->element == 2 ) || ( cur->linElement->element == 5 ))  //if element is voltage source or inductor
    	{
    		k_pos++;
    	}

    	addElementStampSparse ( cur, finalSystem, A, names, rows, k_pos );

    	//Transient Analysis : Build C
    	if (cmdList->_transient == 1)
    	{
    		buildCSp( cur, Cm, names, k_pos);
    	}

    }

    finalSystem->array_A = cs_compress(A);
    cs_dupl(finalSystem->array_A);
    cs_spfree(A);

    //Transient Analysis : Build C
    if (cmdList->_transient == 1)
    {
    	finalSystem->array_C = cs_compress(Cm);
    	cs_dupl(finalSystem->array_C);
    	cs_spfree(Cm);
    }


    return ( finalSystem );
}

void printMNASparse (hashTable *names, mnaSpSystem *systemMna)
{
	int i;
	char *str="matrixG.txt";
	char *str1="matrixC.txt";

	FILE *bfile;

	printf( "\tINFO : Group 2 elements:\n\t\t-Voltage Sources: %d\n\t\t-Inductors: %d\n", vnumber, lnumber );
	printf( "\tINFO : Group 1 elements:\n\t\t-Resistors: %d\n\t\t-Current Sources: %d\n", rnumber, inumber );
	printf( "\tINFO : Number of nodes of sub-Matrix G: %d\n", names->currPos );

	//Sparse
	cs_print(systemMna->array_A,str,0);
	cs_print(systemMna->array_C,str1,0);

	bfile = fopen ( "b_vector.txt", "w" );

	fprintf( bfile, "\tINFO : Group 2 elements:\n\t\t-Voltage Sources: %d\n\t\t-Inductors: %d\n", vnumber, lnumber );
	fprintf( bfile, "\tINFO : Group 1 elements:\n\t\t-Resistors: %d\n\t\t-Current Sources: %d\n", rnumber, inumber );
	fprintf( bfile, "\tINFO : Number of nodes of sub-Matrix G: %d\n", names->currPos );
	fprintf( bfile, "\n\n\n\n" );

	fprintf ( bfile, "\t   -----B Vector----- \n\n" );
	for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
	{
		fprintf ( bfile, "pos(%d)\t\t\t%.10f\n", i, systemMna->vector_B[i] );
	}

	fclose ( bfile );

}
