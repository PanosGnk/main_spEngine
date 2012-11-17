/*
 ============================================================================
 Name        : mainFlow.c
 Author      : Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 Version     : 0.1
 Copyright   : engineeR
 Description : Main flow of the SPICE Simulator
 ============================================================================
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/** User defined libraries **/
#include "types.h"
#include "csparse.h"
#include "sim.h"
#include "linearSolve.h"
#include "transientParse.h"
#include "spiceParse.h"
#include "hashtable.h"
//#include "transientSolve.h"







int main( int argc, char *argv[]) {

	mnaSystem *_systemMna;
	mnaSpSystem *_systemMnaSp;

	zeroCircuit *_mainCircuit;
	hashTable 	*_nameTab;

	traceFileList *_trfiles;


	/*******************************
	 *   Various Initializations   *
	 *******************************/


	vnumber = 0;
	lnumber = 0;
	rnumber = 0;
	inumber = 0;

	initializeCommandList ();

	_nameTab = createHash ( SIZE );


	/*******************************
	 *      Parse Spice File       *
	 *******************************/

	if ( argc == 1 )
	{
		printf ( "ERROR : No filename specified\n" );
		exit ( 0 );
	}
	printf ( "- File Registered : %s\n", argv[1] );
	_mainCircuit = parseNetlist ( argv[1], _nameTab );
	printf ( "- Completed Parsing File : %s\n", argv[1] );


	/*******************************
	 *      Netlist Analysis       *
	 *******************************/

	//estimate memory allocation policy
	if ( argc > 2 )
	{
		int iarg;
		int found = 0;

		for ( iarg = 0; iarg < argc; iarg ++ )
		{
			if ( !strcmp ( argv[iarg], "-sparse" ) ||  !strcmp ( argv[iarg], "-SPARSE" ) )
				found = 1;
		}
		if ( found )
			cmdList->_useSparse = 1;
	}
	else
		setMemoryPolicy ( _nameTab, lnumber, vnumber );


	if ( cmdList->_useSparse )
	{
		printf( "\tINFO : SPARSE Matrix allocation will be used\n" );
		printf("- Starting Sparse Matrix Construction\n");
		fflush(stdout);

		_systemMnaSp = formatMnaTableSparse( _mainCircuit, _nameTab, vnumber, lnumber, rnumber, inumber );
		printMNASparse ( _nameTab, _systemMnaSp);

		printf("- Sparse Matrix Construction Complete\n");
		fflush(stdout);
	}
	else
	{
		printf( "\tINFO : Dense Matrix allocation will be used\n" );
		printf("- Starting Dense Matrix Construction\n");
		fflush(stdout);

		_systemMna = formatMnaTable ( _mainCircuit, _nameTab, vnumber, lnumber );
		printMNA ( _nameTab, _systemMna);

		printf("- Dense Matrix Construction Complete\n");
		fflush(stdout);
	}


	/*******************************
	 *          Simulation         *
	 *******************************/

	printf("- Starting Fixed Point DC Analysis\n");
	fflush(stdout);

	if ( cmdList->_useSparse )
	{
		dcPointSimSparse ( _systemMnaSp );
	}
	else
	{
		dcPointSim ( _systemMna );
	}

	printf("- Fixed Point DC Analysis Complete\n");
	fflush(stdout);


	int i;

		for ( i=0; i < cmdList->_plotnum; i++)
		{
			printf ( "\t\t-node : %s\n", cmdList->_plot[i]);
		}



	if ( cmdList->_plotnum )
	{
		printf ( "- Starting DC Sweep Analysis\n" );
		fflush(stdout);

		_trfiles = initializeTraceFiles();
		initializeFiles ( _trfiles );

		if ( cmdList->_useSparse )
		{
			dcSweepSimSparse ( _systemMnaSp, _trfiles, _nameTab );
		}
		else
		{
			dcSweepSim ( _systemMna, _trfiles, _nameTab );
		}

		printf ( "- DC Sweep Analysis Complete\n" );
		fflush(stdout);
	}


	return EXIT_SUCCESS;
}
