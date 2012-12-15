/*
 ============================================================================
 Name        : spiceParce.h
 Author      : Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 Version     : 0.1
 Copyright   : engineeR
 Description : Parse functions for spice format file parsing
 ============================================================================
 */





#include "hashtable.h"
#include <string.h>
//#include "types.h"

#define BUFSIZE 512



extern commandList *cmdList;

int optionState[2] = { 0, 0 };  //pos0:SPD, pos1:ITER

extern int vnumber;
extern int lnumber;
extern int rnumber;
extern int inumber;

/**
 * detectEmptyComment
 *
 * Detects if the string str is empty or is a SPICE syntax comment
 *
 * @param str : input string
 * @return 0:is empty/is comment 1:non-empty line
*/

int detectEmptyComment ( char *str )
{
	char *cur;

	cur = str;
	if ( *cur == '*')
	{
		return ( 1 );
	}

	if ( str[0] == '\n' )
		return ( 1 );

	while ( *cur != '\0' )
	{
		if ( *cur != ' ' )
		{
			return ( 0 );
		}
		++cur;
	}

	return ( 1 );

}



int detectCommand ( char *str )
{
	//int i;
	int nofPlots = 0;

	double  *dcTable;
	char **dcNames, **plotNames;
	char *token;
	//char *buf, *buf1;

	//char delimiters [] = " \t";

	//char token_options[8] = ".options";
	//char token_dc[3] = ".dc";
	//char token_plot[5] = ".plot";



	if ( str[0] != '.')
	{
		if ( str[0] != '*' )
		{
			*str = toupper ( *str );
			if ( str[0] != 'R' && str[0] != 'V' && str[0] != 'I' && str[0] != 'L' && str[0] != 'C' && str[0] != 'D' && str[0] != 'M' && str[0] != 'Q' && str[0] != 'X' )
			{
				printf ( "ERROR : Wrong Command Declaration : \"%s\"\n", str );
				exit ( 0 );
			}
		}
		return 0;
	}

	str ++;

	token = strtok ( str, " ,=()\t\n" );
	if ( token != NULL )
	{
		if ( !strcmp ( token, "OPTIONS" ) )
		{
			token = strtok ( NULL, " =()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : .OPTIONS command without any parameter specified\n" );
				return 1;
			}
			if ( !strcmp ( token, "METHOD" ) )
			{
				token = strtok ( NULL, " =()\t\n" );
				if ( token == NULL )
				{
					printf ( "WARNING : No specific METHOD Specified. Default METHOD will be used.\n" );
					return 1;
				}
				if ( !configMethod ( token ) )
				{
					printf ( "WARNING : Unknown Method. Default METHOD will be used.\n" );
					return 1;
				}

			}
			if ( !strcmp ( token, "ITOL" ) )
			{
				token = strtok ( NULL, " =()\t\n" );
				if ( token == NULL )
				{
					printf ( "WARNING : No ITOL Value Specified. Default Value will be used.\n" );
					return 1;
				}
				cmdList->_itolValue = atof ( token );
				printf ( "\tINFO : set for ITOL value : %lf\n", cmdList->_itolValue );
				return 1;
			}
			if ( !strcmp ( token, "ITER" ) )
			{
				optionState[1] = 1;
				printf ( "\tINFO : set for ITER\n");
				configOptions ( optionState );
				return 1;
			}
			if ( !strcmp ( token, "SPD" ) )
			{
				optionState[0] = 1;
				configOptions ( optionState );
				printf ( "\tINFO : set for SPD\n");
				return 1;
			}
			return 1;
		}

		if ( !strcmp ( token, "DC" ) )
		{
			dcTable = ( double * )malloc ( 3 * sizeof ( double ) );
			dcNames = ( char ** )malloc ( sizeof ( char *) );

			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : DC Command without any parameters specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}
			dcNames[0] = strdup ( token );
			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : DC Command without <start_value> parameter specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}
			dcTable[0] = atof ( token );
			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : DC Command without <end_value> parameter specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}
			dcTable[1] = atof ( token );
			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : DC Command without <increment> parameter specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}
			dcTable[2] = atof ( token );
			printf ( "\tINFO : set for DC Sweep\n");

			configDcSweep( dcTable, dcNames , 1 );
			return 1;
		}

		if ( !strcmp ( token, "TRAN" ) )
		{
			double timeStep;

			token = strtok ( NULL, " =()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : TRAN Command without any parameters specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}

			timeStep = atof ( token );
			token = strtok ( NULL, " =()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : TRAN Command without <fin_time> parameter specified.\n\t>Command will be ignored -May lead to Simulation Errors.\n" );
				return 1;
			}

			cmdList->_tranTimeStep = timeStep;
			cmdList->_tranFinTime = atof ( token );
			cmdList->_transient = 1;
			printf ( "\tINFO : set for Transient Analysis\n\t\t-time step : %lf\n\t\t-finish time : %lf\n", cmdList->_tranTimeStep, cmdList->_tranFinTime );

			return 1;
		}

		if ( !strcmp ( token, "PLOT" ) )
		{
			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : PLOT Command without any parameters specified.\n\t>No data will be plotted.\n" );
				return 1;
			}
			if ( strcmp ( token, "v" ) )
			{
				printf ( "WARNING : PLOT Command without info what to plot.\n\t>No data will be plotted.\n" );
				return 1;
			}
			token = strtok ( NULL, " ,=()\t\n" );
			if ( token == NULL )
			{
				printf ( "WARNING : PLOT Command without any nodes specified.\n\t>No data will be plotted.\n" );
				return 1;
			}

			cmdList->_plot = ( char ** )malloc ( sizeof ( char *) );

			while ( token != NULL )
			{
				cmdList->_plot[nofPlots] = strdup ( token );
				cmdList->_plot = ( char ** )realloc ( cmdList->_plot, ( nofPlots + 1 ) * sizeof ( char *) );
				++nofPlots;
				cmdList->_plotnum = nofPlots;
				token = strtok ( NULL, " ,=()\t\n" );
			}
			//configPlot ( plotNames, nofPlots );

			return 1;
		}

	}

	return 0;

}



/**
 * classifyElement
 *
 * Parses the string <buf> and decides wether it should create a linear or
 * a non linear element struct and fills it appropriately.
 * Once done, links the structure to the main list structure <element> and updates
 * the type field of the later.
 *
 * @param element: main list structure
 * @param buf: input string
 * @return void
 */

void classifyElement ( zeroCircuit *element, char *buf, hashTable *nameTable  )
{
	char delimiters [] = " \t";
	char *running;
	char *token;
	char *cur;
	struct linear_element *linStruct;
	struct nonlinear_element *nonlinStruct;

	int res;

	running = strdup ( buf );

	token = strtok ( buf, " ()\t\n" );
	cur = token;
	*cur = toupper ( *cur );


	if( *cur == 'R' || *cur == 'V' || *cur == 'I' || *cur == 'C' || *cur == 'L' )
	{
		if(*cur=='V') {
			vnumber++;
		}

		if(*cur=='L') {
			lnumber++;
		}

		if(*cur=='R') {
			rnumber++;
		}

		if(*cur=='I') {
			inumber++;
		}

		//allocate linear structure
		linStruct = ( struct linear_element * )malloc ( sizeof ( struct linear_element ) );
		linStruct->vInfo = NULL;

		linStruct->id = strdup ( token );

		token = strtok ( NULL, " ()\t\n" );
		linStruct->connectors[0] = strdup ( token );
		res = insertHash ( nameTable, token );


		token = strtok ( NULL, " ()\t\n" );
		linStruct->connectors[1] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		linStruct->value = atof ( token );

		token = strtok ( NULL, " ()\t\n" ); //for later use on transient parsing

		if (*cur == 'R')
			linStruct->element = R;
		if (*cur == 'V')
		{
			int result;

			linStruct->element = V;
			result = classifyVsource ( running, linStruct );
		}
		if (*cur == 'I')
		{
			int result;

			linStruct->element = I;
			result = classifyVsource ( running, linStruct );
		}
		if (*cur == 'C')
			linStruct->element = C;
		if (*cur == 'L')
			linStruct->element = L;

		//fix wrapper struct
		element->type = LIN;
		element->linElement = linStruct;

	}
	else if ( *cur == 'D' )
	{
		nonlinStruct = ( struct nonlinear_element * )malloc ( sizeof ( struct nonlinear_element ) );
		nonlinStruct->id = strdup ( token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[0] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[1] = strdup ( token );
		res = insertHash ( nameTable, token );

		nonlinStruct->connectors[2] = NULL;

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->area = atof ( token );

		nonlinStruct->element = D;

		//fix wrapper struct
		element->type = NONLIN;
		element->nonlinElement = nonlinStruct;
	}
	else if ( *cur == 'M' )
	{
		nonlinStruct = ( struct nonlinear_element * )malloc ( sizeof ( struct nonlinear_element ) );
		nonlinStruct->id = strdup ( token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[0] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[1] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[2] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strsep (&running, delimiters);
		nonlinStruct-> length = atof ( token );

		token = strsep (&running, delimiters);
		nonlinStruct-> width = atof ( token );

		nonlinStruct->element = M;

		//fix wrapper struct
		element->type = NONLIN;
		element->nonlinElement = nonlinStruct;
	}
	else if ( *cur == 'Q' )
	{
		nonlinStruct = ( struct nonlinear_element * )malloc ( sizeof ( struct nonlinear_element ) );
		nonlinStruct->id = strdup ( token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[0] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[1] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->connectors[2] = strdup ( token );
		res = insertHash ( nameTable, token );

		token = strtok ( NULL, " ()\t\n" );
		nonlinStruct->area = atof ( token );

		nonlinStruct->element = Q;

		//fix wrapper struct
		element->type = NONLIN;
		element->nonlinElement = nonlinStruct;
	}


}



/**
 *
 */

zeroCircuit *parseNetlist ( char *filename, hashTable *namTable )
{

	zeroCircuit	*netlist;
	zeroCircuit *current;

	char lineBuf [ BUFSIZE ];

	FILE *fileObj;

	//initialize circuit list
	netlist = ( zeroCircuit * )malloc( sizeof ( zeroCircuit ) );
	netlist->prev = netlist;
	netlist->next = netlist;
	netlist->linElement = NULL;
	netlist->nonlinElement = NULL;
	netlist->type = LIN;	//no real meaning


	fileObj =  fopen ( filename , "r" );
	if ( fileObj == NULL )
	{
		printf( "Couldn't read file\n" );
	}


	while ( fgets ( lineBuf, BUFSIZE, fileObj ) )
	{
			if( detectEmptyComment( lineBuf ) )
			{
				continue;
			}


			if ( detectCommand( lineBuf ) )
			{
				continue;
			}


			//create and add node of the circuit
			current = ( zeroCircuit * )malloc( sizeof ( zeroCircuit ) );
			current->next = netlist->next;
			current->prev = netlist;
			current->next->prev = current;
			current->prev->next = current;   //bug fixed!
			classifyElement ( current, lineBuf, namTable );
	}

	fclose ( fileObj );

	return netlist;
}




