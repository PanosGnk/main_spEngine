
#include "sim.h"
#include <stdio.h>
#include "types.h"
#include "mna.h"
#include "linearSolve.h"
#include "mnaSparse.h"



commandList *cmdList;

void initializeCommandList ()
{

	cmdList = ( commandList * )malloc ( sizeof ( commandList ) );

	cmdList->_dcValues = NULL;
	cmdList->_dcId = NULL;
	cmdList->_dcnum = 0;

	cmdList->_options = DEFAULT;
	cmdList->_method = NODEF;
	cmdList->_itolValue = 1.0e-3;	//default itol value (maybe a routine to change default values could be good)

	cmdList->_transient = 0;
	cmdList->_tranFinTime = 0.0;
	cmdList->_tranTimeStep = 0.0;

	cmdList->_plot = NULL;
	cmdList->_plotnum = 0;

	cmdList->_useSparse = 0;

}



void setMemoryPolicy ( hashTable *namTab, int lnum, int vnum )
{
	int dim;

	dim = namTab->currPos + lnum + vnum;

	if ( dim < DIM_THRESHOLD )
		cmdList->_useSparse = 0;
	else
		cmdList->_useSparse = 1;
}



void configOptions ( int *options )
{
	if ( options[0] && options[1] )
	{
		cmdList->_options = ITERSPD;
		return;
	}

	if ( options [0] )
	{
		cmdList->_options = SPD;
		return;
	}

	if ( options [1] )
	{
		cmdList->_options = ITER;
		return;
	}
}




void configDcSweep ( double *values, char **srcs, int num )
{

	cmdList->_dcnum = num;
	cmdList->_dcValues = values;
	cmdList->_dcId = srcs;

	printf ( "\t\t-number of sources : %d\n\t\t-name : %s\n\t\t-start value : %lf\n\t\t-end value : %lf\n\t\t-increment : %lf\n", cmdList->_dcnum, cmdList->_dcId[0], cmdList->_dcValues[0],cmdList->_dcValues[1],cmdList->_dcValues[2] );
}



void configPlot( char **plots, int num )
{
	cmdList->_plot = plots;
	cmdList->_plotnum = num;

	printf ( "\tINFO : Nodes to be plotted:\n");

	int i;

	for ( i=0; i < cmdList->_plotnum; i++)
	{
		printf ( "\t\t-node : %s\n", cmdList->_plot[i]);
	}
}



int configMethod ( char *str )
{
	if ( !strcmp ( str, "TR\0" ) )
	{
		cmdList->_method = TR;
		printf ( "\tINFO : set METHOD : Trapezoidal\n");
		return 1;
	}

	if ( !strcmp ( str, "BE\0" ) )
	{
		cmdList->_method = BE;
		printf ( "\tINFO : set METHOD : Backward Euler\n");
		return 1;
	}

	return 0;
}





/*
void configItol ( double tolerance )
{
	cmdList->_itolValue = tolerance;
}
*/


traceFileList *initializeTraceFiles ()
{
	int nofiles;
	traceFileList *files;

	nofiles = cmdList->_plotnum;

	files = ( traceFileList * )malloc ( sizeof ( traceFileList ) );

	files->trFiles = ( FILE ** )malloc ( nofiles * sizeof ( FILE * ) );
	files->size = nofiles;

	printf ( "Initialized file List...\n" );

	return files;
}

void initializeFiles ( traceFileList *flist )
{
	int i;
	char *str;

	//printf ( "List size :%d...\n", flist->size );

	for ( i = 0; i < flist->size; i++ )
	{
		str = (char *)malloc(32 * sizeof(char));

		sprintf(str, "trace_node(%s).tr", cmdList->_plot[i] );

		//printf ( ">>>>>%s\n",cmdList->_plot[i]);

		flist->trFiles[i] = fopen( str,"w" );
		if ( flist->trFiles[i]== NULL )
		{
			printf("Can't Create File ... \n");
			exit(0);
		}
	}

	//printf ( "Opened trace files...\n" );
}



int *findVarPositions ( char **vec, int size )
{
	int *posTable;
	int i,j;

	posTable = (int *)malloc ( cmdList->_plotnum * sizeof ( int ) );

	for ( i = 0; i < size; i++ )
	{
		if ( vec[i][0] == 'V' )
		{
			for (j = 0; j < cmdList->_plotnum; j++)
			{
				if ( vec[i][7] == *cmdList->_plot[j] )
				{
					posTable[j] = i;
				}
			}
		}
	}

	return posTable;
}



int *findPositions ( hashTable *namTab )
{
	int *positions;
	int i;

	positions = ( int * )malloc ( cmdList->_plotnum * sizeof (int));

	for ( i = 0; i < cmdList->_plotnum; i++ )
	{
		positions[i] = findHash ( namTab, cmdList->_plot[i]);
		if ( positions[i] == -1 )
		{
			printf ( "ERROR : Specified node '%s' does not exhist\n", cmdList->_plot[i] );
			exit  (-1);
		}
		//printf ("_________>%d\n",positions[i]);
	}

	return ( positions );
}


int checkDCPos (  zeroCircuit *circuit )
{
	int pos;
	zeroCircuit *cur;

	for( cur = circuit->next; ( (cur->linElement != NULL) || (cur->nonlinElement != NULL) ); cur = cur->next )
	{
		if ( !strcmp ( cur->linElement->id, cmdList->_dcId[0] ) )
		{
			pos = cur->linElement->k;
			return pos;
		}
		
	}
	
	printf ( "ERROR : Specified node '%s' does not exist\n", cmdList->_dcId[0] );
	exit  (-1);
}

int checkDCPosSp ( zeroCircuit *circuit )
{
	int pos;
	zeroCircuit *cur;

	for( cur = circuit->next; ( (cur->linElement != NULL) || (cur->nonlinElement != NULL) ); cur = cur->next )
	{
		if ( !strcmp ( cur->linElement->id, cmdList->_dcId[0] ) )
		{
			pos = cur->linElement->k;
			return pos;
		}
		
	}
	
	printf ( "ERROR : Specified node '%s' does not exist\n", cmdList->_dcId[0] );
	exit  (-1);
}




///****************************************************************************************
// * 								MAIN SIMULATION FUNCTIONS                               *
// ****************************************************************************************/
//
//
//
///*
// ==================================
// Linear (Normal) Storage Functions
// ==================================
// */



void dcPointSim ( mnaSystem *sys )
{
	double *pVec;
	double *result;
	int *pivotTab;
	double **newPiv;


	int i,j;
	int iters=0;

	switch (cmdList->_options)	//use cholesky method
	{
	case SPD:
		printf("\tINFO : Used Cholesky Factorization\n");
		pVec = ( double * )malloc ( sys->dim * sizeof( double ) );

		choleskyDecomp ( sys->rowVector, sys->dim, pVec );
		result = choleskySolve ( sys->rowVector, sys->dim, pVec, sys->vector_B );
		break;

	case DEFAULT:
	    printf("\tINFO : Used LU Factorization\n");
	    pivotTab = ( int * )calloc( sys->dim , sizeof ( int ) );
		newPiv = luDecomp ( sys->rowVector, sys->dim, pivotTab );

		result = luBackSubst( newPiv, sys->dim, pivotTab, sys->vector_B);
		break;

	case ITER:
		printf("\tINFO : Used Bi-CG Method\n");
		result = bicg (sys->rowVector ,sys->vector_B ,sys->dim,&iters ,cmdList->_itolValue );
		printf("\tINFO : Iterations : %d\n",iters);
		break;

	case ITERSPD:
		printf("\tINFO : Used CG Method\n");
		result = cg (sys->rowVector ,sys->vector_B ,sys->dim,&iters ,cmdList->_itolValue );
		printf("\tINFO : Iterations : %d\n",iters);
		break;
	}

	sys->dcPointRes = result;
		printf("\n\n");

		printf ("\t--- Solution ---\n");

		for ( i=0; i < sys->dim; i++)
		{
			printf("\t    %lf\n", result[i]);
		}

		printf("\n\n");
/*
		printf ("\t--- L / U ---\n");
		for ( i=0; i < sys->dim; i++)
		{
			printf("\n");
			for ( j=0; j < sys->dim; j++ )
			{
				printf("\t    %.10lf", *(newPiv[i] + j));
			}
		}
*/

}



void dcSweepSim ( mnaSystem *sys, traceFileList *flist, hashTable *names, zeroCircuit *cir )
{

	//int nofSims;
	int k;
	double i;
	int pos;

	double *pVec;
	double *result;
	int *pivotTab;
	double **newPiv;
	char **vars;

	int *positions;

	int iters;

	iters = 0;

	switch ( cmdList->_options )
	{
	case DEFAULT:

		printf ("In DEFAULT case....\n\n");
		pivotTab = ( int * )calloc( sys->dim , sizeof ( int ) );
		newPiv = luDecomp ( sys->rowVector, sys->dim, pivotTab );

		positions = findPositions ( names );
		pos = checkDCPos ( cir );

		for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
		{
			sys->vector_B[pos] = i;
 			result = luBackSubst( newPiv, sys->dim, pivotTab, sys->vector_B);
 			if ( cmdList->_plotnum )
 			{
 				int m;

 				for ( m = 0; m < flist->size; m++ )
 				{
 					fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result [positions[m]]  );
 				}
 			}
		}

		printf( "\n\n" );

		break;

	case SPD:

		//printf ("In SPD case....\n\n");
		positions = findPositions ( names );
		pos = checkDCPos ( cir );

		pVec = ( double * )malloc ( sys->dim * sizeof( double ) );
		choleskyDecomp ( sys->rowVector, sys->dim, pVec );



		for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
		{
			sys->vector_B[pos] = i;
			result = choleskySolve ( sys->rowVector, sys->dim, pVec, sys->vector_B );

			if ( cmdList->_plotnum )
			{
				int m;

				for ( m = 0; m < flist->size; m++ )
			 	{
					//printf ( "Editing file...\n" );
			 		fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
				}
			 }
		}

		result = choleskySolve ( sys->rowVector, sys->dim, pVec, sys->vector_B );
		break;

	case ITER:

			//printf ("In ITER case....\n\n");
			positions = findPositions ( names );
			pos = checkDCPos ( cir );

	        printf( "\n\n" );

	        for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
	        {

	            sys->vector_B[pos] = i;

	            result = bicg (sys->rowVector ,sys->vector_B ,sys->dim,&iters ,cmdList->_itolValue );

	            if ( cmdList->_plotnum )
	                  {
	                        int m;

	                        for ( m = 0; m < flist->size; m++ )
	                        {
	                            fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
	                        }
	                  }

	         }

	         printf( "\n\n" );

	         break;

	case ITERSPD:

			positions = findPositions ( names );
			pos = checkDCPos ( cir );

	         //printf ("In ITERSPD case....\n\n");
	         printf( "\n\n" );

	         for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
	         {

	               sys->vector_B[pos] = i;
	               result = cg (sys->rowVector ,sys->vector_B ,sys->dim ,&iters ,cmdList->_itolValue );
	               printf("CG DONE...\n");

	               if ( cmdList->_plotnum )
	                     {
	                           int m;

	                           for ( m = 0; m < flist->size; m++ )
	                           {
	                               printf ( "Editing file...\n" );
	                               fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
	                           }
	                      }


	               printf ("\t--- Solution ---\n");
	               for ( k=0; k < sys->dim; k++)
	               {
	            	   printf("\t    %lf\n", result[k]);
	               }
	         }
	         printf( "\n\n" );

	         break;
	}
}



/*
 ===============================================
 Sparse Storage Functions (csparse library used)
 ===============================================
 */



void dcPointSimSparse ( mnaSpSystem *sys )
{
	double *result;
	int i, status;
	int iters=0;

	if( !( result = (double *)calloc(sys->dim, sizeof(double)) ) )
	{
	    printf("ERROR : Memory full ( allocate failure )\n");
	    exit(-1);
	}

	switch (cmdList->_options)
	{
	case SPD:
		printf("\tINFO : Used Cholesky Factorization\n");
		vecCpy(result, sys->vector_B, sys->dim);
		status = cs_cholsol(1, sys->array_A, result);
		break;

	case DEFAULT:
		printf("\tINFO : Used LU Factorization\n");
                vecCpy(result, sys->vector_B, sys->dim);
		status = cs_lusol(2, sys->array_A, result, 1.0);
		break;

	case ITER:
		printf("\tINFO : Used Bi-CG Method\n");
		result = cs_bicg (sys->array_A ,sys->vector_B ,sys->dim ,&iters ,cmdList->_itolValue );
		printf("\tINFO : Iterations : %d\n",iters);
		break;

	case ITERSPD:
		printf("\tINFO : Used CG Method\n");
		result = cs_cg (sys->array_A ,sys->vector_B ,sys->dim ,&iters ,cmdList->_itolValue );
		printf("\tINFO : Iterations : %d\n",iters);
		break;
	}

	sys->dcPointRes = result;


	FILE *bfile;

	bfile = fopen ( "solution_vector.txt", "w" );

	fprintf ( bfile, "Result vector\nDC Fixed Point Analysis\n\n" );
	for ( i=0; i<( sys->dim ); i++ )
	{
		fprintf ( bfile, "node(%d)\t\t\t%.10f\n", i, result[i] );
	}

	fclose ( bfile );

}



void dcSweepSimSparse ( mnaSpSystem *sys, traceFileList *flist, hashTable *names, zeroCircuit *cir )
{

	//int nofSims;
	int k;
	double i;
	int pos;

	double *pVec;
	double *result;
	int *pivotTab;
	double **newPiv;
	char **vars;

	int *positions;

	int iters, state;

	iters = 0;
        if( !( result = (double *)calloc(sys->dim, sizeof(double)) ) )
        {
            printf("ERROR : Memory full ( allocate failure )\n");
            exit(-1);
        }

	switch ( cmdList->_options )
	{
	case DEFAULT:

		//BUGFIX : The X vector indicating the system variables is not transversed
		positions = findPositions ( names );
		pos = checkDCPosSp ( cir );

		printf( "\n\n" );

		for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
		{
		        vecCpy(result, sys->vector_B, sys->dim);
			result[pos] = i;

 			state = cs_lusol(2, sys->array_A, result, 1.0);

 			if ( cmdList->_plotnum )
 			{
 				int m;

 				for ( m = 0; m < flist->size; m++ )
 				{
 					//printf ( "Editing file...\n" );
 					//fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\t%.10lf\n", cmdList->_dcId[0], i, result[ positions[ m ]  ]  );
 					fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
 				}
 			}


 			//printf ("\t--- Solution ---\n");

 				//	for ( k=0; k < sys->dim; k++)
 					//{
 						//printf("\t    %lf\n", result[k]);
 					//}
		}

		printf( "\n\n" );

		break;

	case SPD:

		//printf ("In SPD case....\n\n");
		positions = findPositions ( names );
		pos = checkDCPosSp ( cir );

		for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
		{

		        vecCpy(result, sys->vector_B, sys->dim);
			result[pos] = i;

			state =  cs_cholsol(1, sys->array_A, result);


			if ( cmdList->_plotnum )
			{
				int m;

				for ( m = 0; m < flist->size; m++ )
			 	{
					//printf ( "Editing file...\n" );
			 		fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
				}
			 }
		}


		result = choleskySolve ( sys->rowVector, sys->dim, pVec, sys->vector_B );
		break;

	case ITER:

			//printf ("In ITER case....\n\n");
			positions = findPositions ( names );
			pos = checkDCPosSp ( cir );

	        printf( "\n\n" );

	        for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
	        {

	            vecCpy(result, sys->vector_B, sys->dim);
	            result[pos] = i;

	            result = cs_bicg (sys->array_A ,result ,sys->dim ,&iters ,cmdList->_itolValue );

	            if ( cmdList->_plotnum )
	                  {
	                        int m;

	                        for ( m = 0; m < flist->size; m++ )
	                        {
	                            //printf ( "Editing file...\n" );
	                            fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
	                        }
	                  }
	         }
	         printf( "\n\n" );

	         break;

	case ITERSPD:

	         //printf ("In ITERSPD case....\n\n");
			positions = findPositions ( names );
			pos = checkDCPosSp ( cir );
	         printf( "\n\n" );

	         for ( i = cmdList->_dcValues[0]; i <= cmdList->_dcValues[1]; i += cmdList->_dcValues[2] )
	         {

	             vecCpy(result, sys->vector_B, sys->dim);
	             result[pos] = i;

	               result = cs_cg (sys->array_A ,result ,sys->dim ,&iters ,cmdList->_itolValue );

	               printf("CG DONE...\n");

	               if ( cmdList->_plotnum )
	                     {
	                           int m;

	                           for ( m = 0; m < flist->size; m++ )
	                           {
	                               printf ( "Editing file...\n" );
	                               fprintf ( flist->trFiles[m], "> %s :%.10lf\t\t\tnode %s :%.10lf\n", cmdList->_dcId[0], i, cmdList->_plot[m], result[ positions[ m ]  ]  );
	                           }
	                      }

	                printf ("\t--- Solution ---\n");
                    	for ( k=0; k < sys->dim; k++)
                    	{
                    		printf("\t    %lf\n", result[k]);
                    	}
	         }

	         printf( "\n\n" );

	         break;


	}
}





