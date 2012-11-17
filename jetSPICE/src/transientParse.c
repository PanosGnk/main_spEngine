/*
 * transientParse.c
 *
 *  Created on: Jan 22, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */

#include "transientParse.h"




//enum _infotype
//{
//	NONDEF 	= 0,
//	EXP 	= 1,
//	SIN 	= 2,
//	PULSE	= 3,
//	PWL 	= 4
//};


int classifyVsource ( char *str, struct linear_element *linElement )
{
	char *token;
	char *buf;
	int res;
	int i;

	double *vData;
	vSourceInfo *srcInf;

	buf = strdup ( str );
	token = strtok ( str, " ()\t\n" );
	for ( i=0; i < 4; i++ )
	{
		token = strtok ( NULL, " ()\t\n" );
	}
	if ( token == NULL )
	{
		//printf ( "INFO : Voltage Source without \"[transient_spec]\" defined" );
		return 0;
	}

	printf ( "\tINFO : Found Source with transient description : %s", buf );
	free ( buf );

		if ( !strcmp ( token, "EXP" ) )
		{
			//parameters to be parsed : i1 i2 td1 tc1 td2 tc2
			vData = ( double * )malloc ( 6 * sizeof ( double ) );
			srcInf = ( vSourceInfo * )malloc ( sizeof ( vSourceInfo ) );

			for ( i=0; i < 6; i++ )
			{
				token = strtok ( NULL, " ()\t\n" );
				if ( token == NULL )
				{
					printf ( "ERROR : Missing transient parameters. Exiting\n" );
					exit ( 0 );
				}
				vData[i] = atof ( token );
			}

			srcInf->pVal = linElement->value;
			srcInf->_data = vData;
			srcInf->_infoType = EXP;
			linElement->vInfo = srcInf;
			linElement->vInfo->_nofData = 6;

		}

		if ( !strcmp ( token, "SIN" ) )
		{
			//parameters to be parsed : i1 ia fr td df ph
			vData = ( double * )malloc ( 6 * sizeof ( double ) );
			srcInf = ( vSourceInfo * )malloc ( sizeof ( vSourceInfo ) );

			for ( i=0; i < 6; i++ )
			{
				token = strtok ( NULL, " ()\t\n" );
				if ( token == NULL )
				{
					printf ( "ERROR : Missing transient parameters. Exiting\n" );
					exit ( 0 );
				}
					vData[i] = atof ( token );
			}

			srcInf->pVal = linElement->value;
			srcInf->_data = vData;
			srcInf->_infoType = SIN;
			linElement->vInfo = srcInf;
			linElement->vInfo->_nofData = 6;
		}

		if ( !strcmp ( token, "PULSE" ) )
		{
			//parameters to be parsed : i1 i2 td tr tf pw per
			vData = ( double * )malloc ( 7 * sizeof ( double ) );
			srcInf = ( vSourceInfo * )malloc ( sizeof ( vSourceInfo ) );

			for ( i=0; i < 7; i++ )
			{
				token = strtok ( NULL, " ()\t\n" );
				if ( token == NULL )
				{
					printf ( "ERROR : Missing transient parameters. Exiting\n" );
					exit ( 0 );
				}
				vData[i] = atof ( token );
			}

			srcInf->pVal = linElement->value;
			srcInf->_data = vData;
			srcInf->_infoType = PULSE;
			linElement->vInfo = srcInf;
			linElement->vInfo->_nofData = 7;
		}

		if ( !strcmp ( token, "PWL" ) )
		{
			int parsed;

			parsed = 1;

			vData = ( double * )malloc ( parsed * sizeof ( double ) );
			srcInf = ( vSourceInfo * )malloc ( sizeof ( vSourceInfo ) );
			//parameters to be parsed : (t1 i1) ... (tn in)

			token = strtok ( NULL, " ()\t\n" );
			if ( token == NULL )
			{
				printf ( "ERROR : Missing transient parameters. Exiting\n" );
				exit ( 0 );
			}

			while( token != NULL )
			{
				vData[parsed -1] = atof ( token );
				token = strtok ( NULL, " ()\t\n" );
				if ( ( token == NULL ) && ( parsed%2 == 1 ) )
				{
					printf ( "ERROR : Missing transient parameters. Exiting\n" );
					exit ( 0 );
				}
				parsed ++;
				vData = ( double * )realloc ( vData, parsed * sizeof ( double ) );
			}

			srcInf->pVal = linElement->value;
			srcInf->_data = vData;
			srcInf->_infoType = PWL;
			linElement->vInfo = srcInf;
			linElement->vInfo->_nofData = parsed - 1;
		}
		return res;
}


