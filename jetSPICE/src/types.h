/*
 * types.h
 *
 *  Created on: Oct 8, 2011
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 *     Version: 0.2
 */

#ifndef TYPES_H_
#define TYPES_H_


#include <stdio.h>



#define DIM_THRESHOLD 2000



int vnumber;
int lnumber;
int rnumber;
int inumber;

enum element_id
{
	R = 1,
	V = 2,
	I = 3,
	C = 4,
	L = 5,
	D = 6,
	M = 7,
	Q = 8
};


enum element_type
{
	LIN = 0,
	NONLIN = 1
};


typedef struct model_definition
{
	//TODO model definition values here

}modeldef;



enum _infotype
{
	NONDEF 	= 0,
	EXP 	= 1,
	SIN 	= 2,
	PULSE	= 3,
	PWL 	= 4

};

struct vSource
{
	int _infoType;
	int	_nofData;
	double	*_data;
	double pVal;
};

typedef struct vSource vSourceInfo;


// Struct to contain the linear circuit elements

struct linear_element
{
	int element;
	char *id;
	char *connectors [2];		//Polarity: 0= node+, 1= node-

	double value;

	int k;
	vSourceInfo *vInfo;
};


struct nonlinear_element
{
	int element;
	char *id;
	char *connectors [3];		//For MOS: 		0=D, 1=G. 2=S
					//For BJT: 		0=C, 1=B, 2=E
					//For Diodes: 	        0= node+, 1= node-

	modeldef *model;		//FIXME: Change type on a more appropriate one

	double area;			//For diodes/BJT
	double width;			//For MOS
	double length;			//For MOS
};



struct elementList
{
	int type;

	struct linear_element *linElement; //typecast to the right type of element
	struct nonlinear_element *nonlinElement;

	struct elementList *prev;
	struct elementList *next;


	//Points to the next/previous element of the same type (element_id). Used for skip list integration
	struct elementList *nextSimilar;
	struct elementList *prevSimilar;
};

typedef struct elementList zeroCircuit;


enum options
{
	DEFAULT = 0,
	SPD 	= 1,
	ITER	= 2,
	ITERSPD = 3
};



enum method
{
	NODEF = 0,
	TR = 1,
	BE = 2
};



struct _commandList
{
	int  	_options;			//keep the option in
	double  _itolValue;			//keep values declared in options
	int	_method;

	double  *_dcValues;			//keep start_values, end_values, increment_steps
	char	**_dcId;			//IDs of the voltage sources for the dc sweep
	int	_dcnum;

	int     _transient;
	double	_tranTimeStep;		//Values necessary for transitional analysis
	double 	_tranFinTime;

	char   	**_plot;			//variable names to be ploted
	int	_plotnum;			//number of variables. If 0 then no sim.

	int		_useSparse;		//change memory policy

};

typedef struct _commandList commandList;


/*********** File Handling Structures **************/


//Optional structures in case of handling multiple files

struct parse_element
{
	char *	ParsedFile;
	FILE *  FileReader;

	int parsedLine;
};


struct _traceFileList
{
	FILE 	**trFiles;
	int	size;
};

typedef struct _traceFileList traceFileList;




#endif /* TYPES_H_ */
