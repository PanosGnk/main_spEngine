/*
 * transientParse.h
 *
 *  Created on: Jan 22, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */

#ifndef TRANSIENTPARSE_H_
#define TRANSIENTPARSE_H_


#include "types.h"
#include <string.h>
#include <stdlib.h>

int classifyVsource ( char *str,struct linear_element *linElement );
int finiteParametersParse ( char *str, int nofParams, struct linear_element *linElement, int type );
int classifyVsource ( char *str, struct linear_element *linElement );


#endif /* TRANSIENTPARSE_H_ */
