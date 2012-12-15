/*
 * mna.c
 *
 *  Created on: Jan 18, 2012
 *      Author: engineer
 */



#include "mna.h"
#include <stdio.h>



extern int vnumber;
extern int lnumber;
extern int rnumber;
extern int inumber;

extern commandList *cmdList;

nodeTable initializeTable ()
{
  nodeTable nodes;

  nodes.size = 0;
  nodes.table = NULL;

  return nodes;
}



/*
   int findPlace(int connector,nodeTable *table )
   {
   int i;

   for (i=0; i<table->size; i++) {

   if(table->table[i]==connector) {
   return i;
   }

   }
   return -1;
   }
 */



/*
   nodeTable  analyseNodes( zeroCircuit *node, nodeTable nodes )
   {
   int *cur;
   int iter;
   int found1, found2, found3;


   found1 = 0;
   found2 = 0;
   found3 = 0;


   switch ( node->type )
   {
   case 0:

   for ( cur = nodes.table , iter = 0; iter < nodes.size; ++iter, ++cur )
   {
   if ( *cur == node->linElement->connectors[0] )
   {
   found1 = 1;
   }
   if ( *cur == node->linElement->connectors[1] )
   {
   found2 = 1;
   }

   }

   if ( !found1 )
   {
   if ( node->linElement->connectors[0] != 0 )
   {
   nodes.table = realloc( nodes.table, ( nodes.size + 1 )*sizeof( int ) );
   nodes.size++;

   nodes.table [nodes.size-1] = node->linElement->connectors[0];
   }
   }

   if ( !found2 )
   {
   if ( node->linElement->connectors[1] != 0 )
   {
   nodes.table = realloc( nodes.table, ( nodes.size + 1 )*sizeof( int ) );
   nodes.size++;

   nodes.table [ nodes.size-1 ] = node->linElement->connectors[1];
   }
   }

   break;
   case 1:

// Code for non-linear elements handling

break;
}

return nodes;

}



void printNodes( nodeTable nodes )
{
int i;

printf( "LookUp table of nodes:\n\t-" );

for ( i=0; i< nodes.size; i++)
{
  printf( "%d ", nodes.table[i] );
}
printf( "\n" );

}
*/


/*
 *  R = 1,
 V = 2,
 I = 3,
 C = 4,
 L = 5,
 D = 6,
 M = 7,
 Q = 8
 */


// #define STRSIZE 32


void initVectorX ( mnaSystem *system, hashTable *table  )
{
  //char *str;
  int i;

  for ( i=0; i<table->currPos; i++ )
  {
    //str = NULL;
    //str = (char *)malloc( STRSIZE * sizeof( char ) );

    //sprintf(str, "V(node %d)", findHas () );

    //system->vector_X[i] = str;
  }
}



void fillVectorX ( zeroCircuit *node, mnaSystem *system, int k, int origsize )
{
  char *str;

  str = (char *)malloc( STRSIZE * sizeof( char ) );

  sprintf( str, "I(%s)", node->linElement->id );

  system->vector_X[ origsize + k -1 ] = str;
}



void addElementStamp( zeroCircuit *node, mnaSystem *system, hashTable *namTab, int row, int kPos )
{
  int posn1;
  int posn2;

  struct linear_element *element;

  element = node->linElement;


  switch ( node->type )
  {
    case 0:
      switch( node->linElement->element )
      {
        case 1:

          if ( ( strcmp( node->linElement->connectors[0], "0" ) ) && ( strcmp( node->linElement->connectors[1], "0" ) ) )
          {
            //code for simple integer mapping
            //posn1 = findPlace ( node->linElement->connectors[0], &nodes );
            //posn2 = findPlace ( node->linElement->connectors[1], &nodes );

            //code for hashtable string mapping
            posn1 = findHash ( namTab, node->linElement->connectors[0] );
            //printf ( "----> %d", posn1 );
            posn2 = findHash ( namTab, node->linElement->connectors[1] );
            //printf ( "----> %d", posn2 );

            system->array_A[ row * posn1 + posn2 ] += (-1.0) * ( 1.0 / node->linElement->value);
            system->array_A[ row * posn2 + posn1 ] += (-1.0) * ( 1.0 / node->linElement->value);


          }

          if ( strcmp( node->linElement->connectors[0], "0\0" ) )
          {
            //posn1 = findPlace ( node->linElement->connectors[0], &nodes );
            posn1 = findHash ( namTab, node->linElement->connectors[0] );
            system->array_A[ row * posn1 + posn1 ] += ( 1.0 / node->linElement->value );
          }

          if ( strcmp( node->linElement->connectors[1], "0\0" ) )
          {
            //posn2 = findPlace ( node->linElement->connectors[1], &nodes );
            posn2 = findHash ( namTab, node->linElement->connectors[1] );
            system->array_A[ row * posn2 + posn2 ] += ( 1.0 / node->linElement->value );

            //printf( "bug! %f\n", system->array_A[ row * posn2 + posn2 ]);
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
          }


          if ( ( strcmp( node->linElement->connectors[0], "0\0" ) ) && ( strcmp( node->linElement->connectors[1], "0\0" ) ) )  // implement with strcmp()
          {
            //posn1 = findPlace ( element->connectors[0], &nodes );
            //posn2 = findPlace ( element->connectors[1], &nodes );

            posn1 = findHash ( namTab, element->connectors[0] );
            posn2 = findHash ( namTab, element->connectors[1] );

            system->array_A[ row * (kPos + namTab->currPos -1 ) + posn1 ] +=  1.0;
            system->array_A[ row * (kPos + namTab->currPos -1 ) + posn2 ] += -1.0;

            system->array_A[ row * posn1 + kPos + namTab->currPos - 1 ] +=  1.0;
            system->array_A[ row * posn2 + kPos + namTab->currPos - 1 ] += -1.0;

            break;
          }

          if ( strcmp( node->linElement->connectors[0], "0\0" ) )
          {            //posn2 = findPlace ( element->connectors[1], &nodes );
            posn2 = findHash ( namTab, element->connectors[0] );

            system->array_A[ row * (kPos + namTab->currPos - 1 ) + posn2 ] += 1.0;
            system->array_A[ row * posn2 + kPos + namTab->currPos - 1 ] += 1.0;

            break;
          }

          if ( strcmp( node->linElement->connectors[1], "0\0" ) )
          {
            //posn1 = findPlace ( element->connectors[0], &nodes );
            posn1 = findHash ( namTab, element->connectors[1] );

            system->array_A[ row * (kPos + namTab->currPos - 1 ) + posn1  ] +=  -1.0;
            system->array_A[ row * posn1 + kPos + namTab->currPos - 1 ] +=  -1.0;

          }

          break;


      }



      break;
    case 1:
      //irrelevant
      break;
  }

}

void buildC ( zeroCircuit *node, mnaSystem *system, hashTable *namTab, int row, int kPos )
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
            system->array_C[ row * posn1 + posn2  ] += tmp;
            system->array_C[ row * posn2 + posn1  ] += tmp;

          }

          if ( strcmp( node->linElement->connectors[0], "0\0" ) )
          {
            //posn1 = findPlace ( node->linElement->connectors[0], &nodes );
            posn1 = findHash ( namTab, node->linElement->connectors[0] );
            tmp = node->linElement->value;
            system->array_C[ row * posn1 + posn1  ] += tmp;
          }

          if ( strcmp( node->linElement->connectors[1], "0\0" ) )
          {
            //posn2 = findPlace ( node->linElement->connectors[1], &nodes );
            posn2 = findHash ( namTab, node->linElement->connectors[1] );
            tmp = node->linElement->value;
            system->array_C[ row * posn2 + posn2  ] += tmp;
          }

          break;

        case 5:

          tmp = -node->linElement->value;
          system->array_C[ row * (kPos + namTab->currPos - 1) + (kPos + namTab->currPos - 1)  ] += tmp;
          break;
      }
      break;

    case 1:
      //irrelevant
      break;
  }

}


void formatRowVec ( mnaSystem *finalSystem )
{
  int i;

  for ( i = 0; i < finalSystem->dim; i++ )
  {
    finalSystem->rowVector[i] = &finalSystem->array_A[ i*finalSystem->dim ];
    finalSystem->rowVectorC[i] = &finalSystem->array_C[i*finalSystem->dim ];
  }
}




mnaSystem *formatMnaTable ( zeroCircuit *circuit, hashTable *nameTab, int v_num, int l_num )
{
  int rows;
  int k_pos;

  zeroCircuit *cur;
  mnaSystem *finalSystem;

  //initialize k position
  k_pos = 0;

  finalSystem = ( mnaSystem * )malloc ( sizeof( mnaSystem ) );

  rows = ( nameTab->currPos ) + v_num + l_num ;

  printf("Megethos: %d\n", rows);

  finalSystem->array_A = (double *)calloc(rows * rows, sizeof(double));  //cols = rows
  finalSystem->array_C = (double *)calloc(rows * rows, sizeof(double));
  finalSystem->vector_B = (double *)calloc(rows, sizeof(double));
  finalSystem->vector_X = (char **)malloc( rows * sizeof (char *) );

  finalSystem->rowVector = (double **)calloc(rows, sizeof(double *));
  finalSystem->rowVectorC = (double **)calloc(rows, sizeof(double *));

  finalSystem->dim = rows;

  // initVectorX ( finalSystem, nameTab );

  //for ( cur = circuit->next; (cur->linElement == NULL && cur->nonlinElement == NULL); cur = cur->next )
  for( cur = circuit->next; ( (cur->linElement != NULL) || (cur->nonlinElement != NULL) ); cur = cur->next )
  {
    if (( cur->linElement->element == 2 ) || ( cur->linElement->element == 5 ))  //if element is voltage source or inductor
    {
      k_pos++;
      fillVectorX ( cur, finalSystem, k_pos, nameTab->currPos );
    }

    addElementStamp ( cur, finalSystem, nameTab, rows, k_pos );

    //Transient Analysis : Build C
    if (cmdList->_transient == 1)
    {
      buildC( cur, finalSystem, nameTab, rows, k_pos);
    }

  }

  formatRowVec ( finalSystem );

  return ( finalSystem );
}



void printMNA (hashTable *names, mnaSystem *systemMna )
{
  int i,j;

  printf( "\tINFO : Group 2 elements:\n\t\t-Voltage Sources: %d\n\t\t-Inductors: %d\n", vnumber, lnumber );
  printf( "\tINFO : Group 1 elements:\n\t\t-Resistors: %d\n\t\t-Current Sources: %d\n", rnumber, inumber );
  printf( "\tINFO : Number of nodes of sub-Matrix G: %d\n", names->currPos );



  FILE *bfile;

  bfile = fopen ( "mna_system.txt", "w" );


  fprintf( bfile, "\tINFO : Group 2 elements:\n\t\t-Voltage Sources: %d\n\t\t-Inductors: %d\n", vnumber, lnumber );
  fprintf( bfile, "\tINFO : Group 1 elements:\n\t\t-Resistors: %d\n\t\t-Current Sources: %d\n", rnumber, inumber );
  fprintf( bfile, "\tINFO : Number of nodes of sub-Matrix G: %d\n", names->currPos );
  fprintf( bfile, "\n\n\n\n" );
  fprintf ( bfile, "\t   -----B Vector----- \n\n" );
  for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
  {
    fprintf ( bfile, "\t\t%f ", systemMna->vector_B[i] );
    fprintf ( bfile, "\n");
  }

  fprintf( bfile, "\n\n\n\n" );
  /*
     printf ( "\t   -----X Vector----- \n\n" );
     for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
     {
     fprintf (bfile, "\t\t%s ", systemMna->vector_X[i] );
     fprintf(bfile, "\n");
     }
   */
  fprintf ( bfile, "MNA Array ... \n " );
  for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
  {

    fprintf(bfile, "\n");

    for( j=0; j<( names->currPos + vnumber + lnumber ); j++)
    {
      fprintf ( bfile, "%f\t", systemMna->array_A[i * (names->currPos + vnumber + lnumber) + j] );
    }
  }

  fprintf(bfile, "\n\n\n\n");
  if (cmdList->_transient)
  {
    fprintf ( bfile, "C Array ... \n " );
    for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
    {
      fprintf(bfile, "\n");
      for( j=0; j<( names->currPos + vnumber + lnumber ); j++)
      {
        fprintf ( bfile, "%f\t", systemMna->array_C[i * (names->currPos + vnumber + lnumber) + j] );
      }
    }
  }

  fprintf(bfile, "\n\n\n\n");
  fprintf ( bfile, "\t   -----Rows Vector----- \n\n" );
  for ( i=0; i<( names->currPos + vnumber + lnumber ); i++ )
  {
    fprintf ( bfile, "\t\t%lf ", *systemMna->rowVector[i] );
    fprintf(bfile, "\n");
  }

  fclose ( bfile );

}
