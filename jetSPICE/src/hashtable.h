/*
 * hashtable.h
 *
 *  Created on: Jan 2, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 *
 *      Description: Hashtable Implementation.
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRSIZE 40
#define SIZE 2048


/********************************************************************************
 *                                                                              *
 *                       DATA STRUCTURES DEFINITIONS                            *
 *                                                                              *
 ********************************************************************************/


struct _hashNode
{
	char *key;
	int pos;
	struct _hashNode *nxt;
};

struct _hashTable
{
	int size;
	int currPos;
	struct _hashNode **nodes;
};

typedef struct _hashNode hashNode;
typedef struct _hashTable hashTable;


/********************************************************************************
 *                                                                              *
 *                            FUNCTION DECLARATIONS                             *
 *                                                                              *
 ********************************************************************************/


/**
 *  Function that creates a hashtable.
 *  @param size size of the hashtable.
 *  @return the hashtable.
 */
hashTable *createHash( int size );


/**
 *  Function to destroy the hashtable.
 *  @param hashtbl the hashtable.
 */
void destroyHash (hashTable *hashtbl);


/**
 *  Function to insert an element in hashtable.
 *  @param hashtable the hashtable.
 *  @param key the node we want to insert.
 *  @result 0 if successful and -1 in case of failure.
 */
int insertHash (hashTable *hashtable, char *key);


/**
 *  Function to find an element and return the position.
 *  @param hastbl the hashtabl.
 *  @param key the name we search for.
 *  @return the position.
 */
int findHash (hashTable *hashtbl, char *key);


/**
 *  Hash Function
 *  @param s string to be hashed.
 *  @return the position in the hash table.
 */
int hashFunc(char *s);


/**
 *  Function to print all the elements of the hashtable.
 *  @param hashtbl the hashtable.
 */
void printHash (hashTable *hashtbl);

#endif /* HASHTABLE_H_ */
