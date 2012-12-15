
/*
 * hashtable.c
 *
 * Created on: Jan 2, 2012
 * Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */
#include "hashtable.h"



hashTable *createHash( int size )
{
  hashTable *hashtbl;

  if( !( hashtbl = malloc( sizeof ( hashTable ) ) ) )
  {
    return NULL;
  }
  if ( !( hashtbl->nodes = calloc( size , sizeof( hashNode ) ) ) ) {
    free(hashtbl);
    return NULL;
  }
  hashtbl->size = size;
  hashtbl->currPos = 0;

  return (hashtbl);

}


void destroyHash (hashTable *hashtbl)
{
  int n;
  hashNode *node, *oldnode;

  for(n = 0; n < hashtbl->size; n++)
  {
    node = hashtbl->nodes[n];
    while ( node )
    {
      free(node->key);
      oldnode = node;
      node = node->nxt;
      free(oldnode);
    }
  }
  free(hashtbl->nodes);
  free(hashtbl);
}


int insertHash (hashTable *hashtbl, char *key)
{
  hashNode *node;
  int hash;

  hash = hashFunc(key) % hashtbl->size;
  node = hashtbl->nodes[hash];

  if ( strcmp(key, "0") == 0 )
    return 0;

  while (node)
  {
    if(!strcmp(node->key, key))
    {
      //printf("Node Already Exists ... \n");
      return 0;
    }
    node = node->nxt;
  }

  if( !(node = malloc( sizeof ( hashNode ) ) ))
  {
    return -1;
  }

  if (!(node->key = strdup(key)))
  {
    free(node);
    printf("Can't save node name ... \n");
    return -1;
  }
  node->pos = hashtbl->currPos;
  hashtbl->currPos++;

  node->nxt = hashtbl->nodes[hash];
  hashtbl->nodes[hash] = node;

  return 0;

}


int findHash (hashTable *hashtbl, char *key)
{
  hashNode *node;
  int hash;

  hash = hashFunc(key) % hashtbl->size;
  node = hashtbl->nodes[hash];

  while(node) {
    if(!strcmp(node->key, key))
    {
      return node->pos;
    }
    node = node->nxt;
  }
  return -1;
}


int hashFunc(char *s)
{
  unsigned int tmp = 0;

  while (*s)
  {
    tmp = tmp * 37 + *s++;
  }

  return tmp%SIZE;
}


void printHash (hashTable *hashtbl)
{
  int n;
  hashNode *node;

  for(n=0; n < hashtbl->size; n++)
  {
    printf("HashNode : %d \n",n);
    node = hashtbl->nodes[n];
    while ( node )
    {
      printf(" Node : %s --- Pos : %d \n", node->key, node->pos);
      node = node->nxt;
    }
  }
}


