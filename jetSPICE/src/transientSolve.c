/*
 * transientSolve.c
 *
 *  Created on: Feb 19, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */
#include "transientSolve.h"

extern commandList *cmdList;


void backwardEuler(mnaSystem *system, zeroCircuit *circuit, hashTable *hashTab)
{

        double h, beta;
        double *xold, *xnew, *c, *e, *tmp, *b;
        double t, tf;
        int iters = 0, k = 0, size, status;
        double *M, **rowVecM;

        //For PULSE and PWL
        int periods = 0, s[3] = {0,0,0};

        size = system->dim;


        if( !( M = (double *)calloc(size * size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( xold = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( xnew = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( c = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( e = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( tmp = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( b = (double *)calloc(size, sizeof(double)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        if( !( rowVecM = (double **)calloc( size, sizeof(double *)) ) )
        {
            printf("Not enough memory...\n");
            exit(-1);
        }
        ;

        vecCpy(xold, system->dcPointRes, size);
        vecCpy(e, system->vector_B, size);

        tf = cmdList->_tranFinTime;
        h = cmdList->_tranTimeStep;
        beta = 1.0 / h;
        addMMScal(M, 1.0, system->array_A, beta, system->array_C, size);
        formatRowVecM(rowVecM, M, size);

        for (t = 0.0; t <= tf; t+=h )
        {
            k++;
            clearVec(tmp, size);
            mvProd(tmp, system->rowVectorC, xold, size, 0);
            vecMult(c, tmp, -beta, size);
            constructE(e, t, s, &periods, hashTab, circuit );
            vecSumScal(b, e, c, 1.0, size);
            if ( (cmdList->_options == ITER ) || (cmdList->_options == ITERSPD) )
            {
                b = bicg(rowVecM, b, size, &iters, cmdList->_itolValue);
                printVec(b, size);
            }
            else
            {
                int *pivotTab;
                double **newPiv;

                pivotTab = ( int * )calloc( size , sizeof ( int ) );
                newPiv = luDecomp ( rowVecM, size, pivotTab );
                b = luBackSubst( newPiv, size, pivotTab, b);

            }
            vecCpy(xnew, b, size);
            vecCpy(xold, xnew, size);
        }

        printVec(xold, size);

        free(xold);
        free(e);
        free(c);
        free(tmp);
        free(b);
}

/*
void trapezoidal(mnaSystem *system, zeroCircuit *circuit, hashTable *hashTab)
{

}
*/

void cs_backwardEuler(mnaSpSystem *system, zeroCircuit *circuit, hashTable *hashTab)
{

    double h, beta;
    double *xold, *xnew, *c, *e, *tmp, *b;
    double t, tf;
    int iters = 0, k = 0, size, status, i;
    cs *M;

    FILE * tr;

    tr = fopen ( "transient_vector.txt", "w" );
    fprintf ( tr, "Result vector\n Transient Analysis\n\n" );;

    //For PULSE and PWL
    int periods = 0, s[3] = {0,0,0};

    size = system->dim;
    M = cs_spalloc(size,size,size * size,1,1);

    if( !( xold = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( xnew = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( c = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( e = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( tmp = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( b = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }

    vecCpy(xold, system->dcPointRes, size);
    vecCpy(e, system->vector_B, size);

    tf = cmdList->_tranFinTime;
    h = cmdList->_tranTimeStep;
    fprintf ( tr, "step = %f \t finTime = %f \n\n",h, tf );
    beta = 1.0 / h;
    M = cs_add(system->array_A, system->array_C, 1.0, beta);
    for (t = 0.0; t <= tf; t+=h )
    {
        k++;
        clearVec(tmp, size);
        status = cs_gaxpy(system->array_C, xold, tmp);
        vecMult(c, tmp, -beta, size);
        constructE(e, t, s, &periods, hashTab, circuit );
        vecSumScal(b, e, c, 1.0, size);
        if ((cmdList->_options == ITER ) || (cmdList->_options == ITERSPD))
        {
            b = cs_bicg (M, b, size, &iters, cmdList->_itolValue );
        }
        else
        {
            status = cs_lusol (2, M, b, 1);
        }
        vecCpy(xnew, b, size);
        vecCpy(xold, xnew, size);

        //Print Outputs
        fprintf ( tr, "\n\n t = %f \t \n\n",t );
        for ( i = 0; i < size; i++ )
        {
                fprintf ( tr, "%.10f\n", xold[i] );
        }
    }

    printVec(xold, size);

    free(xold);
    free(e);
    free(c);
    free(tmp);
    free(b);
    cs_spfree(M);
    fclose ( tr );
}


void cs_trapezoidal(mnaSpSystem *system, zeroCircuit *circuit, hashTable *hashTab)
{
    double h, beta;
    double *xold, *e, *eold, *enew, *tmp, *b;
    double t, tf;
    int temp, iters = 0, k = 0, size, status, i;
    cs *M1, *M2;

    FILE * tr;

    tr = fopen ( "transient_vector.txt", "w" );
    fprintf ( tr, "Result vector\n Transient Analysis\n\n" );;

    //For PULSE and PWL
    int periods = 0, s[3] = {0,0,0};

    size = system->dim;

    M1 = cs_spalloc(size,size,size * size,1,1);
    M2 = cs_spalloc(size,size,size * size,1,1);

    if( !( xold = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( e = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( eold = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( enew = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( b = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }
    if( !( tmp = (double *)calloc(size, sizeof(double)) ) )
    {
        printf("Not enough memory...\n");
        exit(-1);
    }

    vecCpy(xold, system->dcPointRes, size);
    vecCpy(eold , system->vector_B, size);
    vecCpy(enew, system->vector_B, size);

    tf = cmdList->_tranFinTime;
    h = cmdList->_tranTimeStep;
    fprintf ( tr, "step = %f \t finTime = %f \n\n",h, tf );
    beta = 2.0 / h;
    M1 = cs_add(system->array_A, system->array_C, 1.0, beta);
    M2 = cs_add(system->array_A, system->array_C, 1.0, -beta);

    for (t = 0.0; t <= tf; t+=h )
    {
        k++;
        clearVec(tmp, size);
        temp = cs_gaxpy(M2, xold, tmp);
        constructE(enew, t, s, &periods, hashTab, circuit );
        vecSumScal(e, enew, eold, 1.0, size);
        vecSumScal(b, e, tmp, -1.0, size);
        if ((cmdList->_options == ITER ) || (cmdList->_options == ITERSPD))
        {
            b = cs_bicg (M1, b, size, &iters, cmdList->_itolValue );
        }
        else
        {
            status = cs_lusol (2, M1, b, 1);
        }
        vecCpy(xold, b, size);
        vecCpy(eold, enew, size);

        //Print Outputs
        fprintf ( tr, "\n\n t = %f \t \n\n",t );
        for ( i = 0; i < size; i++ )
        {
                fprintf ( tr, "%.10f\n", xold[i] );
        }
    }

    printVec(xold, system->dim);

    free(xold);
    free(e);
    free(eold);
    free(enew);
    free(tmp);
    //free(b);
    cs_spfree(M1);
    cs_spfree(M2);
    fclose(tr);

}


void constructE(double *e, double t, int *s, int *periods, hashTable *hashTab, zeroCircuit *circuit )
{

    double tmp1, tmp2, t1, t2, step1, step2;
    int v1, v2, size, i;
    double res;
    double *tmp;
    zeroCircuit *cur;

    for( cur = circuit->next; ( (cur->linElement != NULL) || (cur->nonlinElement != NULL) ); cur = cur->next )
    {
        if ((cur->linElement->element == V) || (cur->linElement->element == I))
        {
        	if (cur->linElement->vInfo != NULL)
        	{
        		switch (cur->linElement->vInfo->_infoType)
                {
                	case EXP:
                	{
                		tmp = cur->linElement->vInfo->_data;
                		if(t < tmp[2])
                		{
                			res = tmp[0];
                		}
                		else if((t >= tmp[2]) && (t < tmp[4]))
                		{
                			tmp1 = -(t - tmp[2])/tmp[3];
                			res = tmp[0] + (tmp[1] - tmp[0]) * (1.0 - exp(tmp1));
                		}
                		else
                		{
                			tmp1 = -(t - tmp[4])/tmp[5];
                			tmp2 = -(t - tmp[2])/tmp[3];
                			res = tmp[0] + (tmp[1] - tmp[0]) * (exp(tmp1) - exp(tmp2));
                		}
                		break;
                	}
                	case SIN:
                	{
                		tmp = cur->linElement->vInfo->_data;
                    	if(t < tmp[3])
                    	{
                        	res = tmp[0] + tmp[1] * sin(2 * M_PI * tmp[5] / 360);
                    	}
                    	else
                    	{
                        	tmp1 = -(t - tmp[3]) * tmp[4];
                        	res = tmp[0] + tmp[1] * sin((2 * M_PI * tmp[2] * (t - tmp[3])) + (2 * M_PI * tmp[5] / 360)) * exp(tmp1);
                    	}
                		break;
                	}
                	case PULSE:
                	{
                		tmp = cur->linElement->vInfo->_data;
                		if (t > (tmp[2] + tmp[6] + (*periods) * tmp[6] ))
                		{
                			(*periods)++;
                		}

                		tmp = cur->linElement->vInfo->_data;
                		t1 = fabs(tmp[3]/cmdList->_tranTimeStep);
                		t2 = fabs(tmp[4]/cmdList->_tranTimeStep);
                		step1 = (tmp[1] - tmp[0])/t1;
                		step2 = (tmp[0] - tmp[1])/t2;
                		if(t < tmp[2])
                		{
                			s[0] = 0;
                			res = tmp[0];
                		}
                		else if( (t >= (tmp[2] + *periods * tmp[6] ) ) && (t < (tmp[2] + tmp[3] + *periods * tmp[6]) ) )
                		{
                			s[0] +=1;
                			res = tmp[0] + s[0] * step1;

                		}
                		else if((t >= (tmp[2] + tmp[3] + *periods * tmp[6] ))
                				&& ((t < ( tmp[2] + tmp[3] + tmp[5] + *periods * tmp[6] ))))
                		{
                			s[0] = 0;
                			res = tmp[1];
                		}
                		else if ( (t >= (tmp[2] + tmp[3] + tmp[5] + *periods * tmp[6] ))
                				&& ((t < ( tmp[2] + tmp[3] + tmp[5] + tmp[4] + *periods * tmp[6] ))))
                		{
                			s[0] +=1 ;
                			res = tmp[1] + s[0] * step2 ;
                		}
                		else
                		{
                			s[0] = 0;
                			res = tmp[0];
                		}

                		break;
                	}
                	case PWL:
                	{
                		size = cur->linElement->vInfo->_nofData;
                		tmp = cur->linElement->vInfo->_data;
                		i = s[2];
                   		if((i < size-2) && t < tmp[i+2] )
                		{
                   			s[1]+=1;
              				t1 = (tmp[i+2] - tmp[i]) / cmdList->_tranTimeStep;
               				step1 = (tmp[i+3]-tmp[i+1])/t1;
               				res = tmp[i] + s[1] * step1;
               				break;
               			}
                   		else if ((i < size-2))
                   		{
                   			i+=2;
                   			s[2]+=2; s[1]=1;
                   			t1 = (tmp[i+2] - tmp[i]) / cmdList->_tranTimeStep;
                   			step1 = (tmp[i+3]-tmp[i+1])/t1;
                   			res = tmp[i] + s[1] * step1;
                   			break;
                   		}
                   		else
                   		{
                   			printf("Error ... \n");
                   		}
                   		break;
                	}
                	default:
                	{
                		printf("No such type ... \n");
                		break;
                	}
                }

        		switch(cur->linElement->element)
        		{
                	case V :
                	{
               			e[hashTab->currPos- 1 + cur->linElement->k] -= cur->linElement->vInfo->pVal;
               			e[hashTab->currPos - 1 + cur->linElement->k] += res;
               			cur->linElement->vInfo->pVal = res;
               			break;
                	}
                	case I :
                	{
               			if ( strcmp( cur->linElement->connectors[0], "0\0" ))
               			{
               				v1 = findHash ( hashTab, cur->linElement->connectors[0] );
               				e[v1] -= cur->linElement->vInfo->pVal;
               				e[v1] += res;
               			}
               			if ( strcmp( cur->linElement->connectors[1], "0\0" ))
               			{
               				v2 = findHash ( hashTab, cur->linElement->connectors[1] );
               				e[v2] -= cur->linElement->vInfo->pVal;
               				e[v2] += res;
               			}
               			cur->linElement->vInfo->pVal = res;
               			break;
                	}
                	default:
                	{
                		printf("No such type ... \n");
                		break;
                	}
        		}
        	}
        }
    }
}

void formatRowVecM ( double **rV, double *A, int size )
{
        int i;

        for ( i = 0; i < size; i++ )
        {
                rV[i] = &A[ i * size ];
        }
}
