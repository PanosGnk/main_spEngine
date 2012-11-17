/*
 * linearSolve.c
 *
 *  Created on: Jan 9, 2012
 *      Author: Giannakou Panagiotis Taxiarchis, Pinitas Georgios, Gourgounia Sofia, Skourths Anastasios
 */

#include <math.h>
#include <string.h>
#include "linearSolve.h"


/********************************************************************************
 *                                                                              *
 *                            FUNCTION DEFINITIONS                              *
 *                                                                              *
 ********************************************************************************/

void choleskyDecomp ( double **rowVector, int size, double *diag )
{

	int i,j,k;
	double sum;


	for ( i = 0; i < size; i++ )
	{
		for ( j = i; j < size; j++ )
		{
			for ( sum = *(rowVector[i] + j), k = i - 1; k >= 0; k-- )
			{
				sum -=  *(rowVector[i] + k ) *  ( *(rowVector[j] + k ) );
			}

			if ( i == j )
			{
				if ( sum <= 0.0 )
				{
					printf ( "\n Cholesky Decomposition failed. Exiting .....\n" );
					exit ( 0 );
				}

				diag[i] = sqrt ( sum );
			}
			else
			{
				 *(rowVector[j] + i ) = sum / diag[i];
			}

		}
	}
}

double *choleskySolve ( double **rowVector, int size, double *diag, double *bVector )
{
	int i,k;
	double sum;

	double *sol;

	sol = (double *)calloc(size, sizeof(double));


	for ( i = 0; i < size; i++ )
	{
		for ( sum = bVector[i], k = i-1; k >= 0; k-- )
		{
			sum -= *(rowVector[i] + k) * sol[k];
		}

		sol[i] = sum / diag[i];
	}

	for ( i = size - 1; i >= 0; i-- )
	{
		for ( sum = sol[i], k = i + 1; k < size; k++ )
		{
			sum -= (*(rowVector[k] + i)) * sol[k];
		}

		sol[i] = sum / diag[i];
	}

	return ( sol );
}

double **luDecomp ( double **rowVector, int size, int *pivotTable )
{
	int i,j,k;

	int m, temp_piv;

	double xVal;
	double *temp;

	double **pivots;

	pivots = ( double ** )malloc( size * sizeof ( double * ) );
	memcpy( pivots, rowVector, size * sizeof ( double * ) );


	for ( i = 0; i < size; i++)
	{
		pivotTable[i] = i;
	}


	for ( i = 0; i < size; i++ )
	{
		//partial pivoting
		xVal = fabs ( *(pivots[i] + i) );

		//m=0;
		for( j = i; j < size; j++ )
		{
			if ( fabs ( *(pivots[j] + i) ) > xVal )
			{
				m = j;


				//update row table
				temp = pivots[m];
				pivots[m] = pivots[i];
				pivots[i] = temp;

				//update pivot table
				temp_piv = pivotTable[m];
				pivotTable[m] = pivotTable[i];
				pivotTable[i] = temp_piv;
			}
		}




		//factorizations
		k=0;
		j=0;

		for ( k = i + 1; k < size; k++ )
		{
			if ( *(pivots[i] + i) == 0 )
			{
				*(pivots[k] + i) = 0;
			}
			else
			{
				*(pivots[k] + i) = *(pivots[k] + i) /  *(pivots[i] + i); //= l[i][k]
			}
			for ( j = i + 1; j < size; j++ )
			{
				*(pivots[k] + j) -= ( *(pivots[k] + i) ) * ( *(pivots[i] + j) );
			}
		}
	}

	return ( pivots );
}

double *luBackSubst ( double **rowVector, int size, int *pivotTable, double *bVector )
{
	int i;
	int ii = -1;
	int ip, j;

	double sum = 0.0;
	double *bcopy;

	bcopy = ( double * )malloc( size * sizeof ( double ) );

	for ( i=0; i < size; i++ )
	{
		bcopy[i] = bVector[ pivotTable [i] ];
	}

	for ( i = 0; i < size; i++ )
	{
		for ( j = 0; j < i; j++  )
		{
			bcopy[i] -= ( *(rowVector[i] + j) ) * bcopy[j];
		}
	}

	for ( i = size -1; i >= 0; i-- )
	{
		sum = bcopy[i];
		for (j = i + 1; j < size; j++ )
		{
				sum -= *(rowVector[i] + j) * bcopy[j];
		}


		if ( *(rowVector[i] + i) == 0 )
		{
			bcopy[i] = sum / *(rowVector[i] + i);

		}
		else
		{
			bcopy[i] = sum / *(rowVector[i] + i);
		}
	}

	return ( bcopy );
}

void jacobiPrecond (double *M, double **A, int size)
{
  int i;

  for(i = 0; i < size; i++)
    {
      if ( (*(A[i] + i)) == 0.0 )
        {
          M[i] = 1.0;
        }
      else
        {
          M[i] = (*(A[i] + i)) ;
        }
    }
}

void precSolve (double *M, double *r, double *z, int size)
{
  int i;

  for (i=0; i<size; i++)
    {
      z[i] = r[i] * M[i];
    }
}

double *cg (double **rowVector, double *bVector, int size, int *iter, double itol )
{

  double alpha, beta, omega;
  double bnorm2, rnorm2, itol2;
  double rho,rho1;
  double *p, *q, *r, *z, *x, *M;

  p = (double *)calloc(size, sizeof(double));
  q = (double *)calloc(size, sizeof(double));
  r = (double *)calloc(size, sizeof(double));
  z = (double *)calloc(size, sizeof(double));
  x = (double *)calloc(size, sizeof(double));
  M = (double *)calloc(size, sizeof(double));

  *iter = 0;
  itol2 = itol * itol;
  bnorm2 = innerProd(bVector, bVector, size);
  if (bnorm2==0.0)
    {
      bnorm2=1.0;
    }

  vecCpy(r,bVector,size);
  rnorm2 = innerProd(r , r, size);

  jacobiPrecond(M, rowVector, size);
  invVec(M, size);

  while ( (rnorm2 > (bnorm2 * itol2)) &&  (*iter < size) )
    {
      *iter+=1;

      precSolve (M, r, z, size);

      rho = innerProd(r, z, size);

      if ( *iter == 1 )
        {
          vecCpy(p,z,size);
        }
      else
        {
          beta = rho / rho1;
          vecSumScal(p,z,p,beta,size);
        }

      rho1=rho;

      mvProd(q, rowVector, p, size, 0);

      omega = innerProd(p, q, size);
      alpha=rho / omega;

      vecSumScal(x,x,p,alpha,size);
      vecSumScal(r,r,q,-alpha,size);

      rnorm2=innerProd(r,r,size);
    }

  free(p);
  free(q);
  free(r);
  free(z);
  free(M);

  return( x );
}

double *bicg (double **rowVector, double *bVector, int size, int *iter, double itol )
{

  double alpha, beta, omega;
  double bnorm2, rnorm2, itol2;
  double rho,rho1;
  double *p, *pp, *q, *qq, *r, *rr, *z, *zz, *x, *M;

  p = (double *)calloc(size, sizeof(double));
  pp = (double *)calloc(size, sizeof(double));
  q = (double *)calloc(size, sizeof(double));
  qq = (double *)calloc(size, sizeof(double));
  r = (double *)calloc(size, sizeof(double));
  rr = (double *)calloc(size, sizeof(double));
  z = (double *)calloc(size, sizeof(double));
  zz = (double *)calloc(size, sizeof(double));
  x = (double *)calloc(size, sizeof(double));
  M = (double *)calloc(size, sizeof(double));

  *iter = 0;
  itol2 = itol * itol;
  bnorm2 = innerProd(bVector, bVector, size);

  if (bnorm2==0.0)
    {
      bnorm2=1.0;
    }

  vecCpy(r,bVector,size);
  vecCpy(rr,bVector,size);

  rnorm2 = innerProd(r , r, size);

  jacobiPrecond(M, rowVector, size);
  invVec(M, size);

  while ( (rnorm2 > (bnorm2 * itol2)) &&  (*iter < size) )
    {
      *iter+=1;

      precSolve (M, r, z, size);
      precSolve (M, rr, zz, size);

      rho = innerProd(rr, z, size);

      if (rho < EPS)
        {
          printf("Bi-Conjugate Gradient Algorithm failed ... \n");
          exit(0);
        }

      if ( *iter == 1 )
        {
          vecCpy(p,z,size);
          vecCpy(pp,zz,size);
        }
      else
        {
          beta = rho / rho1;
          vecSumScal(p,z,p,beta,size);
          vecSumScal(pp,zz,pp,beta,size);
        }

      rho1=rho;

      mvProd(q, rowVector, p, size, 0);
      mvProd(qq, rowVector, pp, size, 1);

      omega = innerProd(pp, q, size);
      if (omega < EPS)
        {
          printf("Bi-Conjugate Gradient Algorithm failed ... \n");
          exit(0);
        }

      alpha = rho / omega;

      vecSumScal(x,x,p,alpha,size);
      vecSumScal(r,r,q,-alpha,size);
      vecSumScal(rr,rr,qq,-alpha,size);

      rnorm2=innerProd(r,r,size);
    }

  free(p);
  free(pp);
  free(q);
  free(qq);
  free(r);
  free(rr);
  free(z);
  free(zz);
  free(M);

  return( x );
}

void update(double *x, double *u, double **A, double *b,double *h, double *s, int resd,int hsize, int size,int iter)
{
  int i,j;

  double *tmp, *ss, *y;
  double sum;

  ss = (double *)calloc(resd, sizeof(double));
  y = (double *)calloc(resd, sizeof(double));
  tmp = (double *)calloc(size, sizeof(double));

  vecCpy(ss,s,resd);

  for (i = 0; i < resd; i++)
    {
      if (h[i * hsize + i] == 0)
        {
          printf("Matrix Is Singular ... \n");
          exit(0);
        }
    }

  //Solve Hy=ss;
  for (i = resd;  i >= 0; --i)
    {
      sum = 0.0;
      for (j = i+1; i <= resd; ++j)
        {
          sum += h[i * hsize + j] * y[j];
        }
      y[i] = (s[i] - sum) / h[i * hsize + i];
    }

  for (i = 0; i < resd; i++)
    {
      vecMult(tmp,&u[i * size],y[i],size);

      for (j = 0; j < size; j++)
        {
          x[i] += tmp[i];
        }
    }

  for (j = 0; j < size; j++)
    {
      tmp[i] = 0.0;
    }

  mvProd(tmp,A,x,size,0);
  vecDiff(tmp,b,tmp,size);

  s[resd + 1] = norm2(tmp,size);


}

double *gmres(double **rowVector, double *bVector, int size, int *iter,int res, double itol)
{

  int i,k,hsize;

  double bnorm, rnorm;
  double est;
  double *r, *u, *s, *w, *x, *M, *ax, *e, *h;

  r = (double *)calloc(size, sizeof(double));
  w = (double *)calloc(size, sizeof(double));
  x = (double *)calloc(size, sizeof(double));
  M = (double *)calloc(size, sizeof(double));
  ax = (double *)calloc(size, sizeof(double));
  e = (double *)calloc(size, sizeof(double));

  *iter = 0;

  hsize = res + 1;

  u = (double *)calloc(hsize * size, sizeof(double));
  h = (double *)calloc(hsize * hsize, sizeof(double ));
  s = (double *)calloc(hsize, sizeof(double));

  bnorm = norm2(bVector, size);

  e[0]=1.0;
  jacobiPrecond(M, rowVector, size);
  invVec(M, size);

  while ( (est > EPS) &&  (*iter < size) )
    {
      *iter+=1;

      mvProd(ax,rowVector,x,size,0);
      vecDiff(r,bVector,ax,size);
      precSolve (M, r, r, size);

      rnorm = norm2(r, size);

      vecDivScal(&u[0],r,rnorm,size);
      vecMult(s,e,rnorm,size);

      for (i = 0; i < res; i++)
        {

          mvProd(w,rowVector,&u[i * size],size,0);
          precSolve(M,w,w,size);

          for (k = 0; k < i; k++)
            {
              h[k * hsize + i] = innerProd(w,&u[k * size],size);
              vecDiffScal(w,w,&u[k * size],h[k * hsize + i],size);
            }

          h[(i + 1) * hsize + i] = norm2(w,size);
          vecDivScal(&u[(i + 1) * size],w,h[(i + 1) * hsize + i],size);




        }

      //est = ;
      est /= bnorm;
    }

  free(r);
  free(M);

  return( x );
}

void cs_jacobiPrecond (double *M, cs *A, int size)
{
  int j, p, flag;

  for (j = 0; j < size; j++)
    {
	  flag = 0;
	  for (p = A->p[j]; p < A->p[j+1]; p++)
	  {
		  if ( j == A->i[p])
		  {
			flag = 1;
			M[j] = A->x[p];
		  }
	  }

	  if ( flag == 0 )
	  {
		  M[j] = 1.0;
	  }

    }
}

double *cs_cg (cs *A, double *bVector, int size, int *iter, double itol )
{

  int state;
  double alpha, beta, omega;
  double bnorm2, rnorm2, itol2;
  double rho,rho1;
  double *p, *q, *r, *z, *x, *M;

  p = (double *)calloc(size, sizeof(double));
  q = (double *)calloc(size, sizeof(double));
  r = (double *)calloc(size, sizeof(double));
  z = (double *)calloc(size, sizeof(double));
  x = (double *)calloc(size, sizeof(double));
  M = (double *)calloc(size, sizeof(double));

  *iter = 0;
  itol2 = itol * itol;
  bnorm2 = innerProd(bVector, bVector, size);
  if (bnorm2==0.0)
    {
      bnorm2=1.0;
    }

  vecCpy(r,bVector,size);
  rnorm2 = innerProd(r , r, size);

  cs_jacobiPrecond(M, A, size);
  invVec(M, size);

  while ( (rnorm2 > (bnorm2 * itol2)) &&  (*iter < size) )
    {
      *iter+=1;

      precSolve (M, r, z, size);

      rho = innerProd(r, z, size);

      if ( *iter == 1 )
        {
          vecCpy(p,z,size);
        }
      else
        {
          beta = rho / rho1;
          vecSumScal(p,z,p,beta,size);
        }

      rho1=rho;

      clearVec(q, size);
      state = cs_gaxpy(A,p,q);
      if (state == 0)
      {
    	  printf("ERROR : in CG q[i] = A * p[i] multiplication using Sparse Matrices... \n");
    	  exit(0);
      }

      omega = innerProd(p, q, size);
      alpha=rho / omega;

      vecSumScal(x,x,p,alpha,size);
      vecSumScal(r,r,q,-alpha,size);

      rnorm2=innerProd(r,r,size);
    }

  free(p);
  free(q);
  free(r);
  free(z);
  free(M);

  return( x );
}

double *cs_bicg (cs *A, double *bVector, int size, int *iter, double itol )
{

	int state;
	  double alpha, beta, omega;
	  double bnorm2, rnorm2, itol2;
	  double rho,rho1;
	  double *p, *pp, *q, *qq, *r, *rr, *z, *zz, *x, *M;
	  cs *AT;

	  p = (double *)calloc(size, sizeof(double));
	  pp = (double *)calloc(size, sizeof(double));
	  q = (double *)calloc(size, sizeof(double));
	  qq = (double *)calloc(size, sizeof(double));
	  r = (double *)calloc(size, sizeof(double));
	  rr = (double *)calloc(size, sizeof(double));
	  z = (double *)calloc(size, sizeof(double));
	  zz = (double *)calloc(size, sizeof(double));
	  x = (double *)calloc(size, sizeof(double));
	  M = (double *)calloc(size, sizeof(double));
	  AT = cs_transpose(A,1);

	  *iter = 0;
	  itol2 = itol * itol;
	  bnorm2 = innerProd(bVector, bVector, size);

	  if (bnorm2 == 0.0)
	    {
	      bnorm2 = 1.0;
	    }

	  vecCpy(r,bVector,size);
	  vecCpy(rr,bVector,size);

	  rnorm2 = innerProd(r , r, size);

	  cs_jacobiPrecond(M, A, size);
	  invVec(M, size);

	  while ( (rnorm2 > (bnorm2 * itol2)) &&  (*iter < size) )
	    {
	      *iter+=1;

	      precSolve (M, r, z, size);
	      precSolve (M, rr, zz, size);

	      rho = innerProd(rr, z, size);

	      if (fabs(rho) < EPS)
	        {
	          printf("Bi-Conjugate Gradient Algorithm failed ... \n");
	          exit(0);
	        }

	      if ( *iter == 1 )
	        {
	          vecCpy(p,z,size);
	          vecCpy(pp,zz,size);
	        }
	      else
	        {
	          beta = rho / rho1;
	          vecSumScal(p,z,p,beta,size);
	          vecSumScal(pp,zz,pp,beta,size);
	        }

	      rho1=rho;

	      clearVec(q, size);
	      state = cs_gaxpy(A,p,q);
	      if (state == 0)
	      {
	    	  printf("Error in Bi-CG q[i] = A * p[i] multiplication using Sparse Matrices... \n");
	          exit(0);
	      }

	      clearVec(qq, size);
	      state = cs_gaxpy(AT,pp,qq);
	      if (state == 0)
	      {
	      	  printf("Error in Bi-CG qq[i] = A * pp[i] multiplication using Sparse Matrices... \n");
	       	  exit(0);
	      }

	      omega = innerProd(pp, q, size);

	      alpha = rho / omega;

	      vecSumScal(x,x,p,alpha,size);
	      vecSumScal(r,r,q,-alpha,size);
	      vecSumScal(rr,rr,qq,-alpha,size);

	      rnorm2=innerProd(r,r,size);
	    }

	  free(p);
	  free(pp);
	  free(q);
	  free(qq);
	  free(r);
	  free(rr);
	  free(z);
	  free(zz);
	  free(M);
	  cs_spfree(AT);

	  return( x );
}

void solveSystem ()
{
		double *p, *result;
		double **piv;

		int size;

		int i,j;
		int *pivotTab ;

		int itern;
		double tol = 1e-4;
		size = 4;

		itern=0;

		double tab[16] = { 5, 1.2, 0.3, -0.6, 1.2, 6, -0.4, 0.9, 0.3, -0.4, 8, 1.7, -0.6, 0.9, 1.7, 10 };
		//double tab[9] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0 };

		double *rows [4] = { &tab[0], &tab[4], &tab[8], &tab[12] };
		//double *rows [3] = { &tab[0], &tab[3], &tab[6] };

		double bvec[4] = { 3.0, 5.0, 7.0, 9.0 };

		//p = ( double * )malloc ( size * sizeof( double ) );
		//pivotTab = ( int * )calloc( size , sizeof ( int ) );

		//choleskyDecomp ( rows, size, p );
		//result = choleskySolve ( rows, size, p, bvec );

		//piv = luDecomp ( rows, size, pivotTab );

		//luBackSubst( piv, size, pivotTab, bvec);

		result = bicg(rows, bvec, size, &itern, tol);

		printf("\n\n");
		printf("Iteration Number : %d \n",itern);
		printf("\n\n");
		printf ("\t--- Solution ---\n");

		for ( i=0; i < size; i++)
		{
			printf("\t    %lf\n", result[i]);
		}

		printf("\n\n");
		/*
		printf ("\t--- L ---\n");

			for ( i=0; i < size; i++)
			{
				printf("\n");
				for ( j=0; j < size; j++ )
				{
					printf("\t    %.10lf", *(rows[i] + j));
				}

			}
			*/
}

/********************************************************************************
 *                                                                              *
 *                            UTILITY FUNCTIONS                                 *
 *                                                                              *
 ********************************************************************************/


char **transposeVectorX ( int *pivots, char **vec, int size)
{
	char **temp;
	int i;

	temp = ( char ** )malloc( size * sizeof ( char * ) );

	for ( i=0; i < size; i++ )
	{
		temp[i] = strdup (vec[ pivots [i] ]);
	}

	return temp;
}

int findVarPos ( char *str, char **vector , int size)
{
	int pos;
	int i,j;
	int exist[2]={0,0};

	for ( i=0; i < size; i++ )
	{
		exist[0] = 0;
		exist[1] = 0;

		for ( j=0; j < 2; j++ )
		{
			if( str[j]!= vector[i][j+2])
			{
				exist[j]=0;
			}
			else
			{
				exist[j]=1;
			}
		}

		if ( exist[0]==1 && exist[1]==1 )
		{
			pos=i;
			return(pos);
		}

	}

	//return pos;
	return -1;
}

void clearVec (double *x, int size)
{
	int i;

	for (i = 0; i < size; i++)
	{
		x[i] = 0.0;
	}

}

double innerProd(double *x, double *y, int size)
{
  int i;
  double sum = 0.0;

  for( i = 0; i < size; i++)
    {
      sum += x[i] * y[i];
    }

  return sum;

}

double norm2(double *x, int size)
{
  int i;
  double sum = 0.0;

  for( i = 0; i < size; i++)
    {
      sum += x[i] * x[i];
    }

  return (sqrt(sum));

}

void vecCpy (double *x, double *y, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = y[i];
    }
}

void vecCpyScal(double *x, double *y, double s, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] += s * y[i];
    }
}

void vecSumScal(double *x, double *y, double *z, double s, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = y[i] + s * z[i];
    }
}

void vecDivScal(double *x, double *y, double s, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = y[i] / s ;
    }
}

void vecDiff(double *x, double *y, double *z, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = y[i] - z[i];
    }

}

void vecDiffScal(double *x, double *y, double *z, double s, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = y[i] - ( s * z[i] );
    }

}

void vecMult(double *x, double *y, double s, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      x[i] = s * y[i];
    }

}

void invVec (double *x, int size)
{
  int i;

  for (i = 0; i < size; i++)
    {
      if (x[i] == 0.0) {
          printf(" Division by zero \n ");
          exit ( 0 );
      }
      x[i] = 1.0 / x[i];
    }
}

void mvProd (double *Ax, double **A, double *x, int size, int method)
{
  int i,j;

  switch(method)
  {
  case 0:
    {
      for ( i = 0; i < size; i++ )
        {
          Ax[i]=0.0;
          for ( j = 0; j < size; j++)
            {
              Ax[i] += (*(A[i] + j)) * x[j];
            }
        }
      break;
    }
  case 1:
    {
      for ( j = 0; j < size; j++ )
        {
          Ax[j]=0.0;
          for ( i = 0; i < size; i++)
            {
              Ax[j] += (*(A[i] + j)) * x[i];
            }
        }
      break;
    }
  default:
    printf("\n Error... Not Available Method \n");
    break;
  }

}
