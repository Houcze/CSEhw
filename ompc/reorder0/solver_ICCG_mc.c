/*
 * solver_ICCG
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include "solver_ICCG_mc.h"
#include "allocate.h"

extern int
solve_ICCG_mc(int N, int NL, int NU, int *indexL, int *itemL, int *indexU, int *itemU,
		double *D, double *B, double *X, double *AL, double *AU,
		int NCOLORtot, int PEsmpTOT, int *SMPindex, 
		double EPS, int *ITR, int *IER, int *itemLU, double *ALU, double *XX)
{
	double **W;
	double VAL, BNRM2, WVAL, SW, RHO, BETA, RHO1, C1, DNRM2, ALPHA, ERR;
	int i, j, ic, ip, L, ip1;
	int R = 0;
	int Z = 1;
	int Q = 1;
	int P = 2;
	int DD = 3;

/*********
 * INIT. *
 *********/
        W = (double **)allocate_matrix(sizeof(double *), 4, N+128);

#pragma omp parallel for private (i)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
			X[i] = 0.0;
			W[1][i] = 0.0;
			W[2][i] = 0.0;
			W[3][i] = 0.0;
		}
	}


	for(ic=0; ic<NCOLORtot; ic++) {
#pragma omp parallel for private (ip1, i, VAL, j)
		for(ip=0; ip<PEsmpTOT; ip++) {
			ip1 = ip * NCOLORtot + ic;
			for(i=SMPindex[ip1]; i<SMPindex[ip1+1]; i++) {
				VAL = D[i];
				for(j=indexL[i]; j<indexL[i+1]; j++) {
					VAL = VAL - AL[j]*AL[j] * W[DD][itemL[j] - 1];
				}

				W[DD][i] = 1.0 / VAL;
			}
		}
	}

/**************************
 * {r0} = {b} - {A}{xini} *
 **************************/

//double *ALU;
//int *itemLU;
//double *XX;
//ALU = (double *)allocate_vector(sizeof(double),N*6);
//itemLU = (int *)allocate_vector(sizeof(int),N*6);
//XX = (double *)allocate_vector(sizeof(double),N*6);
/*
for(int i=0; i<N; i++)
{
	for(int j=0; j<indexL[i+1]-indexL[i]; j++)
	{
		ALU[6 * i + j] = AL[j];
	}
	for(int j=indexL[i+1]-indexL[i]; j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j++)
	{
		ALU[6 * i + j] = AU[j];
	}
	for (int j=indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j<6; j++)
	{
		ALU[6 * i + j] = 0.;
	}
}

for(int i=0; i<N; i++)
{
	for(int j=0; j<indexL[i+1]-indexL[i]; j++)
	{
		itemLU[6 * i + j] = itemL[j] + 1;
	}
	for(int j=indexL[i+1]-indexL[i]; j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j++)
	{
		itemLU[6 * i + j] = itemU[j] + 1;
	}
	for (int j=indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j<6; j++)
	{
		itemLU[6 * i + j] = 0.;
	}
}

for(int i=0; i<N; i++)
{
	for(int j=0; j<indexL[i+1]-indexL[i]; j++)
	{
		XX[itemLU[6 * i + j]] = X[itemL[j]-1];
	}
	for(int j=indexL[i+1]-indexL[i]; j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j++)
	{
		XX[itemLU[6 * i + j]] = X[itemU[j]-1];
	}
	for (int j=indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j<6; j++)
	{
		XX[6 * i + j] = 0.;
	}
}
*/
/*
int jj=0;
for(int i=0; i<N; i++)
{
	for(int j=0; j<indexL[i+1]-indexL[i]; j++)
	{
		ALU[6 * i + j] = AL[jj];
		itemLU[6 * i + j] = itemL[jj] - 1;
		XX[itemLU[6 * i + j]] = X[itemL[jj]-1];
		jj=jj+1;
	}
	for(int j=indexL[i+1]-indexL[i]; j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j++)
	{
		ALU[6 * i + j] = AU[jj];
		itemLU[6 * i + j] = itemU[jj] - 1;
		XX[itemLU[6 * i + j]] = X[itemU[jj]-1];
		jj=jj+1;
	}
	for (int j=indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j<6; j++)
	{
		ALU[6 * i + j] = 0.;
		itemLU[6 * i + j] = 0.;
		XX[6 * i + j] = 0.;
	}
}
*/
/*
for(i=0; i<N; i++)
{
	for(j=0; j<6; j++)
	{
		if(j<indexL[i+1]-indexL[i])
		{
			ALU[6 * i + j] = AL[indexL[i]+j];
			itemLU[6 * i + j] = itemL[indexL[i]+j] - 1;
			XX[itemLU[6 * i + j]] = X[itemL[indexL[i]+j]-1];
			
		}else if (j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i])
		{
			ALU[6 * i + j] = AU[indexU[i]+j-indexL[i+1]+indexL[i]];
			itemLU[6 * i + j] = itemU[indexU[i]+j-indexL[i+1]+indexL[i]] - 1;
			XX[itemLU[6 * i + j]] = X[itemU[indexU[i]+j-indexL[i+1]+indexL[i]] - 1];
		}
		
	}
}
*/
/*
int *delta;
delta = (int *) allocate_vector(sizeof(int), 2*N);

for(i=0; i<N; i++)
{
	delta[2*i] = indexL[i+1]-indexL[i];
	delta[2*i+1] = indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i];
}
*/
/*
for(i=0; i<N; i++)
{
	for(j=0; j<indexL[i+1]-indexL[i]; j++)
	{
		ALU[6 * i + j] = AL[indexL[i]+j];
		itemLU[6 * i + j] = itemL[indexL[i]+j] - 1;
		XX[itemLU[6 * i + j]] = X[itemL[indexL[i]+j]-1];
		
	}
	for(j=indexL[i+1]-indexL[i];j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i];j++)
	{
		ALU[6 * i + j] = AU[indexU[i]+j-indexL[i+1]+indexL[i]];
		itemLU[6 * i + j] = itemU[indexU[i]+j-indexL[i+1]+indexL[i]] - 1;
		XX[itemLU[6 * i + j]] = X[itemU[indexU[i]+j-indexL[i+1]+indexL[i]] - 1];		
	}
}
*/

#pragma omp parallel for private (i, VAL, j)
/*
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
			VAL = D[i] * X[i];
			for(j=0; j<6; j++) {
				VAL += ALU[6 * i + j] * XX[itemLU[6 * i + j]];
				//printf("itemL[j] is %d\n", itemL[j]);
			}
			//printf("-----------------\n");
			//for(j=indexU[i]; j<indexU[i+1]; j++) {
			//	VAL += AU[j] * X[itemU[j]-1];
				//printf("itemU[j] is %d\n", itemU[j]);
			//}
			//printf("-----------------\n");
			W[R][i] = B[i] - VAL;
		}
	}
	*/
	for(i=0; i<N; i++) {
		VAL = D[i] * X[i];
		
		for(j=0; j<6; j++) {
			VAL += ALU[6 * i + j] * XX[itemLU[6 * i + j]];
			//printf("itemL[j] is %d\n", itemL[j]);
		}
		
		//printf("-----------------\n");
		//for(j=indexU[i]; j<indexU[i+1]; j++) {
		//	VAL += AU[j] * X[itemU[j]-1];
			//printf("itemU[j] is %d\n", itemU[j]);
		//}
		//printf("-----------------\n");
		W[R][i] = B[i] - VAL;
		
	}
	BNRM2 = 0.0;


#pragma omp parallel for private (i) reduction (+:BNRM2)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
		  BNRM2 += B[i]*B[i];
		}
	}


/************************************************************** ITERATION */
	*ITR = N;

	for(L=0; L<(*ITR); L++) {

/*******************
 * {z} = [Minv]{r} *
 *******************/

#pragma omp parallel for private(i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
				W[Z][i] = W[R][i];
			}
		}

#pragma omp parallel private (ic, ip, ip1, i, WVAL, j)
		for(ic=0; ic<NCOLORtot; ic++) {
#pragma omp for
			for(ip=0; ip<PEsmpTOT; ip++) {
				ip1 = ip * NCOLORtot + ic;
				for(i=SMPindex[ip1]; i<SMPindex[ip1+1]; i++) {
					WVAL = W[Z][i];
					for(j=indexL[i]; j<indexL[i+1]; j++) {
						WVAL -= AL[j] * W[Z][itemL[j]-1];
					}
					W[Z][i] = WVAL * W[DD][i];

				}				
			}
		}

#pragma omp parallel private (ic, ip, ip1, i, SW, j)
		for(ic=0; ic<NCOLORtot; ic++) {
#pragma omp for
			for(ip=0; ip<PEsmpTOT; ip++) {
				ip1 = ip * NCOLORtot + NCOLORtot - 1 - ic;
				for(i=SMPindex[ip1]; i<SMPindex[ip1+1]; i++) {
					SW = 0.0;
					for(j=indexU[i]; j<indexU[i+1]; j++) {
						SW += AU[j] * W[Z][itemU[j]-1];
					}
					W[Z][i] -= W[DD][i] * SW;
				}
			}
		}

/////////////////////////////////////////////////////////////////////////////


/****************
 * RHO = {r}{z} *
 ****************/
		RHO = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:RHO)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
				RHO += W[R][i] * W[Z][i];
			}
		}

               

/********************************
 * {p}  = {z} if      ITER=0    *
 * BETA = RHO / RHO1  otherwise *
 ********************************/
		if(L == 0) {
#pragma omp parallel for private (ip, i)
			for(ip=0; ip<PEsmpTOT; ip++) {
				for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
					W[P][i] = W[Z][i];
				}
			}
		} else {
			BETA = RHO / RHO1;
#pragma omp parallel for private (ip, i)
			for(ip=0; ip<PEsmpTOT; ip++) {
				for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
					W[P][i] = W[Z][i] + BETA * W[P][i];
				}
			}
		}

/****************
 * {q} = [A]{p} *
 ****************/
/*
double s1, s2;
s1 = omp_get_wtime();
#pragma omp parallel for private (i, VAL, j)
	for(i=0; i<N; i++)
	{
		VAL = D[i] * W[P][i];
		for(j=0; j<indexL[i+1]-indexL[i]+indexU[i+1]-indexU[i]; j++)
		{

			VAL += ALU[6*i + j] * W[P][itemLU[6*i+j]];
				//printf("i is %d, j is %d, itemLU is %d\n", i, j, itemLU[6*i+j]);
			
		}
		W[Q][i] = VAL;
	}
s2 = omp_get_wtime();
fprintf(stdout, "After: %16.6e sec. (solver)\n", s2 - s1);
*/
#pragma omp parallel for private (ip1, i, VAL, j)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
			VAL = D[i] * W[P][i];
			for(j=indexL[i]; j<indexL[i+1]; j++) {
				VAL += AL[j] * W[P][itemL[j]-1];
			}
			for(j=indexU[i]; j<indexU[i+1]; j++) {
				VAL += AU[j] * W[P][itemU[j]-1];
			}
			W[Q][i] = VAL;
		}
	}
/*
#pragma omp parallel for private (ip1, i, VAL, j)
	for(ip=0; ip<PEsmpTOT; ip++) {
		for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
			VAL = D[i] * W[P][i];
			for(j=indexL[i]; j<indexL[i+1]; j++) {
				//printf("j is %d, indexLU[6*i+j] is %d, itemL[i] is %d\n", j, itemLU[6*i+j], itemL[j]-1);
				VAL += AL[j] * W[P][itemL[j]-1];
			}
			for(j=indexU[i]; j<indexU[i+1]; j++) {
				VAL += AU[j] * W[P][itemU[j]-1];
			}
			W[Q][i] = VAL;
		}
	}
	
*/
/************************
 * ALPHA = RHO / {p}{q} *
 ************************/
		C1 = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:C1)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
				C1 += W[P][i] * W[Q][i];
			}
		}

		ALPHA = RHO / C1;

/***************************
 * {x} = {x} + ALPHA * {p} *
 * {r} = {r} - ALPHA * {q} *
 ***************************/
#pragma omp parallel for private (ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
				X[i]    += ALPHA * W[P][i];
			}
		}

#pragma omp parallel for private (ip, i)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
				W[R][i] -= ALPHA * W[Q][i];
			}
		}

		DNRM2 = 0.0;
#pragma omp parallel for private (ip, i) reduction (+:DNRM2)
		for(ip=0; ip<PEsmpTOT; ip++) {
			for(i=SMPindex[ip*NCOLORtot]; i<SMPindex[(ip+1)*NCOLORtot]; i++) {
			  DNRM2 += W[R][i]*W[R][i];
			}
		}

		ERR = sqrt(DNRM2/BNRM2);
		if( (L+1)%100 ==1) {
			fprintf(stdout, "%5d%16.6e\n", L+1, ERR);
		}



		if(ERR < EPS) {
			*IER = 0;
			goto N900;
		} else {
			RHO1 = RHO;
		}
	}
	*IER = 1;

N900:
	fprintf(stdout, "%5d%16.6e\n", L+1, ERR);
	*ITR = L;

	free(W);

	return 0;
}
