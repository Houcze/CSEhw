/*
 * POI_GEN
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "struct_ext.h"
#include "pcg_ext.h"
#include "mc.h"
#include "cm.h"
#include "rcm.h"
#include "cmrcm.h"
#include "poi_gen.h"
#include "allocate.h"

extern int
POI_GEN(void)
{
	int nn, nnp;
	int ic0, ic01, ic02, icN1, icN2, icN3, icN4, icN5, icN6, icNS, icNE, ik0;
	int icN10, icN20, icN30, icN40, icN50, icN60;
	int i, j, k, ib, ic, ip, icel, icelN, icou, icol, icouG;
	int ii, jj, kk, nn1, num, nr, j0, j1;
	double coef, VOL0, S1t, E1t;
	int isL, ieL, isU, ieU;
	int id1, id2;
	int j_new;


	NL = 6;
	NU = 6;

        BFORCE = (double *)allocate_vector(sizeof(double),ICELTOT);
        D      = (double *)allocate_vector(sizeof(double),ICELTOT);
        PHI    = (double *)allocate_vector(sizeof(double),ICELTOT);
        indexLnew = (int *)allocate_vector(sizeof(int),ICELTOT+1);
        indexUnew = (int *)allocate_vector(sizeof(int),ICELTOT+1);

        INL    = (int *)allocate_vector( sizeof(int),ICELTOT);
        INU    = (int *)allocate_vector( sizeof(int),ICELTOT);
        INLnew = (int *)allocate_vector( sizeof(int),ICELTOT);
        INUnew = (int *)allocate_vector( sizeof(int),ICELTOT);
        IAL    = (int **)allocate_matrix(sizeof(int),ICELTOT,NL);
        IAU    = (int **)allocate_matrix(sizeof(int),ICELTOT,NU);

        indexL = (int *)allocate_vector(sizeof(int),ICELTOT+1);
        indexU = (int *)allocate_vector(sizeof(int),ICELTOT+1);
        indexLnew_org = (int *)allocate_vector(sizeof(int),ICELTOT+1);
        indexUnew_org = (int *)allocate_vector(sizeof(int),ICELTOT+1);

		itemLU = (int *)allocate_vector(sizeof(int),ICELTOT*6);
		XX = (double *)allocate_vector(sizeof(double),ICELTOT*6);
		ALU = (double *)allocate_vector(sizeof(double),ICELTOT*6);
		memset(itemLU, 0, sizeof(int) * ICELTOT*6);
		memset(XX, 0, sizeof(double) * ICELTOT*6);
		memset(ALU, 0, sizeof(double) * ICELTOT*6);

	for(j=0; j<ICELTOT; j++) {
		INL[j] = 0;
		INU[j] = 0;
		INLnew[j] = 0;
		INUnew[j] = 0;
		for(i=0; i<6; i++) {
			IAL[j][i] = 0;
			IAU[j][i] = 0;
		}
	}

	for(i=0;i<=ICELTOT;i++){
		indexL[i]=0;
		indexU[i]=0;
		indexLnew_org[i]=0;
		indexUnew_org[i]=0;
	}


/*********************************
 * INTERIOR & NEUMANN boundary's *
 *********************************/

	for(icel=0; icel<ICELTOT; icel++) {
		icN1 = NEIBcell[icel][0];
		icN2 = NEIBcell[icel][1];
		icN3 = NEIBcell[icel][2];
		icN4 = NEIBcell[icel][3];
		icN5 = NEIBcell[icel][4];
		icN6 = NEIBcell[icel][5];

		if(icN5 != 0) {
			icou = INL[icel] + 1;
			IAL[icel][icou-1] = icN5;
			INL[icel]         = icou;
		}

		if(icN3 != 0) {
			icou = INL[icel] + 1;
			IAL[icel][icou-1] = icN3;
			INL[icel]         = icou;
		}

		if(icN1 != 0) {
			icou = INL[icel] + 1;
			IAL[icel][icou-1] = icN1;
			INL[icel]         = icou;
		}

		if(icN2 != 0) {
			icou = INU[icel] + 1;
			IAU[icel][icou-1] = icN2;
			INU[icel]         = icou;
		}

		if(icN4 != 0) {
			icou = INU[icel] + 1;
			IAU[icel][icou-1] = icN4;
			INU[icel]         = icou;
		}

		if(icN6 != 0) {
			icou = INU[icel] + 1;
			IAU[icel][icou-1] = icN6;
			INU[icel]         = icou;
		}
	}

/*****************
 * MULTICOLORING *
 *****************/

	OLDtoNEW = (int *) allocate_vector(sizeof(int),ICELTOT);
	NEWtoOLD = (int *) allocate_vector(sizeof(int),ICELTOT);
	OLDtoNEWnew = (int *) allocate_vector(sizeof(int),ICELTOT);
	NEWtoOLDnew = (int *) allocate_vector(sizeof(int),ICELTOT);
	COLORindex = (int *) allocate_vector(sizeof(int),ICELTOT+1);

	memset(OLDtoNEW,0,sizeof(int)*ICELTOT);
	memset(NEWtoOLD,0,sizeof(int)*ICELTOT);
	memset(OLDtoNEWnew,0,sizeof(int)*ICELTOT);
	memset(NEWtoOLDnew,0,sizeof(int)*ICELTOT);
	memset(COLORindex,0,sizeof(int)*(ICELTOT+1));

N111:
	fprintf(stdout, "\n\nYou have%8d elements\n", ICELTOT);
	fprintf(stdout, "How many colors do you need ?\n");
	fprintf(stdout, "  #COLOR must be more than 2 and\n");
	fprintf(stdout, "  #COLOR must not be more than%8d\n", ICELTOT);
	fprintf(stdout, "   CM if #COLOR .eq. 0\n");
	fprintf(stdout, "  RCM if #COLOR .eq.-1\n");
	fprintf(stdout, "CMRCM if #COLOR .lt.-1\n");
	fprintf(stdout, "=>\n");

	if(NCOLORtot == 1 || NCOLORtot > ICELTOT) goto N111;

	if(NCOLORtot > 0) {
		MC(ICELTOT, NL, NU, INL, IAL, INU, IAU,
				&NCOLORtot, COLORindex, NEWtoOLD, OLDtoNEW);
		fprintf(stdout, "\n###  MultiColoring\n");
	} else if(NCOLORtot == 0) {
		CM(ICELTOT, NL, NU, INL, IAL, INU, IAU,
				&NCOLORtot, COLORindex, NEWtoOLD, OLDtoNEW);
		fprintf(stdout, "\n###  CM\n");
	} else if(NCOLORtot ==-1) {
		RCM(ICELTOT, NL, NU, INL, IAL, INU, IAU,
				&NCOLORtot, COLORindex, NEWtoOLD, OLDtoNEW);
		fprintf(stdout, "\n###  RCM\n");
	} else if(NCOLORtot < -1) {
		CMRCM(ICELTOT, NL, NU, INL, IAL, INU, IAU,
				&NCOLORtot, COLORindex, NEWtoOLD, OLDtoNEW);
		fprintf(stdout, "\n###  CMRCM\n");
	} 

	fprintf(stdout, "\n### FINAL COLOR NUMBER%8d\n\n", NCOLORtot);


/************************************
* SMPindex, SMPindex_new, SMPindexG *
*************************************/

	SMPindex = (int *) allocate_vector(sizeof(int), NCOLORtot * PEsmpTOT + 1);
        memset(SMPindex, 0, sizeof(int)*(NCOLORtot*PEsmpTOT+1));

	for(ic=1; ic<=NCOLORtot; ic++) {
		nn1 = COLORindex[ic] - COLORindex[ic-1];
		num = nn1 / PEsmpTOT;
		nr  = nn1 - PEsmpTOT * num;
		for(ip=1; ip<=PEsmpTOT; ip++) {
			if(ip <= nr) {
				SMPindex[(ic-1)*PEsmpTOT+ip] = num + 1;
			} else {
				SMPindex[(ic-1)*PEsmpTOT+ip] = num;
			}
		}
	}

	SMPindex_new = (int *) allocate_vector(sizeof(int), NCOLORtot * PEsmpTOT + 1);
        memset(SMPindex_new, 0, sizeof(int)*(NCOLORtot*PEsmpTOT+1));

	for(ic=1; ic<=NCOLORtot; ic++) {
		for(ip=1; ip<=PEsmpTOT; ip++) {
			j1 =(ic-1)*PEsmpTOT + ip;
			j0 =j1-1;
			SMPindex_new[(ip-1)*NCOLORtot+ic] = SMPindex[j1];
			SMPindex[j1] = SMPindex[j0] + SMPindex[j1];
		}
	}

	for(ip=1; ip<=PEsmpTOT; ip++) {
		for(ic=1; ic<=NCOLORtot; ic++) {
			j1 = (ip-1) * NCOLORtot + ic;
			j0 = j1 - 1;
			SMPindex_new[j1] += SMPindex_new[j0];
		}
	}
	for(ip=0; ip<=PEsmpTOT; ip++) {
		for(ic=1; ic<=NCOLORtot; ic++) {
			j1 = (ip-1) * NCOLORtot + ic;
		}
	}



	SMPindexG = (int *) allocate_vector(sizeof(int), PEsmpTOT + 1);
        memset(SMPindexG, 0, sizeof(int)*(PEsmpTOT+1));

	num = ICELTOT / PEsmpTOT;
	nr = ICELTOT - num * PEsmpTOT;
	for(ip=1; ip<=PEsmpTOT; ip++) {
		SMPindexG[ip] = num;
		if(ip <= nr) {
			SMPindexG[ip] += 1;
		}
	}

	for(ip=1; ip<=PEsmpTOT; ip++){
		SMPindexG[ip] = SMPindexG[ip-1] +SMPindexG[ip];
	}

/***************************
* OLDtoNEWnew, NEWtoOLDnew *
****************************/
	if(NFLAG == 0){
		for(i=0; i<ICELTOT; i++) {
			OLDtoNEWnew[i] = 0;
			NEWtoOLDnew[i] = 0;
		}
	}else {
#pragma omp parallel for private (icel, j)
		for(ip=1; ip<=PEsmpTOT; ip++){
			for(icel = SMPindex_new[(ip-1)*NCOLORtot]+1; icel<=SMPindex_new[ip*NCOLORtot]; icel++){
				OLDtoNEWnew[icel-1] = 0;
				NEWtoOLDnew[icel-1] = 0;
			}
		}
	}

	for(ip=0; ip<PEsmpTOT; ip++){
		for(ic=0; ic<NCOLORtot; ic++){
			icNS = SMPindex_new[ip*NCOLORtot + ic];
			ic01 = SMPindex[ic*PEsmpTOT + ip];
			ic02 = SMPindex[ic*PEsmpTOT + ip+1];
			icou = 0;
			for(k=ic01; k<ic02; k++){
				icel=NEWtoOLD[k];
				icou = icou +1;
				icelN=icNS+icou;
				OLDtoNEWnew[icel-1] = icelN;
				NEWtoOLDnew[icelN-1]= icel;


			}
		}
	}


/********************************************
* 1D ordering: indexL, indexU, itemL, itemU *
*********************************************/


        for(i=1; i<=ICELTOT; i++){
                indexL[i]=indexL[i-1]+INL[i-1];
                indexU[i]=indexU[i-1]+INU[i-1];
        }

        NPL  = indexL[ICELTOT];
        NPU  = indexU[ICELTOT];

        itemL = (int *)allocate_vector(sizeof(int),NPL);
        itemU = (int *)allocate_vector(sizeof(int),NPU);
        itemLnew    = (int *)allocate_vector(sizeof(int),NPL);
        itemUnew    = (int *)allocate_vector(sizeof(int),NPU);
        ALnew       = (double *)allocate_vector(sizeof(double),NPL);
        AUnew       = (double *)allocate_vector(sizeof(double),NPU);

	memset(itemL, 0, sizeof(int)*NPL);
	memset(itemU, 0, sizeof(int)*NPU);

        for(i=1; i<=ICELTOT; i++){
                for(k=1;k<=INL[i-1];k++){
                        kk=k+indexL[i-1];
                        itemL[kk-1]=IAL[i-1][k-1];
                }
                for(k=1;k<=INU[i-1];k++){
                        kk=k+indexU[i-1];
                        itemU[kk-1]=IAU[i-1][k-1];
                }
        }


	for(ip=1; ip<=PEsmpTOT; ip++) {
        id1 = ip    *NCOLORtot;
        id2 = (ip-1)*NCOLORtot;

		for(icel=SMPindex_new[id2]+1; icel<=SMPindex_new[id1]; icel++) {
			ic0  = NEWtoOLDnew[icel-1];
			ik0  = OLDtoNEW[ic0-1];

			INLnew[icel-1] =INL[ik0-1];
			INUnew[icel-1] =INU[ik0-1];

		}
	}

        for(i=1; i<=ICELTOT; i++){
                indexLnew_org[i]=indexLnew_org[i-1]+INLnew[i-1];
                indexUnew_org[i]=indexUnew_org[i-1]+INUnew[i-1];

        }

	icou=0;
	for(icel=0; icel<ICELTOT; icel++){
		if(INL[icel] == 0){
			icou=icou+1;
		}
	}
	fprintf(stdout, "incompatible nodes:%d\n", icou);


        free(INL);
        free(INU);
       // free(INLnew);
       // free(INUnew);
        free(IAL);
        free(IAU);



/************
* ARRAY init.
*************/
	if(NFLAG == 0){
		for(i=0; i<ICELTOT; i++) {
			BFORCE[i] = 0.0;
			D[i]      = 0.0;
			PHI[i]    = 0.0;
		}
		for(i=0; i<=ICELTOT; i++) {
			indexLnew[i] = indexLnew_org[i];
			indexUnew[i] = indexUnew_org[i];

		}
		for(i=0; i<NPL; i++) {
			itemLnew[i] = 0;
			ALnew[i] = 0.0;
		} 
 		for(i=0; i<NPU; i++) {
			itemUnew[i] = 0;
			AUnew[i] = 0.0;
		}

	}else {

		indexLnew[0]=0;
		indexUnew[0]=0;
#pragma omp parallel for private (icel, j)
		for(ip=1; ip<=PEsmpTOT; ip++){
			for(icel = SMPindex_new[(ip-1)*NCOLORtot]+1; icel<=SMPindex_new[ip*NCOLORtot]; icel++){
				BFORCE[icel-1] = 0.0;
				PHI[icel-1] = 0.0;
				D[icel-1] = 0.0;
				indexLnew[icel]=indexLnew_org[icel];
				indexUnew[icel]=indexUnew_org[icel];

				for(j=indexLnew_org[icel-1];j<indexLnew_org[icel];j++){
					itemLnew[j]=0;
					ALnew[j] = 0.0;
				}
				for(j=indexUnew_org[icel-1];j<indexUnew_org[icel];j++){
					itemUnew[j]=0;
					AUnew[j] = 0.0;
				}
			}
		}
	}




/*************************************
* INTERIOR & NEUMANN BOUNDARY CELLs *
**************************************/

	S1t = omp_get_wtime();
#pragma omp parallel for private (icel,id1, id2, ic0,icN1, icN2, icN3, icN4, icN5, icN6, coef, j, ii, jj, j_new, kk, isL, ieL, isU, ieU, ik0, icN10, icN20, icN30, icN40, icN50, icN60)
	for(ip=1; ip<=PEsmpTOT; ip++) {
		id1= ip    *NCOLORtot;
		id2= (ip-1)*NCOLORtot;
		for(icel=SMPindex_new[id2]+1; icel<=SMPindex_new[id1]; icel++) {

			ic0  = NEWtoOLDnew[icel-1];
			ik0  = OLDtoNEW[ic0-1];


			icN10 = NEIBcell[ic0-1][0];
			icN20 = NEIBcell[ic0-1][1];
			icN30 = NEIBcell[ic0-1][2];
			icN40 = NEIBcell[ic0-1][3];
			icN50 = NEIBcell[ic0-1][4];
			icN60 = NEIBcell[ic0-1][5];

			isL=indexL[ik0-1];
			ieL=indexL[ik0  ];
			isU=indexU[ik0-1];
			ieU=indexU[ik0  ];


			if(icN50 != 0) {
				icN5 = OLDtoNEW[icN50-1];
				coef = RDZ * ZAREA;
				D[icel-1] -= coef;

				if(icN5 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN5) {
							j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new]    = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN50-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN5) {
							j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN50-1];
							break;
						}
					}
				}
			}

			if(icN30 != 0) {
				icN3 = OLDtoNEW[icN30-1];
				coef = RDY * YAREA;
				D[icel-1] -= coef;

				if(icN3 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN3) {
                                                        j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new] = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN30-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN3) {
                                                        j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN30-1];
							break;
						}
					}
				}
			}

			if(icN10 != 0) {
				icN1 = OLDtoNEW[icN10-1];
				coef = RDX * XAREA;
				D[icel-1] -= coef;

				if(icN1 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN1) {
                                                        j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new] = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN10-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN1) {
                                                        j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN10-1];
							break;
						}
					}
				}
			}

			if(icN20 != 0) {
				icN2 = OLDtoNEW[icN20-1];
				coef = RDX * XAREA;
				D[icel-1] -= coef;

				if(icN2 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN2) {
                                                        j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new] = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN20-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN2) {
                                                        j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN20-1];
							break;
						}
					}
				}
			}

			if(icN40 != 0) {
				icN4 = OLDtoNEW[icN40-1];
				coef = RDY * YAREA;
				D[icel-1] -= coef;

				if(icN4 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN4) {
                                                        j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new] = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN40-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN4) {
                                                        j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN40-1];
							break;
						}
					}
				}
			}

			if(icN60 != 0) {
				icN6 = OLDtoNEW[icN60-1];
				coef = RDZ * ZAREA;
				D[icel-1] -= coef;

				if(icN6 < ik0) {
					for(j=isL; j<ieL; j++) {
						if(itemL[j] == icN6) {
                                                        j_new=indexLnew[icel-1]+j-isL;
							ALnew[j_new] = coef;
							itemLnew[j_new] = OLDtoNEWnew[icN60-1];
							break;
						}
					}
				} else {
					for(j=isU; j<ieU; j++) {
						if(itemU[j] == icN6) {
                                                        j_new=indexUnew[icel-1]+j-isU;
							AUnew[j_new] = coef;
							itemUnew[j_new] = OLDtoNEWnew[icN60-1];
							break;
						}
					}
				}
			}


			ii = XYZ[ic0-1][0];
			jj = XYZ[ic0-1][1];
			kk = XYZ[ic0-1][2];

			BFORCE[icel-1] = - (double)(ii + jj + kk) * VOLCEL[ic0-1];
		}
	}



	E1t = omp_get_wtime();
	fprintf(stdout, "%16.6e sec. (assemble)\n", E1t - S1t);

/****************************
 * DIRICHLET BOUNDARY CELLs *
 ****************************/
/* TOP SURFACE */
	for(ib=0; ib<ZmaxCELtot; ib++) {
		ic0  = ZmaxCEL[ib];
		coef = 2.0 * RDZ * ZAREA;
		icel = OLDtoNEWnew[ic0-1];
		D[icel-1] -= coef;
	}



	return 0;
}
