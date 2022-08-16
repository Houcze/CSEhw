#ifndef __H_PCG
#define __H_PCG

	static int N2 = 256;
	int NUmax, NLmax, NCOLORtot, NCOLORk, NU, NL, NLU;
	int METHOD, ORDER_METHOD, NFLAG;
	int NPL, NPU, NPL2, NPU2;

	double EPSICCG;

	double *D, *PHI, *BFORCE;
	//double *AL, *AU;
	double *ALnew, *AUnew;

	int *INL, *INU;
	int *INLnew, *INUnew;
	int *indexL, *indexU;
	int *indexLnew_org, *indexUnew_org;
	int *indexLnew, *indexUnew;
	int *SMPindex, *SMPindexG, *COLORindex;
	int *OLDtoNEW, *NEWtoOLD;

	int **IAL, **IAU;
	int **IALnew, **IAUnew;

	int *SMPindex_new;
        int *itemL, *itemU;
        int *itemLnew, *itemUnew;
	int *OLDtoNEWnew, *NEWtoOLDnew;
	// New modified

	double *AMAT;
	int *INLU, *indexLU, *itemLU;
	int NPLU;
	int *indexLU;

#endif /* __H_PCG */
