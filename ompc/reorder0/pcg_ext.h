#ifndef __H_PCG
#define __H_PCG

	extern int N2;
	extern int NUmax, NLmax, NCOLORtot, NCOLORk, NU, NL;
	extern int METHOD, ORDER_METHOD, NFLAG;
	extern int NPL, NPU, NPL2, NPU2;

	extern double EPSICCG;

	extern double *D, *PHI, *BFORCE;
	//extern double **AL, **AU;
        double *ALnew, *AUnew;

        extern int *INL, *INU;
        extern int *SMPindex, *SMPindexG, *COLORindex;
        extern int *OLDtoNEW, *NEWtoOLD;

	extern int **IAL, **IAU;
	extern int *SMPindex_new;
	extern int *INLnew, *INUnew;

        extern int *indexL, *indexU;
        extern int *indexLnew, *indexUnew;
        extern int *indexLnew_org, *indexUnew_org;

        extern int *itemL, *itemU;
        extern int *itemLnew, *itemUnew;
        extern int *OLDtoNEWnew, *NEWtoOLDnew;
        // New modified
  

        extern int *INLU, *indexLU;
        extern int NPLU;
        extern double *XX;
        extern double *ALU;
        extern int *itemLU;

#endif /* __H_PCG */
