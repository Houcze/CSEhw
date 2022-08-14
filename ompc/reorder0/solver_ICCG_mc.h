#ifndef __H_SOLVER_ICCG_MC
#define __H_SOLVER_ICCG_MC

extern int
solve_ICCG_mc(int N, int NL, int NU, int *indexLnew, int *itemLnew, int *indexUnew, int *itemUnew,
		double *D, double *B, double *X, double *ALnew, double *AUnew,
		int NCOLOR, int PEsmpTOT, int *COLOR, 
		double EPS, int *ITR, int *IER);

#endif /* __H_SOLVER_ICCG_MC */
