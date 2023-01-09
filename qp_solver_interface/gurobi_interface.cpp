#include "gurobi_interface.hpp"
#include "step_param.h"

extern "C"
{
#include "acado_common.h"
#include <stdio.h>
#include <time.h>

#include <stdlib.h>
#include "gurobi_c.h"

}

#define atoa(x) #x

int solver_gurobi() {

  double *H = acadoWorkspace.H;
  double *g = acadoWorkspace.g;
  double *A = acadoWorkspace.A;
  double *lb = acadoWorkspace.lb;
  double *ub = acadoWorkspace.ub;
  double *lbA = acadoWorkspace.lbA;
  double *ubA = acadoWorkspace.ubA;
  double *evH = acadoWorkspace.evH;


  GRBenv   *env   = NULL;
  GRBenv   *env2   = NULL;
  GRBmodel *model = NULL;
  GRBmodel *model2 = NULL;
  GRBmodel *model3 = NULL;

	int       error = 0;
  double    sol[Nstep*2 + Nobs*2];
  double    sol2[Nstep*2];
  double    sol3[Nstep*2];
	double    solpi[(Ntot+Nstep*2)*2];

	int       indA[Nstep*2];
	double    valA[Nstep*2];
	int       ind[1];
	double    val[1];

	int       qrow[Nstep*2*Nstep*2];
	int       qcol[Nstep*2*Nstep*2];
	double    qval[Nstep*2*Nstep*2];

	int       optimstatus;
	double    objval;

	/* Create environment */
  error = GRBloadenv(&env, NULL); //"qp.log");
	if (error) goto QUIT;

	/* Create an empty model */
	error = GRBnewmodel(env, &model, "qp", 0, NULL, NULL, NULL, NULL, NULL);
	if (error) goto QUIT;

	/* Add variables */
	error = GRBaddvars(model, Nstep*2, 0, NULL, NULL, NULL, NULL, lb, ub, NULL,
	                  NULL);
	if (error) goto QUIT;

	char vtype[Nobs*2];
	for (int i=0; i<Nobs*2; i++)	vtype[i] = GRB_BINARY;
	error = GRBaddvars(model, Nobs*2, 0, NULL, NULL, NULL, NULL, NULL, NULL, vtype,
                     NULL);
	if (error) goto QUIT;

	/* silence output */
	error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_OUTPUTFLAG, 0); //1);
	if (error) goto QUIT;

	/* Quadratic objective terms */
	for (int i = 0; i < Nstep*2; i++) {
		for (int j = 0; j < Nstep*2; j++) {
			int nb = i*Nstep*2+j;
			qrow[nb] = i;
			qcol[nb] = j;
			qval[nb] = H[nb]*0.5;
		}
	}
	error = GRBaddqpterms(model, Nstep*2*Nstep*2, qrow, qcol, qval);
	if (error) goto QUIT;

	/* Linear objective term */
	for (int i = 0; i < Nstep*2; i++) {
		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, i, g[i]);
		if (error) goto QUIT;
	}

	/* constraints */
	for (int i = 0; i < Ntot; i++) {
		for (int j = 0; j < Nstep*2; j++) {
			indA[j] = j;
			valA[j] = A[i*Nstep*2+j];
		}
			error = GRBaddconstr(model, Nstep*2, indA, valA, GRB_LESS_EQUAL, ubA[i], NULL);
			if (error) goto QUIT;
			error = GRBaddconstr(model, Nstep*2, indA, valA, GRB_GREATER_EQUAL, lbA[i], NULL);
			if (error) goto QUIT;
	}
	for (int i = 0; i < Nstep*2; i++) {
		ind[0] = i;
		val[0] = 1;
		error = GRBaddconstr(model, 1, ind, val, GRB_LESS_EQUAL, ub[i], NULL);
		if (error) goto QUIT;
		error = GRBaddconstr(model, 1, ind, val, GRB_GREATER_EQUAL, lb[i], NULL);
		if (error) goto QUIT;
	}

  int       indC[Nstep*2+2];
	double    valC[Nstep*2+2];
	int       indT[2];
	double    valT[2];

  for (int i = 0; i < Nobs; i++) {
		int idx;

		for(int s=5; s<6; s++) {
			idx = (Nstep*8)+i*15+2+s;
			for (int j = 0; j < Nstep*2; j++) {
				indC[j] = j;
				valC[j] = A[idx*Nstep*2+j];
			}
			indC[Nstep*2] = i*2 +0 +Nstep*2;
			valC[Nstep*2] = -1e+6 * 1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
			valC[Nstep*2+1] = -1e+6 * 1;
			error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_LESS_EQUAL, -evH[idx-(Nstep*8)], NULL);
			if (error) goto QUIT;
		}

		for(int s=6; s<7; s++) {
			idx = (Nstep*8)+i*15+2+s;
			for (int j = 0; j < Nstep*2; j++) {
				indC[j] = j;
				valC[j] = A[idx*Nstep*2+j];
			}
      indC[Nstep*2] = i*2 +0 +Nstep*2;
      valC[Nstep*2] = +1e+6 * -1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
      valC[Nstep*2+1] = +1e+6 * -1;
			error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_GREATER_EQUAL, -evH[idx-(Nstep*8)] - +1e+6*2, NULL);
			if (error) goto QUIT;
		}

		for(int s=7; s<9; s++) {
			idx = (Nstep*8)+i*15+2+s;
			for (int j = 0; j < Nstep*2; j++) {
				indC[j] = j;
				valC[j] = A[idx*Nstep*2+j];
			}
			indC[Nstep*2] = i*2 +0 +Nstep*2;
			valC[Nstep*2] = -1e+6 * -1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
			valC[Nstep*2+1] = -1e+6 * +1;
			error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_LESS_EQUAL, -evH[idx-(Nstep*8)] - -1e+6, NULL);
			if (error) goto QUIT;
		}
    for(int s=9; s<10; s++) {
      idx = (Nstep*8)+i*15+2+s;
      for (int j = 0; j < Nstep*2; j++) {
        indC[j] = j;
        valC[j] = A[idx*Nstep*2+j];
      }
      indC[Nstep*2] = i*2 +0 +Nstep*2;
      valC[Nstep*2] = +1e+6 * -1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
      valC[Nstep*2+1] = +1e+6 * +1;
      error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_GREATER_EQUAL, -evH[idx-(Nstep*8)] - +1e+6, NULL);
      if (error) goto QUIT;
    }

		for(int s=10; s<11; s++) {
			idx = (Nstep*8)+i*15+2+s;
			for (int j = 0; j < Nstep*2; j++) {
				indC[j] = j;
				valC[j] = A[idx*Nstep*2+j];
			}
      indC[Nstep*2] = i*2 +0 +Nstep*2;
      valC[Nstep*2] = +1e+6 * +1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
      valC[Nstep*2+1] = +1e+6 * -1;
			error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_GREATER_EQUAL, -evH[idx-(Nstep*8)] - +1e+6, NULL);
			if (error) goto QUIT;
		}
    for(int s=11; s<12; s++) {
      idx = (Nstep*8)+i*15+2+s;
      for (int j = 0; j < Nstep*2; j++) {
        indC[j] = j;
        valC[j] = A[idx*Nstep*2+j];
      }
      indC[Nstep*2] = i*2 +0 +Nstep*2;
      valC[Nstep*2] = -1e+6 * +1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
      valC[Nstep*2+1] = -1e+6 * -1;
      error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_LESS_EQUAL, -evH[idx-(Nstep*8)] - -1e+6, NULL);
      if (error) goto QUIT;
    }
    for(int s=12; s<13; s++) {
      idx = (Nstep*8)+i*15+2+s;
      for (int j = 0; j < Nstep*2; j++) {
        indC[j] = j;
        valC[j] = A[idx*Nstep*2+j];
      }
      indC[Nstep*2] = i*2 +0 +Nstep*2;
      valC[Nstep*2] = +1e+6 * +1;
      indC[Nstep*2+1] = i*2 +1 +Nstep*2;
      valC[Nstep*2+1] = +1e+6 * -1;
      error = GRBaddconstr(model, Nstep*2+2, indC, valC, GRB_GREATER_EQUAL, -evH[idx-(Nstep*8)] - +1e+6, NULL);
      if (error) goto QUIT;
    }
	}

	/* Optimize model */
	error = GRBoptimize(model);
	if (error) goto QUIT;

	/* Capture solution information */
	error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
	if (error) goto QUIT;

	if (optimstatus != GRB_OPTIMAL) {
		printf("%s\n", "INFEASIBLE!");
		printf("%d\n", optimstatus);
	}

  //int sC;
  //error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &sC);
	//if (error) goto QUIT;
  //error = GRBsetintattr(model, GRB_INT_PAR_SOLUTIONNUMBER, sC-1);
  //if (error) goto QUIT;

	error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, Nstep*2 + Nobs*2, sol);
	if (error) goto QUIT;

  for (int i = 0; i < Nobs; i++) {
		printf("STEP %2d result - ", i);
		for (int j=0; j<2; j++) printf("%d ", (int)sol[i*2+j+Nstep*2]);
		printf(" | ");

    int c1 = sol[i*2+0+Nstep*2];
    int c2 = sol[i*2+1+Nstep*2];

		if ( (c1==0) && (c2==0) )
			printf("%s ", "S-");
		if ( (c1==1) && (c2==1) )
			printf("%s ", "S+");
		if ( (c1==1) && (c2==0) )
			printf("%s ", "Ey-");
		if ( (c1==0) && (c2==1) )
			printf("%s ", "Ey+");
		printf("\n");
	}

  model2 = GRBfixedmodel(model);
  error = GRBoptimize(model2);
  if (error) goto QUIT;
  /* Capture solution information */
	error = GRBgetintattr(model2, GRB_INT_ATTR_STATUS, &optimstatus);
	if (error) goto QUIT;

	if (optimstatus != GRB_OPTIMAL) {
		printf("%s\n", "INFEASIBLE!");
		printf("%d\n", optimstatus);
    goto QUIT;
	}

  /*
  int cNum;
  error = GRBgetintattr(model2, GRB_INT_ATTR_NUMCONSTRS, &cNum);
  if (error) goto QUIT;
  printf("cNum: %d\n", cNum);
  double slack[2000];
  error = GRBgetdblattrarray(model2, GRB_DBL_ATTR_SLACK, 0, cNum, slack);
  if (error) goto QUIT;
  for (int i=0; i<cNum; i++)
    printf("slack: %f\n", slack[i]);
  */

  error = GRBgetdblattrarray(model2, GRB_DBL_ATTR_PI, 0, (Ntot+Nstep*2)*2, solpi);
  if (error) goto QUIT;
  error = GRBgetdblattrarray(model2, GRB_DBL_ATTR_X, 0, Nstep*2, sol2);
	if (error) goto QUIT;

	if (optimstatus == GRB_OPTIMAL) {
    /* Write out optimization result */
		for (int i = 0; i < Nstep*2; i++) {
      acadoWorkspace.x[i] = sol2[i];
      // printf("acadoWorkspace.x[%d]: %f\n", i, acadoWorkspace.x[i]);
    }
		for (int i = 0; i < Ntot; i++) {
			if ((solpi[i*2]>1e-12) || (solpi[i*2]<-1e-12)) {
				acadoWorkspace.y[i+Nstep*2] = solpi[i*2];
				continue;
			}
			if ((solpi[i*2+1]>1e-12) || (solpi[i*2+1]<-1e-12)) {
				acadoWorkspace.y[i+Nstep*2] = solpi[i*2+1];
				continue;
			}
			acadoWorkspace.y[i+Nstep*2] = 0.0; //solpi[i*2];
		}
		for (int i = Ntot; i < Ntot+Nstep*2; i++) {
			if ((solpi[i*2]>1e-12) || (solpi[i*2]<-1e-12)) {
				acadoWorkspace.y[i-Ntot] = solpi[i*2];
				continue;
			}
			if ((solpi[i*2+1]>1e-12) || (solpi[i*2+1]<-1e-12)) {
				acadoWorkspace.y[i-Ntot] = solpi[i*2+1];
				continue;
			}
			acadoWorkspace.y[i-Ntot] = 0.0;
		}

	} else if (optimstatus == GRB_INF_OR_UNBD) {
		printf("Model is infeasible or unbounded\n");
	} else {
		printf("Optimization was stopped early\n");
	}

	QUIT:

	/* Error reporting */
	if (error) {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		//exit(1);
	}

	/* Free model */
  GRBfreemodel(model);
  GRBfreemodel(model2);
	/* Free environment */
	GRBfreeenv(env);
};
