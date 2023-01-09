#ifndef SOLVER_SELECTOR
#define SOLVER_SELECTOR

#include "acado_common.h"
#include "acado_qpoases_interface.hpp"
#include "gurobi_interface.hpp"

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#define OPT_GUROBI_QP_SOLVER 0
#define OPT_QPOASES_QP_SOLVER 1

EXTERNC int acado_solve( void );

extern int solver_opt;
void select_qp_solver(int s);

#endif // SOLVER_SELECTOR
