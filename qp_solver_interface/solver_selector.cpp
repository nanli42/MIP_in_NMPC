#include "solver_selector.hpp"

extern "C"
{
#include "acado_common.h"
}

int solver_opt;
void select_qp_solver(int s) {
  solver_opt = s;
}

int acado_solve() {
  if (solver_opt == OPT_QPOASES_QP_SOLVER) return acado_solve_qpoases();
  if (solver_opt == OPT_GUROBI_QP_SOLVER) return solver_gurobi();
};
