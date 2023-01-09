#include "MIP.h"
#include "plot.h"
#include "shadow_vehicle.h"
#include "aux_func.h"
#include "step_param.h"
#include "solver_selector.hpp"

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

namespace nmpc4rc {

std::vector<std::vector<State>> MIP::performMPC(Track track, State &state0, int total_steps, double start_s_LV, const char* file_name) {

  FILE *fp_state_record;
  std::string fn = std::string(file_name)+"-state_record-h2h"+".txt";
  const char *cstr = fn.c_str();
  fp_state_record = fopen(cstr, "w+");

	int    i, j, step;

  // Init position for EV and a fixed step-progress
  double init_pos = state0.s;

  // Get LV's state at time 0.0
  Shadow_vehicle shadow_vehicle;
  State state_LV = shadow_vehicle.get_LV_state(start_s_LV, 0);
  state_LV.s = start_s_LV;

	/* Initialize the solver. */
	acado_initializeSolver();

  select_qp_solver(OPT_QPOASES_QP_SOLVER);
  prep_init_states(state0, track, init_pos, Nstep, step_s);
  select_qp_solver(OPT_GUROBI_QP_SOLVER);

  acado_printDifferentialVariables();

  Plotting plot;
  std::vector<std::vector<State>> state_series_EV;
  std::vector<std::vector<State>> state_series_LV;


  acado_timer t;
  acado_tic( &t );

  // Step for EV
  for(step = 0; step < total_steps; ++step) {

    std::cout <<      "=========================="        << std::endl;
    std::cout << "[-----STEP "<< step << " ------------]" << std::endl;
    std::cout <<      "=========================="        << std::endl;

    double prec = 1.0;
    int sqp = 0;

    acadoVariables.od[ 0*7 + 1 ] = acadoVariables.x0[ 5 ];
    acadoVariables.od[ 0*7 + 2 ] = 0.0;
    acadoVariables.od[ 0*7 + 3 ] = state_LV.s;
    acadoVariables.od[ 0*7 + 4 ] = 0.0;
    acadoVariables.od[ 0*7 + 5 ] = state_LV.ey;
    acadoVariables.od[ 0*7 + 6 ] = calculate_curvature(track.sp, state_LV.s);

    for(int i=1; i<Nstep+1; i++) {
      acadoVariables.od[ i*7 + 1 ] = acadoVariables.x[ i*9 + 5 ];

      std::vector<double> lp = shadow_vehicle.get_LV_linear_interpolation(state_LV.s, acadoVariables.x[ i*9 + 5 ]-acadoVariables.x[ 0*9 + 5 ]);

      acadoVariables.od[ i*7 + 2 ] = lp[0];
      acadoVariables.od[ i*7 + 3 ] = lp[1];
      acadoVariables.od[ i*7 + 4 ] = lp[2];
      acadoVariables.od[ i*7 + 5 ] = lp[3];
      acadoVariables.od[ i*7 + 6 ] = calculate_curvature(track.sp, lp[1]);
    }

    double last_kkt = 1e+6;

    // Iterate until the precision is enough or it reaches the max iteration limit
    while ( ( (prec>1e-4)||prec==0.0 ) &&(sqp<20)  ) {
      /* Prepare step */
      acado_preparationStep();
      /* Feedback step */
      acado_feedbackStep();

      sqp ++;
      prec = acado_getKKT();
      printf("KKT Tolerance = %.3e\n", acado_getKKT());

      if (std::isnan(prec)||prec>1e+1) {
        printf("KKT Tolerance = %.3e\n", acado_getKKT());
        printf("No feasible solution. Try sol!\n");
        break;
      }
    }

    printf("KKT Tolerance = %.3e\n", acado_getKKT());
    if (std::isnan(acado_getKKT())) {
      printf("No feasible solution. Try another initialization setting!\n");
      exit(1);
    }

    std::vector<State> state_vec_LV;
    for(int i=0; i<Nobs+1; i++) {
        State state;
        state.s = acadoVariables.od[ i*7 + 3 ]
          + acadoVariables.od[ i*7 + 2 ]*(acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[ i*7 + 1 ]);
        state.ey = acadoVariables.od[ i*7 + 5 ]
          + acadoVariables.od[ i*7 + 4 ]*(acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[ i*7 + 1 ]);
        state_vec_LV.push_back(state);
    }
    state_series_LV.push_back(state_vec_LV);
    std::vector<State> state_vec_EV = stock_state_and_plot(plot, state_LV, shadow_vehicle, track, Nstep, Nobs, fp_state_record);

    // Update the LV's pos for next step
    double tmp_LV_s = state_LV.s;

    state_LV.s = acadoVariables.od[ 1*7 + 3 ] + acadoVariables.od[ 1*7 + 2 ]*(state_vec_EV[1].t-acadoVariables.od[ 1*7 + 1 ]);
    state_LV.ey = acadoVariables.od[ 1*7 + 5 ] + acadoVariables.od[ 1*7 + 4 ]*(state_vec_EV[1].t-acadoVariables.od[ 1*7 + 1 ]);

    // Save up EV's state and show plot
    state_series_EV.push_back(state_vec_EV);

    plot.drawPlot();

    acado_printDifferentialVariables();
    acado_printControlVariables();

    // Prepare EV's curvature information for next step
    for (i = 0; i < Nstep+1; ++i)
      acadoVariables.od[ i*7 + 0 ] = calculate_curvature(track.sp,init_pos+(step+i+1)*step_s);
    // Update EV's intial state for next step
    for (i = 0; i < 9; ++i)
      acadoVariables.x0[ i ] = acadoVariables.x[ 9 + i ];

    // Update EV's state and control for next step by shifting
    acado_shiftStates(2, 0, 0);
    acado_shiftControls(0);

    std::vector<double> res1 = get_pos_x_y(track.sp, acadoVariables.x0[8], acadoVariables.x0[0] );
    std::cout << "x_EV, y_EV: " << res1[0] << ", " << res1[1] << std::endl;
    std::vector<double> res2  = get_pos_x_y(track.sp, state_LV.s, state_LV.ey);
    std::cout << "x_LV, y_LV: " << res2[0] << ", " << res2[1] << std::endl;
    std::cout << "state_LV.s: " << state_LV.s << std::endl;

    std::cout << sqrt((res1[0]-res2[0])*(res1[0]-res2[0])+(res1[1]-res2[1])*(res1[1]-res2[1])) << std::endl;

    if (sqrt((res1[0]-res2[0])*(res1[0]-res2[0])+(res1[1]-res2[1])*(res1[1]-res2[1]))<0.03444*2) {
      std::cout << "\n\n\n\n\nTOO CLOSE!!!\n\n\n\n\n" << std::endl;
      exit(1);
    }

  }

  real_t te = acado_toc( &t );
  std::cout << "\n\n\n\n\n===== Timing: " << te / total_steps * 1e3 << std::endl;

  replay_trajectory(track, state_series_EV, state_series_LV, Nstep, Nobs, init_pos, start_s_LV);

  fclose(fp_state_record);

  return state_series_EV;
}

std::vector<std::vector<State>> MIP::performMPC_for_generate_LV_traj(Track track, State &state0, int total_steps) {
  select_qp_solver(OPT_QPOASES_QP_SOLVER);

  int    i, j, step;
  acado_timer t;

  FILE *fp_time_record;
  fp_time_record = fopen("time_record.txt", "w+");

  double init_pos = state0.s;

  /* Initialize the solver. */
  acado_initializeSolver();

  /* Initialize the states and controls. */
  for (i = 0; i < Nstep; ++i) {
    /* states 0-8: */
    /* ey, ephi, vx, vy, w, t, d, delta, s */
    acadoVariables.x[ i*9 + 0 ] = 0.0;
    acadoVariables.x[ i*9 + 1 ] = 0.0;
    acadoVariables.x[ i*9 + 2 ] = 0.8; // 1.2 track1 ; 1.5 tack2
    acadoVariables.x[ i*9 + 3 ] = 0.0;
    acadoVariables.x[ i*9 + 4 ] = 0.0;
    acadoVariables.x[ i*9 + 5 ] = 0.0;

    double k = calculate_curvature(track.sp, init_pos+i*step_s);
    acadoVariables.x[ i*9 + 6 ] = 1.0;
    acadoVariables.x[ i*9 + 7 ] = atan(k*0.06);

    acadoVariables.x[ i*9 + 8 ] = init_pos+i*step_s;
  }

  for (i = 0; i < Nstep; ++i) {
    acadoVariables.u[ i*2 ] = 0;
    acadoVariables.u[ i*2+1 ] = 0;
  }
  for (i = 0; i < Nstep; ++i) {
    acadoVariables.od[ i*1 + 0 ] = calculate_curvature(track.sp, init_pos+i*step_s);
  }

  /* MPC: initialize the current state feedback. */
  /* states 0-8: */
  /* ey, ephi, vx, vy, w, t, d, delta, s */
  acadoVariables.x0[ 0 ] = state0.ey;
  acadoVariables.x0[ 1 ] = state0.epsi;
  acadoVariables.x0[ 2 ] = state0.vx;
  acadoVariables.x0[ 3 ] = state0.vy;
  acadoVariables.x0[ 4 ] = 0.0;
  acadoVariables.x0[ 5 ] = 0.0;
  acadoVariables.x0[ 6 ] = 0.0;
  acadoVariables.x0[ 7 ] = 0.0;
  acadoVariables.x0[ 8 ] = init_pos;

  for (i = 0; i < 9; ++i) {
    std::cout << "acadoVariables.x0: " << acadoVariables.x0[ i ] << std::endl;
  }

  Plotting plot;
  std::vector<std::vector<State>> state_series_EV;

  for(step = 0; step < total_steps; ++step) {
  //for(step = 0; step < 78; ++step) {

    std::cout <<      "=========================="        << std::endl;
    std::cout << "[-----STEP "<< step << " ------------]" << std::endl;
    std::cout <<      "=========================="        << std::endl;

    double prec = 1.0;
    int sqp = 0;

    acado_tic( &t );

    while ( (prec>1e-6)&&(sqp<20) ) {
      /* Prepare step */
      acado_preparationStep();
      /* Feedback step */
      acado_feedbackStep();

      sqp ++;
      prec = acado_getKKT();
    }

    real_t te = acado_toc( &t );
    fprintf(fp_time_record,"%d %d %.3g %.3g %.3g\n",   step, sqp, 1e6 * te, 1e6 * te / sqp, prec);

    printf("\nTotal time:   %.3g microseconds\n", 1e6 * te );
    printf("Average time:   %.3g microseconds\n\n", 1e6 * te / sqp);
    if (std::isnan(acado_getKKT())) {
      printf("KKT Tolerance = %.3e\n", acado_getKKT());
      printf("No feasible solution. Try another initialization setting!\n");
      exit(1);
    }

    plot.plotTrack(track.loadTrack());

    std::vector<State> state_vec;
    for(int i=0; i<Nstep+1; i++) {
      State state = {
        acadoVariables.x[ i*9 + 0 ],
        acadoVariables.x[ i*9 + 1 ],
        acadoVariables.x[ i*9 + 2 ],
        acadoVariables.x[ i*9 + 3 ],
        acadoVariables.x[ i*9 + 4 ],
        acadoVariables.x[ i*9 + 5 ],
        acadoVariables.x[ i*9 + 6 ],
        acadoVariables.x[ i*9 + 7 ],
        acadoVariables.x[ i*9 + 8 ],
        0.0,
        0.0,
        acadoVariables.u[ i*2 + 0 ],
        acadoVariables.u[ i*2 + 1 ],
        acadoVariables.od[ i*1 + 0 ],
      };
      plot.plotVehicle_inferred_in_xy_coordinate(track, state);
      state_vec.push_back(state);
    }

    state_series_EV.push_back(state_vec);
    plot.drawPlot();

    acado_printDifferentialVariables();
    acado_printControlVariables();

    for (i = 0; i < Nstep+1; ++i) {
      acadoVariables.od[ i*1 + 0 ] = calculate_curvature(track.sp, init_pos+(step+i+1)*step_s);
    }
    for (i = 0; i < 9; ++i) {
      acadoVariables.x0[ i ] = acadoVariables.x[ 9 + i ];
    }
    acado_shiftStates(2, 0, 0);
    acado_shiftControls(0);
  }

  fclose(fp_time_record);

  // User Interface to reshow the i th step
  std::cout << "You could now close the online preview window and replay specific step. \n" << std::endl;
  plot.finishPlot();

  std::cout << "[Which step you want to replay?] \n(enter an integer to select step; enter \"stop\" to exit) " << std::endl;
  std::string input;
  std::cin>>input;

  while (input.compare("stop")!=0){
    std::cout << "You selected the step " << input << ".\n" << std::endl;

    plot.plotTrack(track.loadTrack());
    for(int i=0; i<Nstep+1; i++) {
      // Plot EV's predicted traj at step "input"
      plot.plotVehicle_inferred_in_xy_coordinate(track, state_series_EV[std::stoi(input)][i]);
    }

    plot.finishPlot();

    std::cout << "Enter another step value to have an another try (\"stop\" to stop): " << std::endl;
    std::cin>>input;
    if (input.compare("\x03")==0) break;
  }
  return state_series_EV;
}

}
