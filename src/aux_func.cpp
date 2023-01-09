#include "aux_func.h"

namespace nmpc4rc {

  std::vector<double> get_pos_x_y(ArcLengthSpline sp, double s, double ey) {
    const Eigen::Vector2d dpos_ref = sp.getDerivative(s);
    const double dx_ref = dpos_ref(0);
    const double dy_ref = dpos_ref(1);
    const double theta_ref = atan2(dy_ref,dx_ref);

    Eigen::Vector2d pos = sp.getPostion(s);
    double x = pos(0) - ey * std::sin( theta_ref );
    double y = pos(1) + ey * std::cos( theta_ref );

    std::vector<double> xy;
    xy.push_back(x);
    xy.push_back(y);
    xy.push_back(theta_ref);
    return xy;
  }

  double calculate_curvature(ArcLengthSpline sp, double s) {
    double dx = sp.getDerivative(s)(0);
    double dy = sp.getDerivative(s)(1);
    double ddx = sp.getSecondDerivative(s)(0);
    double ddy = sp.getSecondDerivative(s)(1);
    double k = (dx*ddy-dy*ddx)/std::pow((dx*dx+dy*dy), 3/2);

    if (abs(k)<1e-5) {
      return 1e-5;
    }

    if (abs(k)>5) {
      printf("original ks: %f\n", k);
      if (k < 0) k = -5;
      else k = 5;
    }

    return k;
  };

  void verify_EV_LV_relative_position(State state_EV, State state_LV, int step) {
    double delta_s_0 = 0.029+0.033;
    double delta_ey_0 = 2*0.015;
    double delta_s = (delta_s_0 * cos(state_EV.epsi) * (1/(1-state_EV.ey*state_EV.curvature)))/2;
    double delta_ey = (delta_ey_0 * cos(state_EV.epsi) + delta_s_0 * sin(abs(state_EV.epsi)))/2;

    double tol = 1e-3;
    std::cout << ">>> step: " << step << " | ";
    if (
      (state_EV.s + delta_s*2 <= state_LV.s + tol) &&
      (state_LV.ey - tol - (state_LV.s-state_EV.s)/delta_s_0*delta_ey_0 <= state_EV.ey) &&
      (state_EV.ey <= state_LV.ey + (state_LV.s-state_EV.s)/delta_s_0*delta_ey_0 + tol)
    )
      std::cout << "S- ";

    if (
      (state_EV.s - delta_s*2 >= state_LV.s - tol) &&
      (state_LV.ey - tol - (-state_LV.s+state_EV.s)/delta_s_0*delta_ey_0 <= state_EV.ey) &&
      (state_EV.ey <= state_LV.ey + (-state_LV.s+state_EV.s)/delta_s_0*delta_ey_0 + tol)
    )
      std::cout << "S+ ";

    if (
      (state_EV.ey + delta_ey*2 <= state_LV.ey + tol) &&
      (state_LV.s - tol - (state_LV.ey-state_EV.ey)/delta_ey_0*delta_s_0 <= state_EV.s) &&
      (state_EV.s <= state_LV.s + (state_LV.ey-state_EV.ey)/delta_ey_0*delta_s_0 + tol)
    )
      std::cout << "Ey- ";

    if (
      (state_EV.ey - delta_ey*2 >= state_LV.ey - tol) &&
      (state_LV.s - tol - (-state_LV.ey+state_EV.ey)/delta_ey_0*delta_s_0 <= state_EV.s) &&
      (state_EV.s <= state_LV.s + (-state_LV.ey+state_EV.ey)/delta_ey_0*delta_s_0 + tol)
    )
      std::cout << "Ey+ ";

    std::cout << std::endl;


    std::cout
      << "state_EV.s: "<< state_EV.s
      << ", state_LV.s: "<< state_LV.s
      << ", delta_s: "<< delta_s
      << std::endl
      << "state_EV.ey: "<< state_EV.ey
      << ", state_LV.ey: "<< state_LV.ey
      << ", delta_ey: "<< delta_ey
    << std::endl;
  };


  bool check_collision(ArcLengthSpline sp, int silence, int N) {
    double delta_s_0 = 0.029+0.033;
    double delta_ey_0 = 2*0.015;
    bool flag=true;

    if (silence==0) std::cout << ">>>>>> check collision result <<<<<<" << std::endl;

    for (int i=1; i<N; i++) {
      double ey = acadoVariables.x[i*9+0];
      double ephi = acadoVariables.x[i*9+1];
      double s = acadoVariables.x[i*9+8];
      double ks = calculate_curvature(sp, s);

      double delta_s = (delta_s_0 * cos(ephi) * (1/(1-ey*ks)))/2;
      double delta_ey_pos = (delta_ey_0 * cos(ephi) + abs(delta_s_0 * sin(ephi)))/2;

      double s_LV = acadoVariables.od[i*7+3]
        + (acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[i*7+1]) * acadoVariables.od[i*7+2];
      double ey_LV = acadoVariables.od[i*7+5]
        + (acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[i*7+1]) * acadoVariables.od[i*7+4];

      if (silence==0) {
        std::cout << ">>> step " << i << ": ";
        if (s+delta_s*2<=s_LV) std::cout << "S- ";
        if (s-delta_s*2>=s_LV) std::cout <<  "S+ ";
        if (ey+delta_ey_pos*2<=ey_LV) std::cout << "RO ";
        if (ey-delta_ey_pos*2>=ey_LV) std::cout << "LO ";
        std::cout << std::endl;

        std::cout
          << "state_EV.s: "<< s
          << ", state_LV.s: "<< s_LV
          << ", delta_s*2: "<< delta_s*2
          << std::endl
          << "state_EV.ey: "<< ey
          << ", state_LV.ey: "<< ey_LV
          << ", delta_ey*2: "<< delta_ey_pos*2
        << std::endl;
      }

      std::cout << s+delta_s*2-s_LV << std::endl;
      std::cout << ey+delta_ey_pos*2-ey_LV << std::endl;

      double tol = 1e-4;
      if (s+delta_s*2<=s_LV+tol || s-delta_s*2>=s_LV-tol || ey+delta_ey_pos*2<=ey_LV+tol || ey-delta_ey_pos*2>=ey_LV-tol)
        continue;
      else flag=false;
    }
    return flag;
  }

  void replay_trajectory(Track track, std::vector<std::vector<State>> state_series_EV,
    std::vector<std::vector<State>> state_series_LV,
    int Nstep, int Nobs, double init_pos, double start_s_LV) {

    Plotting plot;

    // User Interface to reshow the i th step
    std::cout << "You could now close the online preview window and replay specific step. \n" << std::endl;
    plot.finishPlot();

    std::cout << "[Which step you want to replay?] \n(enter an integer to select step; enter \"stop\" to exit) " << std::endl;
    std::string input;
    std::cin>>input;

    while (input.compare("stop")!=0){
      std::cout << "You selected the step " << input << ".\n" << std::endl;

      plot.plotTrack(track.loadTrack());
      std::vector<State> states_EV;
      for(int i=0; i<Nstep+1; i++) {
        // Plot EV's predicted traj at step "input"
        if (i%5==0) plot.plotVehicle_inferred_in_xy_coordinate(track, state_series_EV[std::stoi(input)][i], "b.--", "bo");
        else plot.plotVehicle_inferred_in_xy_coordinate(track, state_series_EV[std::stoi(input)][i], "b..", "b.");

        states_EV.push_back(state_series_EV[std::stoi(input)][i]);
      }
      plot.plot_line(track, states_EV, "b-");

      std::vector<State> states_LV;
      for(int i=0; i<Nobs+1; i++) {
        State state = state_series_LV[std::stoi(input)][i];

        if (i%5==0) plot.plotVehicle_inferred_in_xy_coordinate(track, state, "g.--", "go");
        else plot.plotVehicle_inferred_in_xy_coordinate(track, state, "g..", "g.");

        plot.writeTitle("EV start from progess: " + std::to_string(init_pos).substr(0,3)
          + " [m]; LV start from progess: " + std::to_string(start_s_LV).substr(0,4) + " [m]; STEP: " + input);

        verify_EV_LV_relative_position(state_series_EV[std::stoi(input)][i], state_series_LV[std::stoi(input)][i], i);

        states_LV.push_back(state);
      }
      plot.plot_line(track, states_LV, "g-");

      plot.finishPlot();

      std::cout << "Enter another step value to have an another try (\"stop\" to stop): " << std::endl;
      std::cin>>input;
      if (input.compare("\x03")==0) break;
    }
  }

  std::vector<State> stock_state_and_plot(Plotting plot, State state_LV,
    Shadow_vehicle shadow_vehicle, Track track, int Nstep, int Nobs,
    FILE* fp_state_record) {

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
        acadoVariables.od[i*7 + 0 ],
      };
      for (int j=0; j<9; j++) fprintf(fp_state_record,"%f ", acadoVariables.x[ i*9 + j ]);
      for (int j=0; j<2; j++) fprintf(fp_state_record,"%f ", acadoVariables.u[ i*2 + j ]);
      fprintf(fp_state_record,"\n");

      //plot.plotVehicle_inferred_in_xy_coordinate(track, state);
      if (i%5==0) plot.plotVehicle_inferred_in_xy_coordinate(track, state, "b-", "bx");
      else plot.plotVehicle_inferred_in_xy_coordinate(track, state, "b..", "b.");

      state_vec.push_back(state);
    }
    plot.plot_line(track, state_vec, "b-");

    std::vector<State> state_vec_LV;
    for(int i=0; i<Nobs+1; i++) {
        State state = shadow_vehicle.get_LV_state(state_LV.s, state_vec[i].t-state_vec[0].t);
        state.s = acadoVariables.od[i*7+3]
          + (acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[i*7+1]) * acadoVariables.od[i*7+2];
        state.ey = acadoVariables.od[i*7+5]
          + (acadoVariables.x[ i*9 + 5 ]-acadoVariables.od[i*7+1]) * acadoVariables.od[i*7+4];

        //plot.plotVehicle_inferred_in_xy_coordinate(track, state, "g.-", "gx");
        if (i%5==0) plot.plotVehicle_inferred_in_xy_coordinate(track, state, "g-", "gx");
        else plot.plotVehicle_inferred_in_xy_coordinate(track, state, "g..", "g.");

        state_vec_LV.push_back(state);
    }
    plot.plot_line(track, state_vec_LV, "g-");

    return state_vec;
  }

  void prep_init_states(State state0, Track track, double init_pos, int Nstep, double step_s) {
    int i;

    /* Initialize the states and controls. */
    for (i = 0; i < Nstep+1; ++i) {
      /* states 0-8: */
      /* ey, ephi, vx, vy, w, t, d, delta, s */
      acadoVariables.x[ i*9 + 0 ] = 0.0;
      acadoVariables.x[ i*9 + 1 ] = 0.0;
      acadoVariables.x[ i*9 + 2 ] = 1.0;
      acadoVariables.x[ i*9 + 3 ] = 0.0;
      acadoVariables.x[ i*9 + 4 ] = 0.0;
      acadoVariables.x[ i*9 + 5 ] = 0.0;
      acadoVariables.x[ i*9 + 6 ] =1.0; // 1.0;// 0.0; - N15 //1.0; - N30
      double k = calculate_curvature(track.sp, init_pos+i*step_s);
      acadoVariables.x[ i*9 + 7 ] = atan(k*0.06);
      acadoVariables.x[ i*9 + 8 ] = init_pos+i*step_s;
    }

    for (i = 0; i < Nstep; ++i) {
      acadoVariables.u[ i*2 ] = 0;
      acadoVariables.u[ i*2+1 ] = 0;
    }
    for (i = 0; i < Nstep+1; ++i) {
      acadoVariables.od[i*7 + 0 ] = calculate_curvature(track.sp, init_pos+i*step_s);
      acadoVariables.od[i*7 + 1 ] = 0.0;
      acadoVariables.od[i*7 + 2 ] = 0.0;
      acadoVariables.od[i*7 + 3 ] = 0.0;
      acadoVariables.od[i*7 + 4 ] = 0.0;
      acadoVariables.od[i*7 + 5 ] = 0.0;
      acadoVariables.od[i*7 + 6 ] = 1e-6;
    }

    /* MPC: initialize the current state feedback. */
    /* states 0-8: */
    /* ey, ephi, vx, vy, w, t, d, delta, s */
    acadoVariables.x0[ 0 ] = state0.ey;
    acadoVariables.x0[ 1 ] = state0.epsi;
    acadoVariables.x0[ 2 ] = state0.vx;
    acadoVariables.x0[ 3 ] = state0.vy;
    acadoVariables.x0[ 4 ] = state0.w;
    acadoVariables.x0[ 5 ] = 0.0;
    acadoVariables.x0[ 6 ] = state0.d; //1.0;
  	acadoVariables.x0[ 7 ] = state0.delta;
  	acadoVariables.x0[ 8 ] = init_pos;

    printf("===============================================\n");
    printf("<<< Prepare init states with QPOASES solver >>>\n");
    printf("===============================================\n");

    double prec = 1.0;
    int sqp = 0;
    while ( (prec>1e-6)&&(sqp<20) ) {
      printf("Iter %d, KKT Tolerance = %.3e\n", sqp, prec);

      acado_preparationStep();

      acado_feedbackStep();

      sqp ++;
      prec = acado_getKKT();
      printf("Iter %d, KKT Tolerance = %.3e\n", sqp, prec);
    }
    acado_printDifferentialVariables();

    if (std::isnan(acado_getKKT())) {
      printf("KKT Tolerance = %.3e\n", acado_getKKT());
      printf("No feasible solution. Try another initialization setting!\n");
      exit(1);
    }
  }

  bool check_overtake(ArcLengthSpline sp, double s_LV) {
    double delta_s_0 = 0.029+0.033;
    double delta_ey_0 = 2*0.015;
    int i = 0;

    double ey = acadoVariables.x[i*9+0];
    double ephi = acadoVariables.x[i*9+1];
    double s = acadoVariables.x[i*9+8];
    double ks = calculate_curvature(sp, s);

    double delta_s = (delta_s_0 * cos(ephi)  * (1/(1-ey*ks)))/2;

    if (s-delta_s*2           *2             >=s_LV) {  // two times far away, avoid the close overtaking
      std::cout << "CHECK! s: " << s << std::endl;
      std::cout << "CHECK! s_LV: " << s_LV << std::endl;
      std::cout << "CHECK! delta_s: " << delta_s << std::endl;
      return true;
    }
    else return false;
  }
}
