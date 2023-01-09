#include "shadow_vehicle.h"
#include "LV_traj.h"
#include "MIP.h"
#include "step_param.h"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#if (TRACK_OPT==1)
  #define creation_steps 200
  #define TOTAL_LEN 200
#elif (TRACK_OPT==2)
  #define creation_steps 360
  #define TOTAL_LEN 360
#endif

namespace nmpc4rc {

  //////////////////////////////////////////////////////////////////////////////
  ////////                        Boost integration                     ////////
  //////////////////////////////////////////////////////////////////////////////
  using namespace boost::numeric::odeint;
  double m = 0.041;double Iz = 27.8e-6;double lf = 0.029;double lr = 0.033;
  double Cm1=0.287;double Cm2=0.0545;double Cr0=0.0518;double Cr2=0.00035;
  double Br = 3.3852;double Cr = 1.2691;double Dr = 0.1737;
  double Bf = 2.579;double Cf = 1.2;double Df = 0.192;

  typedef std::vector<double> state_type;

  void obs_odes_in_term_of_s( const state_type &x , state_type &dxds , double s) {
    double dsdt = (x[2]*cos(x[1])-x[3]*sin(x[1])) / (1-x[0]*x[11]); //x[13]);
    double alphaF = -atan((x[4]*lf+x[3])/x[2]) + x[7];
    double alphaR = atan((x[4]*lr-x[3])/x[2]);

    double Frx = (Cm1-Cm2*x[2])*x[6] - Cr0 - Cr2*x[2]*x[2];
    double Fry = Dr*sin(Cr*atan(Br*alphaR));
    double Ffy = Df*sin(Cf*atan(Bf*alphaF));

    dxds[0] = 1/dsdt * (x[2]*sin(x[1]) + x[3]*cos(x[1]));
    dxds[1] = 1/dsdt * x[4] - x[11];
    dxds[2] = 1/dsdt * (x[4]*x[3] + 1/m*(Frx-Ffy*sin(x[7])));
    dxds[3] = 1/dsdt * (-x[4]*x[2] + 1/m*(Fry+Ffy*cos(x[7])));
    dxds[4] = 1/dsdt * (1/Iz * (lf*Ffy*cos(x[7])-lr*Fry));
    dxds[5] = 1/dsdt * 1;
    dxds[6] = 1/dsdt * x[9];
    dxds[7] = 1/dsdt * x[10]; 

    dxds[8] = 1;

    dxds[9] = 0;
    dxds[10] = 0;

    dxds[11] = 0;
  }
  State Shadow_vehicle::integrate_in_term_of_s(State state, double delta_s) {
    state_type x(12); // = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, // 9 states
              //   1.0, 0.0, // 2 ctrls
              //   0.0 // 1 curvature
              // }; // initial conditions // 9 states + 2 ctrls + 1 curvature

    x[0] = state.ey;
    x[1] = state.epsi;
    x[2] = state.vx;
    x[3] = state.vy;
    x[4] = state.w;
    x[5] = state.t;
    x[6] = state.d;
    x[7] = state.delta;
    x[8] = state.s;

    x[9] = state.dd;
    x[10] = state.ddelta;

    x[11] = state.curvature;

    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    size_t steps =  integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ) , // abs_err = 1.0e-10 , rel_err = 1.0e-6
                      obs_odes_in_term_of_s, x, 0.0, delta_s, 1e-6 );  //...start_time, end_time, default_dt...

    State result_state;
    result_state.ey = x[0];
    result_state.epsi = x[1];
    result_state.vx = x[2];
    result_state.vy = x[3];
    result_state.w = x[4];
    result_state.t = x[5];
    result_state.d = x[6];
    result_state.delta = x[7];
    result_state.s = x[8];

    return result_state;
  };

  void obs_odes( const state_type &x , state_type &dxdt , double t) {
    double dsdt = (x[2]*cos(x[1])-x[3]*sin(x[1])) / (1-x[0]*x[11]); //x[13]);
    double alphaF = -atan((x[4]*lf+x[3])/x[2]) + x[7];
    double alphaR = atan((x[4]*lr-x[3])/x[2]);

    double Frx = (Cm1-Cm2*x[2])*x[6] - Cr0 - Cr2*x[2]*x[2];
    double Fry = Dr*sin(Cr*atan(Br*alphaR));
    double Ffy = Df*sin(Cf*atan(Bf*alphaF));

    dxdt[0] = (x[2]*sin(x[1]) + x[3]*cos(x[1]));
    dxdt[1] = x[4] - x[11]*dsdt;
    dxdt[2] = (x[4]*x[3] + 1/m*(Frx-Ffy*sin(x[7])));
    dxdt[3] = (-x[4]*x[2] + 1/m*(Fry+Ffy*cos(x[7])));
    dxdt[4] = (1/Iz * (lf*Ffy*cos(x[7])-lr*Fry));
    dxdt[5] = 1;
    dxdt[6] = x[9];
    dxdt[7] = x[10];

    dxdt[8] = dsdt;

    dxdt[9] = 0;
    dxdt[10] = 0;

    dxdt[11] = 0;
  }
  State Shadow_vehicle::integrate_in_term_of_t(State state, double delta_t) {
    state_type x(12); // = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, // 9 states
              //   1.0, 0.0, // 2 ctrls
              //   0.0 // 1 curvature
              // }; // initial conditions // 9 states + 2 ctrls + 1 curvature

    x[0] = state.ey;
    x[1] = state.epsi;
    x[2] = state.vx;
    x[3] = state.vy;
    x[4] = state.w;
    x[5] = state.t;
    x[6] = state.d;
    x[7] = state.delta;
    x[8] = state.s;

    x[9] = state.dd;
    x[10] = state.ddelta;

    x[11] = state.curvature;

    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    size_t steps =  integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ) , // abs_err = 1.0e-10 , rel_err = 1.0e-6
                      obs_odes, x, 0.0, delta_t, 1e-6 );  //...start_time, end_time, default_dt...

    State result_state;
    result_state.ey = x[0];
    result_state.epsi = x[1];
    result_state.vx = x[2];
    result_state.vy = x[3];
    result_state.w = x[4];
    result_state.t = x[5];
    result_state.d = x[6];
    result_state.delta = x[7];
    result_state.s = x[8];

    return result_state;
  };
  //////////////////////////////////////////////////////////////////////////////

  void Shadow_vehicle::creation() {
    MIP mpc;

    // State x0 = {0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    State x0 = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Track track;
    TrackPos track_xy = track.loadTrack();
    track.fitTrackIntoSpline(track_xy);

    std::cout <<      "=========================="        << std::endl;
    std::cout << "[-- Generating shadow vehicle trajectory !--]" << std::endl;
    std::cout <<      "==========================\n"        << std::endl;
    shadow_vehicle_trajectory = mpc.performMPC_for_generate_LV_traj(track, x0, creation_steps);
  };

  FILE *fp_s;
  void Shadow_vehicle::save_to_file() {
    fp_s = fopen("../include/LV_traj_tmp.h", "w+");

    fprintf(fp_s,"#ifndef LV_TRAJ\n");
    fprintf(fp_s,"#define LV_TRAJ\n\n");

    fprintf(fp_s,"double state_LV[%d][%d][9] = {\n", creation_steps, Nstep);
    for(int i=0; i<creation_steps; i++) {
      fprintf(fp_s,"{\n");
      for(int j=0; j<Nstep; j++) {
        fprintf(fp_s,"{%f, %f, %f, %f, %f, %f, %f, %f, %f}, \n",
          shadow_vehicle_trajectory[i][j].ey,
          shadow_vehicle_trajectory[i][j].epsi,
          shadow_vehicle_trajectory[i][j].vx,
          shadow_vehicle_trajectory[i][j].vy,
          shadow_vehicle_trajectory[i][j].w,
          shadow_vehicle_trajectory[i][j].t,
          shadow_vehicle_trajectory[i][j].d,
          shadow_vehicle_trajectory[i][j].delta,
          shadow_vehicle_trajectory[i][j].s
        );
      }
      fprintf(fp_s,"},\n");
    }
    fprintf(fp_s,"};\n\n");

    fprintf(fp_s,"double ctrl_LV[%d][%d][2] = {\n", creation_steps, Nstep);
    for(int i=0; i<creation_steps; i++) {
      fprintf(fp_s,"{\n");
      for(int j=0; j<Nstep; j++) {
        fprintf(fp_s,"{%f, %f}, \n",
          shadow_vehicle_trajectory[i][j].dd,
          shadow_vehicle_trajectory[i][j].ddelta
        );
      };
      fprintf(fp_s,"},\n");
    }
    fprintf(fp_s,"};\n\n");

    fprintf(fp_s,"double curvature_LV[%d][%d][1] = {\n", creation_steps, Nstep);
    for(int i=0; i<creation_steps; i++) {
      fprintf(fp_s,"{\n");
      for(int j=0; j<Nstep; j++) {
        fprintf(fp_s,"{%f}, \n",
          shadow_vehicle_trajectory[i][j].curvature
        );
      };
      fprintf(fp_s,"},\n");
    }
    fprintf(fp_s,"};\n\n");

    fprintf(fp_s,"#endif //LV_TRAJ\n");
    fclose(fp_s);
  };

  std::string exec(const char* cmd) {
      char buffer[128];
      std::string result = "";
      FILE* pipe = popen(cmd, "r");
      if (!pipe) throw std::runtime_error("popen() failed!");
      try {
          while (fgets(buffer, sizeof buffer, pipe) != NULL) {
              result += buffer;
          }
      } catch (...) {
          pclose(pipe);
          throw;
      }
      pclose(pipe);
      return result;
  }

  int find_LV_idx_start(double s) {
    double ss = s;
    if (ss>=TOTAL_LEN*step_s) {
      std::cout << "Exceed MAX length of LV's trajectory!" << std::endl;
      exit(1);
    }

    for (int i=0; i<TOTAL_LEN; i++) {
      double tmp1 = state_LV[i][0][8];
      double tmp2 = state_LV[i+1][0][8];
      if ( (ss>=tmp1)&&(ss<tmp2) ) return i;
    }
  return 0;
  }

  double find_LV_init_t(int idx_start, double s) {
    // using linear interpolation
    double k_s = (state_LV[idx_start+1][0][8]-state_LV[idx_start][0][8])/(state_LV[idx_start+1][0][5]-state_LV[idx_start][0][5]);
    return (s-state_LV[idx_start][0][8])/k_s;
  };

  int find_LV_idx_end(int idx_start, double delta_t, double s) {
    int idx = idx_start+1;
    double init_t = find_LV_init_t(idx_start, s);
    while ( (state_LV[idx][0][5]-state_LV[idx_start][0][5]-init_t)<delta_t && (idx<TOTAL_LEN) )
      idx++;
    return idx;
  }

  // state_LV_0
  State Shadow_vehicle::get_LV_state(double s, double delta_t) {
    int idx_start = find_LV_idx_start(s);
    int idx_end = find_LV_idx_end(idx_start, delta_t, s);
    State state = {
      state_LV[idx_end-1][0][0],
      state_LV[idx_end-1][0][1],
      state_LV[idx_end-1][0][2],
      state_LV[idx_end-1][0][3],
      state_LV[idx_end-1][0][4],
      state_LV[idx_end-1][0][5],
      state_LV[idx_end-1][0][6],
      state_LV[idx_end-1][0][7],
      state_LV[idx_end-1][0][8],
      0.0,
      0.0,
      ctrl_LV[idx_end-1][0][0],
      ctrl_LV[idx_end-1][0][1],
      curvature_LV[idx_end-1][0][0]
    };

    double init_t = find_LV_init_t(idx_start, s);
    return integrate_in_term_of_t(state, delta_t-(state_LV[idx_end-1][0][5]-state_LV[idx_start][0][5]-init_t));
  };

  std::vector<double> Shadow_vehicle::get_LV_linear_interpolation(double s, double delta_t) {
    int idx_start = find_LV_idx_start(s);
    int idx_end = find_LV_idx_end(idx_start, delta_t, s);

    double ds = state_LV[idx_end][0][8]-state_LV[idx_end-1][0][8];
    double dey = state_LV[idx_end][0][0]-state_LV[idx_end-1][0][0];
    double dt = state_LV[idx_end][0][5]-state_LV[idx_end-1][0][5];

    double init_t = find_LV_init_t(idx_start, s);
    double current_dt = delta_t - (state_LV[idx_end-1][0][5]-state_LV[idx_start][0][5]-init_t);

    std::vector<double> res;
    res.push_back(ds/dt);
    res.push_back(state_LV[idx_end-1][0][8] + ds/dt * current_dt);
    res.push_back(dey/dt);
    res.push_back(state_LV[idx_end-1][0][0] + dey/dt * current_dt);

    return res;
  };

}
