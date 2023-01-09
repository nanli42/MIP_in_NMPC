#ifndef NMPC4RacingCar_AUX_FUNC_H
#define NMPC4RacingCar_AUX_FUNC_H
#include "track.h"
#include "plot.h"
#include "solver_selector.hpp"
#include "shadow_vehicle.h"
#include "shadow_vehicle.h"

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

namespace nmpc4rc {

  double calculate_curvature(ArcLengthSpline sp, double s);
  void verify_EV_LV_relative_position(State state_EV, State state_LV, int step);
  bool check_collision(ArcLengthSpline sp, int silence, int N);
  void replay_trajectory(Track track, std::vector<std::vector<State>> state_series_EV,
    std::vector<std::vector<State>> state_series_LV,
    int Nstep, int Nobs, double init_pos, double start_s_LV);
  std::vector<State> stock_state_and_plot(Plotting plot, State state_LV, Shadow_vehicle shadow_vehicle, Track track, int Nstep, int Nobs, FILE* fp_state_record);
  void prep_init_states(State state0, Track track, double init_pos, int Nstep, double step_s);
  bool check_overtake(ArcLengthSpline sp, double s_LV);

  std::vector<double> get_pos_x_y(ArcLengthSpline sp, double s, double ey);
}

#endif // NMPC4RacingCar_LOROFO_H
