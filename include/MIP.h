#ifndef NMPC4RacingCar_MIP_H
#define NMPC4RacingCar_MIP_H

#include "system_config.h"
#include "track.h"

namespace nmpc4rc {
  class MIP {
    public:
      std::vector<std::vector<State>> performMPC(Track track, State &state, int total_steps, double start_s_LV, const char* file_name);
      std::vector<std::vector<State>> performMPC_for_generate_LV_traj(Track track, State &state0, int total_steps);
  };
}

#endif // NMPC4RacingCar_MIP_H
