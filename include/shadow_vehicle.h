#ifndef NMPC4RacingCar_SHADOW_VEHICLE_H
#define NMPC4RacingCar_SHADOW_VEHICLE_H

#include "system_config.h"

namespace nmpc4rc {

  class Shadow_vehicle {
    public:
      std::vector<std::vector<State>> shadow_vehicle_trajectory;

      void creation();
      void save_to_file();
      void fit_and_generate_header_files();
      State get_state(double s, double t);

      State get_LV_state(double s, double t);
      std::vector<double> get_LV_linear_interpolation(double s, double delta_t);

      void test_header_file();

      State integrate_in_term_of_t(State state, double delta_t);
      State integrate_in_term_of_s(State state, double delta_s);

  };

}

#endif // NMPC4RacingCar_SHADOW_VEHICLE_H
