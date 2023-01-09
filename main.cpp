#include "system_config.h"

#include "track.h"
#include "shadow_vehicle.h"

#include "MIP.h"
#include "solver_selector.hpp"

typedef std::vector<double> state_type;

#define GENERATE_LV_TRAJ 0
#define singRC 0
#define METHOD 1

int main(int argc, char const *argv[]) {

  using namespace nmpc4rc;

#if GENERATE_LV_TRAJ == 1
  Shadow_vehicle shadowVehicle;
  shadowVehicle.creation();
  shadowVehicle.save_to_file();
#else

  double start_s_EV, start_s_LV;
  int total_steps;
  double ey0,epsi0,vx0,vy0;
  double w0,d0,delta0;
  int opt;
  std::string s_pre = "";

  std::cout << "argc: " << argc << std::endl;
  if (argc==12) {
    start_s_EV = std::stod(argv[1]);
    start_s_LV = std::stod(argv[2]);
    total_steps = std::stod(argv[3]);

    ey0 = std::stod(argv[4]);
    epsi0 = std::stod(argv[5]);
    vx0  = std::stod(argv[6]);
    vy0 = std::stod(argv[7]);
    w0 = std::stod(argv[8]);
    d0 = std::stod(argv[9]);
    delta0 = std::stod(argv[10]);

    s_pre = argv[11];
    std::cout << "s_pre: " << s_pre << std::endl;
  } else {
    std::cout << "[Enter the initial position: (range: 0.0->18.0)] " << std::endl;
    std::cin >> start_s_EV;
    std::cout << "[Enter the initial position of LV (Leading Vehicle): (range: 0.0->18.0)] " << std::endl;
    std::cin >> start_s_LV;
    std::cout << "You have set the initial progress as " << start_s_EV << ".\n" << std::endl;
    std::cout << "[Enter the total simulation steps:] " << std::endl;
    std::cin >> total_steps;
    std::cout << "You have set the total_steps as " << total_steps << ".\n" << std::endl;
    std::cout << "[Enter the inital state {ey0,epsi0,vx0,vy0, w0, d0, delta0}: (seperated by space)] " << std::endl;
    std::cin >> ey0 >> epsi0 >> vx0 >> vy0 >> w0 >> d0 >> delta0;
    std::cout << "You have set the state {ey0,epsi0,vx0,vy0} as "
      << ey0 << ", " << epsi0 << ", " << vx0 << ", " << vy0 << ".\n" << std::endl;
    std::cout << "[Enter you choice for QP solver: 1-qpQASES 2-GUROBI] " << std::endl;
  }

  State x0 = {ey0,epsi0,vx0,vy0, w0, 0.0, d0, delta0,  start_s_EV};

  Track track;
  TrackPos track_xy = track.loadTrack();
  track.fitTrackIntoSpline(track_xy);

  std::string s = s_pre + "test-" + std::to_string((int)(ey0*100))
                    + "-LV_s_" + std::to_string((int)(start_s_LV*100));
  std::cout << "file_name: " << s << std::endl;
  char file_name[s.size() + 1];
	strcpy(file_name, s.c_str());

  MIP mpc;

  #if (METHOD==1)
    mpc.performMPC(track, x0, total_steps, start_s_LV, file_name);
  #endif

#endif

  return 0;
}
