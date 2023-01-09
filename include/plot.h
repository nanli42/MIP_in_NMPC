#ifndef NMPC4RacingCar_PLOT_H
#define NMPC4RacingCar_PLOT_H

#include "system_config.h"
#include "track.h"

namespace nmpc4rc {
  namespace plt = matplotlibcpp;

  class Plotting {
    public:
      void plotTrack(const TrackPos &track_xy) const;
      void plotVehicle_inferred_in_xy_coordinate(const Track track, const State &state, std::string boxcolor="b.-", std::string centercolor="bx", int f=0) const;
      void plotVehicle_inferred_in_cuvilinear_coordinate(const Track track, const State &state) const;
      void finishPlot();
      void drawPlot();
      void writeTitle(std::string s);
      void savePlot(std::string s);

      void plot_line(Track track, std::vector<State> states, std::string col);
  };
}

#endif // NMPC4RacingCar_PLOT_H
