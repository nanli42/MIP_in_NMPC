#ifndef NMPC4RacingCar_TRACK_H
#define NMPC4RacingCar_TRACK_H

#include "system_config.h"
#include "arc_length_spline.h"

namespace nmpc4rc {
  using json = nlohmann::json;
  namespace plt = matplotlibcpp;

  struct TrackPos {
    const Eigen::VectorXd X;
    const Eigen::VectorXd Y;

    const Eigen::VectorXd X_inner;
    const Eigen::VectorXd Y_inner;

    const Eigen::VectorXd X_outer;
    const Eigen::VectorXd Y_outer;
  };

  struct CurvInfo {
    const Eigen::VectorXd kappa;
    const Eigen::VectorXd length;
  };

  class Track {
    public:
      ArcLengthSpline sp;

      TrackPos loadTrack();
      void fitTrackIntoSpline(const TrackPos &track_xy);
  };

}

#endif // NMPC4RacingCar_TRACK_H
