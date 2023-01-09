#include "track.h"

#define TRACK_OPT 2

namespace nmpc4rc {

  TrackPos Track::loadTrack() {

  #if (TRACK_OPT==1)
    std::ifstream trackInput("../param/track1.json");
  #elif (TRACK_OPT==2)
    std::ifstream trackInput("../param/track2.json");
  #endif

    json trackInJsonFormat;
    trackInput >> trackInJsonFormat;

    std::vector<double> X = trackInJsonFormat["X"];
    std::vector<double> Y = trackInJsonFormat["Y"];

    std::vector<double> X_inner = trackInJsonFormat["X_i"];
    std::vector<double> Y_inner = trackInJsonFormat["Y_i"];

    std::vector<double> X_outer = trackInJsonFormat["X_o"];
    std::vector<double> Y_outer = trackInJsonFormat["Y_o"];

    return
      {
        Eigen::Map<Eigen::VectorXd>(X.data(), X.size()),
        Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size()),
        Eigen::Map<Eigen::VectorXd>(X_inner.data(), X_inner.size()),
        Eigen::Map<Eigen::VectorXd>(Y_inner.data(), Y_inner.size()),
        Eigen::Map<Eigen::VectorXd>(X_outer.data(), X_outer.size()),
        Eigen::Map<Eigen::VectorXd>(Y_outer.data(), Y_outer.size())
      };

  }

  void Track::fitTrackIntoSpline(const TrackPos &track_xy) {
      sp.gen2DSpline(track_xy.X,track_xy.Y);
  }
}
