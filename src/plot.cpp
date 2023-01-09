#include "plot.h"
#include "arc_length_spline.h"

#define TRACK_OPT 2

namespace nmpc4rc {

  void Plotting::plotTrack(const TrackPos &track_xy) const {
    std::vector<double> plot_xc(track_xy.X.data(),track_xy.X.data() + track_xy.X.size());
    std::vector<double> plot_yc(track_xy.Y.data(),track_xy.Y.data() + track_xy.Y.size());

    std::vector<double> plot_xi(track_xy.X_inner.data(),track_xy.X_inner.data() + track_xy.X_inner.size());
    std::vector<double> plot_yi(track_xy.Y_inner.data(),track_xy.Y_inner.data() + track_xy.Y_inner.size());
    std::vector<double> plot_xo(track_xy.X_outer.data(),track_xy.X_outer.data() + track_xy.X_outer.size());
    std::vector<double> plot_yo(track_xy.Y_outer.data(),track_xy.Y_outer.data() + track_xy.Y_outer.size());

    plt::clf();
    plt::figure(1);
    plt::plot(plot_xc,plot_yc,"r--");
    plt::plot(plot_xi,plot_yi,"k-");
    plt::plot(plot_xo,plot_yo,"k-");

    plt::axis("equal");
  #if (TRACK_OPT==1)
    plt::xlim(-1.0,1.5);
    plt::ylim(-2.0,0.5);
  #elif (TRACK_OPT==2)
    plt::xlim(-1.5,2.0);
    plt::ylim(-2.0,2.0);
  #endif
  }

  std::vector<double> get_x_y(ArcLengthSpline sp, double s, double ey) {
    const Eigen::Vector2d dpos_ref = sp.getDerivative(s);
    const double dx_ref = dpos_ref(0);
    const double dy_ref = dpos_ref(1);
    const double theta_ref = atan2(dy_ref,dx_ref);

    //std::cout << "theta_ref: " << theta_ref/3.1415926*180 << std::endl;

    Eigen::Vector2d pos = sp.getPostion(s);
    double x = pos(0) - ey * std::sin( theta_ref );
    double y = pos(1) + ey * std::cos( theta_ref );

    std::vector<double> xy;
    xy.push_back(x);
    xy.push_back(y);
    xy.push_back(theta_ref);
    return xy;
  }

  void Plotting::plot_line(Track track, std::vector<State> states, std::string col) {
    ArcLengthSpline sp = track.sp;

    std::vector<double> plot_x;
    std::vector<double> plot_y;

    for (int i=0; i<states.size(); i++) {
      State state = states[i];
      std::vector<double> xy = get_x_y(sp, state.s, state.ey);

      double x = xy[0];
      double y = xy[1];
      double theta_ref = xy[2];
      plot_x.push_back(x);
      plot_y.push_back(y);
    }

    plt::plot(plot_x,plot_y, col);
  }

  void Plotting::plotVehicle_inferred_in_xy_coordinate(Track track, const State &state, std::string boxcolor, std::string centercolor, int f) const {
    ArcLengthSpline sp = track.sp;

    std::vector<double> plot_x;
    std::vector<double> plot_y;

    std::vector<double> xy = get_x_y(sp, state.s, state.ey);

    double x = xy[0];
    double y = xy[1];
    double theta_ref = xy[2];
    plot_x.push_back(x);
    plot_y.push_back(y);

    std::vector<double> corner_x;
    std::vector<double> corner_y;

    const double lf = 0.033;
    const double lr = 0.029;
    const double dside = 0.015;

    double body_xlf =  std::cos(theta_ref + state.epsi)*lf;
    double body_xlr =  std::cos(theta_ref + state.epsi)*lr;
    double body_xw =  std::sin(theta_ref + state.epsi)*dside;
    double body_ylf =  std::sin(theta_ref + state.epsi)*lf;
    double body_ylr =  std::sin(theta_ref + state.epsi)*lr;
    double body_yw = -std::cos(theta_ref + state.epsi)*dside;

    corner_x.push_back(x + body_xlf + body_xw);
    corner_x.push_back(x + body_xlf - body_xw);
    corner_x.push_back(x - body_xlr - body_xw);
    corner_x.push_back(x - body_xlr + body_xw);
    corner_x.push_back(x + body_xlf + body_xw);

    corner_y.push_back(y + body_ylf + body_yw);
    corner_y.push_back(y + body_ylf - body_yw);
    corner_y.push_back(y - body_ylr - body_yw);
    corner_y.push_back(y - body_ylr + body_yw);
    corner_y.push_back(y + body_ylf + body_yw);

    plt::plot(corner_x,corner_y, boxcolor);
    plt::plot(plot_x,plot_y, centercolor);

    //if (f)
    //plt::fill(corner_x, corner_y, "b");
  }

  Eigen::Vector2d s_ey_2_x_y(ArcLengthSpline sp, double s, double ey) {
    Eigen::Vector2d dpos_ref = sp.getDerivative(s);
    double dx_ref = dpos_ref(0);
    double dy_ref = dpos_ref(1);
    double theta_ref = std::atan2(dy_ref,dx_ref);

    Eigen::Vector2d xy;
    Eigen::Vector2d pos = sp.getPostion(s);
    xy(0) = pos(0) - ey * std::sin( theta_ref );
    xy(1) = pos(1) + ey * std::cos( theta_ref );

    return xy;
  }

  void Plotting::drawPlot() {
    plt::draw();
    plt::pause(0.001);
  }
  void Plotting::finishPlot() {
    plt::show();
  }

  void Plotting::writeTitle(std::string s) {
    plt::title(s);
  }

  void Plotting::savePlot(std::string s) {
    plt::save(s);
  }

}
