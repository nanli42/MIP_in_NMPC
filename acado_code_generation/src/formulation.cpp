#include <acado_code_generation.hpp>

#define GENERATE_LV_TRAJ 0
#define singRC 0
#define METHOD 1

int main( )
{
	USING_NAMESPACE_ACADO

	DifferentialState         ey, ephi, vx, vy, w, t, d, delta, s;
	IntermediateState         dsdt, Frx, Fry, Ffy;
	IntermediateState        	alphaF, alphaR;

	Control                   dd, ddelta;
	DifferentialEquation      f;


	double m = 0.041;double Iz = 27.8e-6;double lf = 0.029;double lr = 0.033;
	double Cm1=0.287;double Cm2=0.0545;double Cr0=0.0518;double Cr2=0.00035;
	double Br = 3.3852;double Cr = 1.2691;double Dr = 0.1737;
	double Bf = 2.579;double Cf = 1.2;double Df = 0.192;


#if ((GENERATE_LV_TRAJ == 1) || (singRC == 1))
	OnlineData ks;
#else
	OnlineData ks,
		t_est,
		s_LV_k, s_LV_b,
		ey_LV_k, ey_LV_b,
		ks_LV;
#endif

  const double s_start =  0.0;
  const double s_end = 1.80;

  dsdt = (vx*cos(ephi)-vy*sin(ephi)) / (1-ey*ks);
  alphaF = -atan((w*lf+vy)/vx) + delta;
  alphaR = atan((w*lr-vy)/vx);

  Frx = (Cm1-Cm2*vx)*d - Cr0 - Cr2*vx*vx;
  Fry = Dr*sin(Cr*atan(Br*alphaR));
  Ffy = Df*sin(Cf*atan(Bf*alphaF));

  f << dot(ey) == 1/dsdt * (vx*sin(ephi) + vy*cos(ephi));
  f << dot(ephi) == 1/dsdt * w - ks;
  f << dot(vx) == 1/dsdt * (w*vy + 1/m*(Frx-Ffy*sin(delta)));
  f << dot(vy) == 1/dsdt * (-w*vx + 1/m*(Fry+Ffy*cos(delta)));
  f << dot(w) == 1/dsdt * (1/Iz * (lf*Ffy*cos(delta)-lr*Fry));
  f << dot(t) == 1/dsdt;
  f << dot(d) == 1/dsdt * dd;
  f << dot(delta) == 1/dsdt * ddelta;
  f << dot(s) == 1;

  OCP ocp( s_start, s_end, 30 );
	ocp.minimizeMayerTerm( t );

  ocp.setModel( f );

	double ey_bound = (0.17);
  ocp.subjectTo( -ey_bound <= ey <= ey_bound );
  ocp.subjectTo( -1.5 <= ephi <= 1.5 );

#if GENERATE_LV_TRAJ == 1
	ocp.subjectTo( 0.05 <= vx <= 1.2 );
#elif singRC == 1
	ocp.subjectTo( 0.05 <= vx <= 1.6 );
#else
	ocp.subjectTo( 0.05 <= vx <= 1.6 );
#endif

	ocp.subjectTo( -1.0 <= vy <= 1.0 );

  ocp.subjectTo( -8.0 <= w <= 8.0 );

  ocp.subjectTo( 0.0 <= t <=100.0);

  ocp.subjectTo( -1.0 <= d <= 1.0  );
  ocp.subjectTo( -0.6 <= delta <= 0.6 );

  ocp.subjectTo( -10.0 <= dd <= 10.0 );
  ocp.subjectTo( -10.0 <= ddelta <= 10.0);

#if GENERATE_LV_TRAJ == 1
	ocp.subjectTo( vx*vx+vy*vy <= 1.2*1.2 );
#elif singRC == 1
	ocp.subjectTo( vx*vx+vy*vy <= 1.6*1.6 );
#else
	ocp.subjectTo( vx*vx+vy*vy <= 1.6*1.6 );
#endif

	const double e_long =  0.9;
  const double e_eps  =  0.95;
	ocp.subjectTo( (e_long*Frx)*(e_long*Frx) + Fry*Fry - (e_eps*Dr)*(e_eps*Dr) <= 0.0 ); //1e10); //
	const double max_alphaF =  0.6;
	ocp.subjectTo( -max_alphaF <= alphaF <= +max_alphaF );

	double delta_s_0 = lf+lr;
	double delta_ey_0 = 2*0.015;
	IntermediateState delta_s = (delta_s_0 * cos(ephi) * (1/(1-ey*ks)))/2;
	IntermediateState delta_ey_pos, delta_ey_neg;
	delta_ey_pos = (delta_ey_0 * cos(ephi) + delta_s_0 * sin(ephi))/2;
	delta_ey_neg = (delta_ey_0 * cos(ephi) + delta_s_0 * -sin(ephi))/2;

	double tol = 1e-3;
	double xx = 0.035;

#if ((GENERATE_LV_TRAJ == 1) || (singRC == 1))
	ocp.subjectTo( ey + xx <= 0.17);
	ocp.subjectTo( -0.17 <= ey - xx );

	ocp.setNOD(1);
#else

	IntermediateState s_LV = s_LV_k*(t-t_est)+s_LV_b;
	IntermediateState ey_LV = ey_LV_k*(t-t_est)+ey_LV_b;

	#if (METHOD == 0)
		ocp.subjectTo( ey + xx <= 0.17);
		ocp.subjectTo( -0.17 <= ey - xx );

		ocp.subjectTo( s + 1/ks*asin(xx/(1/ks-ey)) - (s_LV - 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) +0.01 <= 1e+6 );

		ocp.subjectTo( ey + xx*2 - ey_LV +0.01 <= 1e+6);
		ocp.subjectTo( -1e+6 <= ey - xx*2 - ey_LV -0.01);
	#else
		ocp.subjectTo( ey + xx <= 0.17);
		ocp.subjectTo( -0.17 <= ey - xx );

		ocp.subjectTo( -0.9 <= xx/(1/ks-ey) <= 0.9);
		ocp.subjectTo( -0.9 <= xx/(1/ks_LV-ey_LV) <= 0.9);

		ocp.subjectTo( s + 1/ks*asin(xx/(1/ks-ey)) - (s_LV - 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) <= 1e+6 );

		ocp.subjectTo( -1e+6 <= s - 1/ks*asin(xx/(1/ks-ey)) - (s_LV + 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) );

		ocp.subjectTo( ey + xx*2 - ey_LV +tol <= 1e+6);
		ocp.subjectTo( s - 1/ks*asin(xx/(1/ks-ey)) - (s_LV + 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) <= 1e+6 );
		ocp.subjectTo( -1e+6 <= s + 1/ks*asin(xx/(1/ks-ey)) - (s_LV - 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) );

		ocp.subjectTo( -1e+6 <= ey - xx*2 - ey_LV -tol );
		ocp.subjectTo( s - 1/ks*asin(xx/(1/ks-ey)) - (s_LV + 1/ks_LV*asin(xx/(1/ks_LV-ey_LV)))  <= 1e+6 );
		ocp.subjectTo( -1e+6 <= s + 1/ks*asin(xx/(1/ks-ey)) - (s_LV - 1/ks_LV*asin(xx/(1/ks_LV-ey_LV))) );
	#endif

	ocp.setNOD(7);
#endif

	OCPexport mpc( ocp );

	mpc.set( HESSIAN_APPROXIMATION,       EXACT_HESSIAN    );
	mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING );
	mpc.set( INTEGRATOR_TYPE,             INT_RK4         );
	mpc.set( NUM_INTEGRATOR_STEPS,        100              );
	mpc.set( QP_SOLVER,                   QP_QPOASES      );
	mpc.set( HOTSTART_QP,                 NO             		);
	mpc.set( GENERATE_TEST_FILE,          NO             );
	mpc.set( GENERATE_MAKE_FILE,          NO             );
	mpc.set( GENERATE_MATLAB_INTERFACE,   NO             );
	mpc.set( SPARSE_QP_SOLUTION, 		  FULL_CONDENSING_N2	);
	mpc.set( DYNAMIC_SENSITIVITY, 		  	FORWARD_OVER_BACKWARD				);
	mpc.set( GENERATE_SIMULINK_INTERFACE, NO             );
	mpc.set( CG_HARDCODE_CONSTRAINT_VALUES, NO 					);
	mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES 				);

	mpc.set( LEVENBERG_MARQUARDT, 1e-6 );

	if (mpc.exportCode( "../../external/acado_generated" ) != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

	mpc.printDimensionsQP( );

	return EXIT_SUCCESS;
}
