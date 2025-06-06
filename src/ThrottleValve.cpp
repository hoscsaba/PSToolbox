#define _USE_MATH_DEFINES

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

#include "ThrottleValve.h"
#include "my_tools.h"

ThrottleValve::ThrottleValve(const string _name,
    const double _A,
    const double _Cd,
    const double _ro,
    const bool _save_data) : Units() {
  save_data = _save_data;
  name = _name;
   A_flowthrough  = _A;
  Cd = _Cd;
  ro = _ro;

  is_Gas = false;
  ini_done = false;
  fname = name + ".dat";

  mp_nevl = Get_MassFlowRate(2.e5, 293., 1.e5, 293.);
}

ThrottleValve::ThrottleValve(const string _name,
    const double _A,
    const double _Cd,
    Gas* _gas,
    const bool _save_data) : Units() {
  save_data = _save_data;
  name = _name;
   A_flowthrough  = _A;
  Cd = _Cd;
  gas = _gas;

  is_Gas = true;
  ini_done = false;
  fname = name + ".dat";

  mp_nevl = Get_MassFlowRate(2.e5, 293., 1.e5, 293.);

}

ThrottleValve::~ThrottleValve() {}

string ThrottleValve::GetName() {
  return name;
}

double ThrottleValve::Get_MassFlowRate(double p_upstream, double T_upstream, double p_downstream, double T_downstream){
  double mp_;
  if (is_Gas)
    mp_ = Get_MassFlowRate_Compressible(p_upstream, T_upstream, p_downstream, T_downstream);
  else
    mp_ = Get_MassFlowRate_InCompressible(p_upstream, p_downstream, ro);
  mp=mp_; // mp is used only for storing lates mass flow rate for Get_dprop("mp");
  return mp_;
}

double ThrottleValve::Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream, double T_downstream) {

  double massflux = gas->Get_MassFlux(p_upstream, T_upstream, p_downstream, T_downstream);
  return Cd*A_flowthrough*massflux;
}

double ThrottleValve::Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream) {
  double rho_upstream = gas->Get_rho(p_upstream, T_upstream);
  double kappa = gas->Get_kappa_Tp();

  double exponent = (kappa + 1.) / (kappa - 1.);
  double tmp = kappa * pow(2. / (kappa + 1.), exponent);
  return Cd * A_flowthrough * sqrt(rho_upstream * p_upstream * tmp);
}

double ThrottleValve::Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream) {
  double rho_upstream = gas->Get_rho(p_upstream, T_upstream);
  double kappa = gas->Get_kappa_Tp();

  double exp1 = 1. / kappa;
  double exp2 = (kappa + 1.) / kappa;
  double tmp = 2.*rho_upstream * p_upstream * kappa / (kappa - 1);
  double pr = p_downstream / p_upstream;
  return Cd * A_flowthrough * sqrt(tmp * (pow(pr, exp1) - pow(pr, exp2)));
}

double ThrottleValve::Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho) {
  double dp = p_upstream - p_downstream;
  double DP_MIN = 0.01e5;
  double G0;
  if (fabs(dp) < DP_MIN) {
    G0=sqrt(2. * ro * DP_MIN)*dp / DP_MIN;
    mp = Cd * A_flowthrough * G0;
  }
  else{
    G0=sqrt(2. * ro * fabs(dp)) * signum(dp);
    mp = Cd * A_flowthrough * G0;
  }

  return mp;
}

double ThrottleValve::signum(double x) {
  // if (x > 0.) return 1.;
  // if (x < 0.) return -1.;
  // return 0.;
  if (x > 0.)
    return 1.;
  else
    return -1.;
}
