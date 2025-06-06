#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include "Units.h"
#include "IdealGas.h"
//#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"

using namespace std;
using namespace Eigen;

//! Class for modeling Damper
/*!
  The class contains data for Damper
  */
class Damper : public Units {
private:
    bool save_data;
    string name;
    double vol, a, p, t, rho, T, mp;
    double K, dia;
    int num_of_orifices;
    double beta;

    //  bool is_Gas;
    VectorXd DamperODE(const double x, const VectorXd &y, const VectorXd &pars);

    void rk45(double x_0,
              double x_max,
              VectorXd y_0,
              vector<double> &x, vector<vector<double> > &y,
              const VectorXd &pars);

    vector<vector<double> > data;
    vector<double> tmpvec;
    bool Ini_done;
    double n_poly;
    bool DEBUG;

public:
    bool is_Gas;

    Damper(string _nev,
              double _vol,
              double _a,
              double _rho,
              int _num_of_orifices,
                 double diameter,
              bool save_data);

    Damper(string _name,
              double _vol,
              IdealGas *gas,
              double n_poly,
              int _num_of_orifices,
                 double diameter,
              bool _save_data);

    Damper(string _name,
              double _vol,
              IdealGas *gas,
              double n_poly,
              int _num_of_orifices,
                 double diameter,
              bool _save_data,
              double Tt);

    Damper(string _name,
              double _vol,
              Gas *gas,
              double n_poly,
              int _num_of_orifices,
                 double diameter,
              bool _save_data);

    ~Damper();

    void SetDebug(bool DEBUG);

    string GetName();

    string Info();

    void Update(double t_target, double mp_in);

    void UpdateDimlessPars(double pref, double mp_nevl, double omega);

    //void Ini();
    void Ini(double _pstart);

    void Ini(double _pstart, double _Tstart);

    double Get_dprop(string prop_string);

    void Set_dprop(string prop_string, double val);

    double GetP() { return p; }

    void List_data();

    void Save_data();

    vector<double> Get_dvprop(string prop_string);

    //void Plot();
    string fname;
    Gas *gas;

    void Plot(bool show, bool save);

    double signedSqrt(double x);

    double GetMassFlowRate(double p_connect);

    double GetInternalPressureDrop();
};
