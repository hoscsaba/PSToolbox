#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Units.h"
#include "IdealGas.h"

using namespace std;
using namespace Eigen;

/*! \brief A class for handling both liquid and fluid valves
  A class for storing data regarding a valve and allowing for a solution of the differential equation
  describing the movement of the spring, mass damper system.
  */
class ThrottleValve : public Units {
private:
    bool save_data; //!< whether to save the data of the valve or not
    string name; //!< name of the valve
    double A_flowthrough;
    double Cd; /*!<discharge coefficient of the valve*/
    double ro /*Desity of fluid in the valve */, Dbore /*Smallest cross section in valve*/, xe
            /*!< spring pre-compression */, p_set /* set pressure */;
    double mp_nevl; /*Design pressure */
    double mp; // actual mass flow rate
    double t/*!< current time of the simulation*/, x /*!< displacement of the valve during the simulation */, v
            /*!< velocity of the valve head during the simulation*/, p /*inlet side pressure in the valve */;
    bool ini_done; //!< whether the system ahs been initialized

    vector<vector<double> > data;
    vector<double> tmpvec;
    Gas *gas;
    bool is_Gas;

    double Get_MassFlowRate_Compressible(double p_upstream, double T_upstream, double p_downstream,
                                         double T_downstream);

    double Get_MassFlowRate_Compressible_Choked(double p_upstream, double T_upstream);

    double Get_MassFlowRate_Compressible_UnChoked(double p_upstream, double T_upstream, double p_downstream);

    double Get_MassFlowRate_InCompressible(double p_upstream, double p_downstream, double rho);

    double signum(double x);

public:
    ThrottleValve(const string _nev,
                  const double _A,
                  const double _Cd,
                  const double _ro,
                  const bool save_data);

    ThrottleValve(const string _nev,
                  const double _A,
                  const double _Cd,
                  Gas *_gas,
                  const bool save_data);

    ~ThrottleValve();

    string GetName();

    string fname;

    double Get_MassFlowRate(double p_upstream, double T_upstream, double p_downstream, double T_downstream);
};
