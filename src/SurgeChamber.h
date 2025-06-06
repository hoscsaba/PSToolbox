#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "PSToolboxBaseEdge.h"
#include "Units.h"

using namespace std;

class SurgeChamber: public PSToolboxBaseEdge, Units
{

public:
	double rho, Vgas, Vgas_old, pgas_old, Vmax, t_old;
	double zeta_open, zeta_closed;
	double Tmax;
	vector< vector<double> > data;
	vector<double> tmpvec;
	string fname;
	int state;
	double last_statechange;
	bool is_closed;
	double Get_Q(double dh);
	double H_fun(double Q);
	void H_fun_lin(double Q, double& a0, double& a1);
	double t_out, dt_out;

public:
	SurgeChamber(const string _name,
//			const string _cspe_name,
			const string _cspv_name,
			double _rho,
			double Vmax, double Vgas_ini, double pgas_ini, bool save_data);
	~SurgeChamber(){};
	void Ini(int);
	void Ini();
	void Ini(double);
	void Ini(double vini, double pstart){};

	string Info();

	void GetBetaAtFront(double, double &, double &);
	void GetAlphaAtEnd(double, double &, double &);
	void GetEdgeEquationCoeffs(double, bool, double &, double &, double &, double &);

	void Set_BC_Left(string type, double val);
	void Set_BC_Right(string type, double val);

	void GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&);
	void GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&);
	void Save_data();
	void Write_data(string);
	double Get_dprop(string);
	void UpdateInternal(double);
	void UpdateTime(double);
	void Set_string_prop(string,string){};
	void Add_transient(string,double,double);

};
