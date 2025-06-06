#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "PSToolboxBaseEdge.h"
#include "Units.h"

using namespace std;

class Pump: public PSToolboxBaseEdge, Units
{

public:
	vector<double> coeff_H;
	vector<double> coeff_P;
	double n_nom;
	double theta;
	double Q_nom, H_nom;
	double H_fun(double Q, double n);
	double P_fun(double Q, double n);
	void H_fun_lin(double Q, double n, double& a0, double& a1);
	double rho;
	double Tmax;
	double Get_Q(double dh);
	vector< vector<double> > data;
	vector<double> tmpvec;
	string fname;
	bool use_true_H_curve; // false=linear approx., true=H curve
	int state;
	double last_statechange;

public:
	Pump(const string _name,
			const string _cspe_name,
			const string _cspv_name,
			double _rho,
			double _Q_nom, double _H_nom,
			vector<double> &coeff_H,
			vector<double> &coeff_P,
			double n_nom, double n_act,
			double theta, bool save_data);
	~Pump(){};
	void Ini(int);
	void Ini();
	void Ini(double);
	void Ini(double vini, double pstart, double dt_target);

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
