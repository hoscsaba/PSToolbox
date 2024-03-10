#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "PSToolboxBaseEdge.h"
#include "Units.h"

using namespace std;

class Pump: public PSToolboxBaseEdge, Units
{

private:
	vector<double> coeff_H;
	vector<double> coeff_P;
	double n_nom, n_act, n_act_new;
	double theta;
	double Q_nom;
	double H_fun(double Q, double n);
	double P_fun(double Q, double n);
	void H_fun_lin(double Q, double n, double& a0, double& a1);
	double p1=0.;
	double p2=0.;
	double rho;
	double Q, Tmax;
	double Get_Q(double dh);
    vector< vector<double> > data;
    vector<double> tmpvec;
string fname;

public:
	Pump(const string _name,
			const string _cspe_name,
			const string _cspv_name,
			double _rho,
			double _Q,
			vector<double> &coeff_H,
			vector<double> &coeff_P,
			double n_nom, double n_act,
			double theta, bool save_data);
	~Pump(){};
	void Ini(int);
	void Ini();

	string Info();

	void GetBetaAtFront(double, double &, double &);
	void GetAlphaAtEnd(double, double &, double &);

	void Set_BC_Left(string type, double val);
	void Set_BC_Right(string type, double val);

	void GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&);
	void GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&);
	void Save_data();
	void Write_data();
	double Get_dprop(string);
	void UpdateInternal();
	void UpdateTime(double);



};
