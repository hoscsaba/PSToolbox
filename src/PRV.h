#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "PSToolboxBaseEdge.h"
#include "Units.h"

using namespace std;

class PRV: public PSToolboxBaseEdge, Units
{

private:
	double H_fun(double Q);
	void H_fun_lin(double Q, double& a0, double& a1);
	double rho;
	double zeta;
	double pset;
	double Get_Q(double dh);
	vector< vector<double> > data;
	vector<double> tmpvec;
	string fname;
	double Tmax;

public:
	PRV(const string _name,
			const string _cspe_name,
			const string _cspv_name,
			double _rho,
			double _pset, bool save_data);
	~PRV(){};
	void Ini(int);
	void Ini();
	void Ini(double vini, double pstart){};

	string Info();
	void calculateZeta(double);

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

