#define _USE_MATH_DEFINES
#include "Tank.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include "my_tools.h"
#include "PSToolboxBaseEdge.h"
#include "Tank.h"

#include <cmath>

Tank::Tank(const string _name,
		const string _n1,
		const string _n2,
		double _rho,
		double _p_total,
		bool _save_data): PSToolboxBaseEdge("Tank",_name,"not used",_n2), Units() {

	rho = _rho;
	p_total=_p_total;
	save_data=_save_data;

	Tmax = 1.;
	t=0.;
	Q=28./3600;
	dt=0.01;
	fname = "data/" + name + ".dat";

	is_rigid_element=true;
	p1=5.e5;
	p2=5.e5;
}

double Tank::H_fun(double QQ){
	return 0.;
}

void Tank::H_fun_lin(double QQ, double& a0, double& a1){
	a0=0.;
	a1=0.;
	cout<<endl<<"ERROR!!! Tank::GetBetaAtFront() - This function should not be called..."<<endl;
	exit(-1);

}

void Tank::GetBetaAtFront(double t_target, double & LHS, double & coeff_Q){
	cout<<endl<<"ERROR!!! Tank::GetBetaAtFront() - This function should not be called..."<<endl;
	exit(-1);

}

void Tank::GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q){
	cout<<endl<<"ERROR!!! Tank::GetBetaAtFront() - This function should not be called..."<<endl;
	exit(-1);

}

void Tank::GetEdgeEquationCoeffs(double t_target, bool is_front, double & LHS, double & coeff_Q, double & coeff_p1, double & coeff_p2){

	LHS=p_total;
	coeff_Q=0.;
	coeff_p1=0.;
	coeff_p2=1.;
}

void Tank::Set_BC_Left(string type, double val){
	cout<<endl<<"ERROR!!! Tank::Set_BC_Left - This function should not be called..."<<endl;
		exit(-1);
}

void Tank::Set_BC_Right(string type, double val){
	cout<<endl<<"ERROR!!! Tank::Set_BC_Right - This function should not be called..."<<endl;
		exit(-1);
}

double Tank::Get_Q(double dh){
	cout<<endl<<"ERROR!!! Tank::Get_Q(double dh) - This function should not be called..."<<endl;
		exit(-1);

	return 0.;
}


void Tank::UpdateInternal(double t_target){
}


void Tank::UpdateTime(double _t){

	t=_t;
}

void Tank::Ini(int){

};
void Tank::Ini(){
	Ini(1);
};

string Tank::Info(){

	std::ostringstream oss;
	oss << "\n\n Tank name     : " << name;
	oss << "\n================================";
	oss << "\n        node_from : " << node_from;
	oss << "\n          node_to : " << node_to;
	oss << "\n              rho : " << rho<<" kg/m3";
	oss << "\n          p_total : " << p_total<<" Pa";
	oss << "\n                t : " << rho<<" s";
	oss << "\n                Q : " << t<<" m3/s = "<<Q*3600<<" m3/h (current)";
	oss << "\n             dt : "<<dt<<" s";
	oss << endl;
	return oss.str();
}

void Tank::GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void Tank::GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void Tank::Save_data(){}

void Tank::Write_data(string){}

void Tank::Add_transient(string _type, double _when, double _val){}

double Tank::Get_dprop(string){return 0.;};
