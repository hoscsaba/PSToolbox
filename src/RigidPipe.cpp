#define _USE_MATH_DEFINES
#include "RigidPipe.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include "my_tools.h"
#include "PSToolboxBaseEdge.h"


#include <cmath>

RigidPipe::RigidPipe(const string _name,
		const string _n1,
		const string _n2,
		double _rho,
		double _L,
		double _D,
		double _lambda,
		bool _save_data): PSToolboxBaseEdge("RigidPipe",_name,_n1,_n2), Units() {

	rho = _rho;
	L=_L;
	D=_D;
	lambda=_lambda;
	A=D*D*M_PI/4;
	zeta_open = lambda*L/D/2/g/A/A;
	zeta_closed = zeta_open;
	save_data=_save_data;

	Tmax = 1.;
	t=0.;
	Q= 1.0 / 3600;
	dt=0.1;
	fname = name + ".dat";

	is_rigid_element=true;
	p1=3.e5;
	p2=3.e5;
}

double RigidPipe::H_fun(double QQ){
	double HH=0.;
	if (QQ<0)
		HH=zeta_closed*QQ*fabs(QQ);
	else
		HH=zeta_open*QQ*fabs(QQ);

	return HH;
}

void RigidPipe::H_fun_lin(double QQ, double& a0, double& a1){

	if (QQ<0)
		a1=2.*zeta_closed*fabs(QQ);
	else
		a1=2.*zeta_open*fabs(QQ);

	a0=H_fun(QQ)-a1*QQ;
if (isnan(QQ)){
		cout<<endl;
		cout<<endl<<"\t RigidPipe "<<name<<": H_fun_lin evaluated at Q="<<QQ*3600<<" m3/h";
		cout<<endl<<"\t\t p1                = "<<p1/rho/g<<" mwc";
		cout<<endl<<"\t\t p2                = "<<p2/rho/g<<" mwc";
		cout<<endl<<"\t\t dp                = "<<(p1-p2)/rho/g<<" mwc";
		cout<<endl<<"\t\t H_fun(QQ)         = "<<H_fun(QQ)<<" mwc";
		cout<<endl<<"\t\t a1                = "<<a1;
		cout<<endl<<"\t\t a1*QQ             = "<<a1*QQ;
		cout<<endl<<"\t\t H_fun(QQ,n)-a1*QQ = "<<H_fun(QQ)-a1*QQ;
		cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1/rho/g<<" + "<<a0<<" + "<<a1<<" * Q (mwc)";
		cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1<<" + "<<a0*rho*g<<" + "<<a1*rho*g<<" * Q (Pa)";
		cout<<endl;

cin.get();
}
}

void RigidPipe::GetBetaAtFront(double t_target, double & LHS, double & coeff_Q){
	// Equation: p1-p2 = rho*g*H_fun(Q) = rho*g*(a0 + a1*Q)
	// p2+rho*g*a0=p1-rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,a0,a1);
	LHS=p2+rho*g*a0;
	coeff_Q=-rho*g*a1;

}

void RigidPipe::GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q){
	// Equation: p1-p2 = rho*g*H_fun(Q) = rho*g*(a0 + a1*Q)
	// p1-rho*g*a0=p2+rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,a0,a1);
	LHS=p1-rho*g*a0;
	coeff_Q=rho*g*a1;
}

void RigidPipe::GetEdgeEquationCoeffs(double t_target, bool is_front, double & LHS, double & coeff_Q, double & coeff_p1, double & coeff_p2){
	// Equation: LHS = coeff_p1*p1-coeff_p2*p2+coeff_Q*Q
	// p1-p2=rho*g*a0+rho*g*a1*Q
	// -rho*g*a0=-p1+p2+rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,a0,a1);
	//cout<<endl<<"RigidPipe "<<name<<" evaluated at Q="<<Q*3600<<"m3/h, H="<<H_fun(Q)<<", a1="<<a1;

	LHS=rho*g*a0;
	coeff_Q=-rho*g*a1;
	coeff_p1=1.;
	coeff_p2=-1.;
}

void RigidPipe::Set_BC_Left(string type, double val){
	p1=val;

	bool is_ok=false;
	if (type.compare("Pressure")==0){
		p1=val;
		double a0,a1;
		H_fun_lin(Q,a0,a1);
		Q=-((p2+rho*g*a0)-val)/a1/rho/g;
		//Q=Get_Q((p2-p1)/rho/g);
		is_ok=true;
	}
	if (type.compare("FlowRate")==0){
		Q=val;
		//		double a0,a1;
		//
		//	H_fun_lin(Q,a0,a1);
		p2=p1+H_fun(Q)*rho*g;
		is_ok=true;
	}
	if (!is_ok){
		cout<<endl<<endl<<"!!!!!!!!!!!!!! RigidPipe::Set_BC_Left !!!!!!!!!!!!!!";
		cin.get();
	}

	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"RigidPipe::Set_BC_Left("<<type<<","<<val<<")";
		cout<<endl<<"\t Edge "<<name<<" BC set up:";
		cout<<endl<<"\t p1   = "<<p1<<" Pa (p2="<<p2<<" Pa)";
		cout<<endl<<"\t Q    = "<<Q<<" m3/s";
		cout<<endl<<"\t H(Q) = "<<H_fun(Q)<<" mwc, (p2-p1)="<<(p2-p1)/rho/g<<" mwc"<<endl;
	}
}

void RigidPipe::Set_BC_Right(string type, double val){
	p2=val;

	bool is_ok=false;

	if (type.compare("Pressure")==0){
		p2=val;
		//Q=Get_Q((p2-p1)/rho/g);
		double a0,a1;
		H_fun_lin(Q,a0,a1);
		Q=((p1-rho*g*a0)-val)/a1/rho/g;
		is_ok=true;
	}
	if (type.compare("FlowRate")==0){
		Q=val;
		p2=p1+H_fun(Q)*rho*g;
		is_ok=true;
	}
	if (!is_ok){
		cout<<endl<<endl<<"!!!!!!!!!!!!!! RigidPipe::Set_BC_Right !!!!!!!!!!!!!!";
		cin.get();
	}

	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"RigidPipe::Set_BC_Right("<<type<<","<<val<<")";
		cout<<endl<<"\t Edge "<<name;
		cout<<endl<<"\t type:  "<<type;
		cout<<endl<<"\t val:   "<<val;
		cout<<endl<<"\t p2   = "<<p2<<" Pa = "<<p2/rho/g<<" (p1="<<p1<<" Pa)";
		cout<<endl<<"\t Q    = "<<Q<<" m3/s";
		cout<<endl<<"\t H(Q) = "<<H_fun(Q)<<" mwc, (p2-p1)="<<(p2-p1)/rho/g<<" mwc"<<endl;
	}

}

double RigidPipe::Get_Q(double dh){

	double QQ;
	if (dh<0)
		QQ=sqrt(fabs(dh)/zeta_open);
	else
		QQ=-sqrt(dh/zeta_closed);

	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"RigidPipe "<<name<<" : entering DEBUG mode, Get_Q() :";
		printf("\n\t QQ=%5.3f, H=%5.3f, dh=%5.3f",QQ,H_fun(QQ),dh);
	}

	if (DEBUG)
		cin.get();

	return QQ;
}


void RigidPipe::UpdateInternal(double t_target){
	t=t_target;
}


void RigidPipe::UpdateTime(double _t){

	t=_t;
}

void RigidPipe::Ini(int){
	//cout<<endl<<"Pump::Ini";
	if (save_data){
		tmpvec.push_back(t);
		tmpvec.push_back(p1);
		tmpvec.push_back(p2);
		tmpvec.push_back(Q);
		tmpvec.push_back(H_fun(Q));
		data.clear();
		data.reserve(100);
		data.push_back(tmpvec);
	}
	//cout<<" done"<<endl;
};

void RigidPipe::Ini(){
	Ini(1);
};


void RigidPipe::Ini(double){
	Ini(1);
};

void RigidPipe::Ini(double vini, double pstart, double dt_target){
	Ini(1);
}


string RigidPipe::Info(){

	std::ostringstream oss;
	oss << "\n\n RigidPipe name     : " << name;
	oss << "\n================================";
	oss << "\n        node_from : " << node_from;
	oss << "\n          node_to : " << node_to;
	oss << "\n              rho : " << rho<<" kg/m3";
	oss << "\n                t : " << rho<<" s";
	oss << "\n                Q : " << t<<" m3/s = "<<Q*3600<<" m3/h (current)";
	oss << "\n      H (current) : " << H_fun(Q)<<" mwc";
	oss << "\n     p1 (current) : " << p1<<" Pa = "<<p1/rho/g<<" mwc";
	oss << "\n     p2 (current) : " << p2<<" Pa = "<<p2/rho/g<<" mwc"<<endl;
	oss << "\n        zeta_open : " << zeta_open<<" m/(m3/s)^2";
	oss << "\n      zeta_closed : " << zeta_closed<<" m/(m3/s)^2";
	oss << "\n             dt : "<<dt<<" s";
	oss << endl;
	return oss.str();
}

void RigidPipe::GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void RigidPipe::GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void RigidPipe::Save_data(){

	tmpvec.at(0) = t;
	tmpvec.at(1) = p1;
	tmpvec.at(2) = p2;
	tmpvec.at(3) = Q;
	tmpvec.at(4) = H_fun(Q);

	data.push_back(tmpvec);

}

void RigidPipe::Write_data(string folder){

	if (!save_data) {
		cout << endl << "WARNING! RigidPipe: " << name;
		cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
	}
	else {
		cout << endl << "Saving to " << fname.c_str() << " ... ";
		FILE * pFile;
		pFile = fopen ((folder + "/" + fname).c_str(), "w");
		fprintf(pFile, "t (s); p1 (bar); p2 (bar); Q m3/h; dH (mwc);\n");
		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e\n",
					data.at(i).at(0),
					data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
					data.at(i).at(3)*3600, data.at(i).at(4));
		fclose (pFile);
		cout << " done. ";
	}
}

void RigidPipe::Add_transient(string _type, double _when, double _val){
	//	transient_type.push_back(_type);
	//	transient_time.push_back(_when);
	//	transient_val.push_back(_val);
	//	transient_triggered.push_back(false);
}

double RigidPipe::Get_dprop(string prop_string){
	double out=-1;

	if (prop_string == "dt")
		out = dt;
	else if (prop_string == "zeta_closed")
		out = zeta_closed;
	else if (prop_string == "zeta_open")
		out = zeta_open;
	else if (prop_string == "p_front")
		out = data.back().at(1);
	else if (prop_string == "v_front")
		out = data.back().at(3);
	else {
		cout << endl
				<< "ERROR! RigidPipe::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl;

		cout << endl << "Name of RigidPipe: " << name << endl;
		cin.get();
	}
	return out;
};
