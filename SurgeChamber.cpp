
#define _USE_MATH_DEFINES
#include "CheckValve.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include "my_tools.h"
#include "PSToolboxBaseEdge.h"
#include "SurgeChamber.h"

#include <cmath>

SurgeChamber::SurgeChamber(const string _name,
//		const string _n1,
		const string _n2,
		double _rho,
		double _Vmax,
		double _Vgas_ini,
		double _pgas_ini,
		bool _save_data): PSToolboxBaseEdge("SurgeChamber",_name,"not used",_n2), Units() {

	rho=1000.;
Vmax=_Vmax;
Vgas=_Vgas_ini;
//pgas=_pgas_ini;
Vgas_old=_Vgas_ini;
pgas_old=_pgas_ini;
p1=_pgas_ini;
p2=_pgas_ini;
Q=0.;
t_out=0;
dt_out=10.;

	save_data=_save_data;

	t=0.;
	t_old=0.;
	dt=0.1;//Tmax/100.;
	fname = name + ".dat";

	is_rigid_element=true;

	double dp = 0.1e5;
  double Q_szelep = 1.; // 3600 m3/h=1 m3/s mellett legyen 0.1bar pressure loss
  double k_open = dp / Q_szelep / Q_szelep;
  zeta_open = k_open / 1000 / 9.81;
  zeta_closed = 1.0e6*zeta_open;

  is_closed=1;

}

double SurgeChamber::H_fun(double QQ){
	double HH=0.;
if (is_closed)
		HH=zeta_closed*QQ*fabs(QQ)/rho/g;
	else
		HH=zeta_open*QQ*fabs(QQ)/rho/g;

	return HH;
}

void SurgeChamber::H_fun_lin(double QQ, double& a0, double& a1){

	if (is_closed)
		a1=2.*zeta_closed*fabs(QQ)/rho/g;
	else
		a1=2.*zeta_open*fabs(QQ)/rho/g;

	a0=H_fun(QQ)-a1*QQ;

	// cout<<endl;
	// cout<<endl<<"\t CheckValve "<<name<<": H_fun_lin evaluated at Q="<<QQ*3600<<"m3/h";
	// cout<<endl<<"\t\t p1                = "<<p1/rho/g<<" mwc";
	// cout<<endl<<"\t\t p2                = "<<p2/rho/g<<" mwc";
	// cout<<endl<<"\t\t dp                = "<<(p1-p2)/rho/g<<" mwc";
	// cout<<endl<<"\t\t H_fun(QQ)         = "<<H_fun(QQ)<<" mwc";
	// cout<<endl<<"\t\t a1                = "<<a1;
	// cout<<endl<<"\t\t a1*QQ             = "<<a1*QQ;
	// cout<<endl<<"\t\t H_fun(QQ,n)-a1*QQ = "<<H_fun(QQ)-a1*QQ;
	// cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1/rho/g<<" + "<<a0<<" + "<<a1<<" * Q (mwc)";
	// cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1<<" + "<<a0*rho*g<<" + "<<a1*rho*g<<" * Q (Pa)";
	// cout<<endl;

}



void SurgeChamber::GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q){
	// Equation: p1-p2 = rho*g*H_fun(Q) = rho*g*(a0 + a1*Q)
	// p1-rho*g*a0=p2+rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,a0,a1);
	LHS=p1-rho*g*a0;
	coeff_Q=rho*g*a1;
}

void SurgeChamber::GetEdgeEquationCoeffs(double t_target, bool is_front, double & LHS, double & coeff_Q, double & coeff_p1, double & coeff_p2){
	// Equation: LHS = coeff_p1*p1-coeff_p2*p2+coeff_Q*Q
	// p1-p2=rho*g*a0+rho*g*a1*Q
	// -rho*g*a0=-p1+p2+rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,a0,a1);
	//cout<<endl<<"CheckValve "<<name<<" evaluated at Q="<<Q*3600<<"m3/h, H="<<H_fun(Q)<<", a1="<<a1;

	LHS=rho*g*a0-p1;
	coeff_Q=-rho*g*a1;
	coeff_p1=0.;
	coeff_p2=-1.;
}



void SurgeChamber::Set_BC_Right(string type, double val){
	
	bool is_ok=false;

	if (type.compare("Pressure")==0){
		p2=val;
		//printf("\n\t p2=%5.2e",p2);
		//Q=Get_Q((p2-p1)/rho/g);
		double a0,a1;
		H_fun_lin(Q,a0,a1);
		if (fabs(a1)<1.e-8)
			Q=0.;
		else
			Q=((p1-rho*g*a0)-val)/a1/rho/g;
		is_ok=true;
	}
	if (type.compare("FlowRate")==0){
		Q=val;
		p2=p1+H_fun(Q)*rho*g;
		is_ok=true;
	}
	if (!is_ok){
		cout<<endl<<endl<<"!!!!!!!!!!!!!! SurgeChamber::Set_BC_Right !!!!!!!!!!!!!!";
		cin.get();
	}


	if (DEBUG){
		cout<<endl<<"SurgeChamber::Set_BC_Right("<<type<<","<<val<<")";
 		cout<<endl<<"\t Edge "<<name;
 		cout<<endl<<"\t type:  "<<type;
 		cout<<endl<<"\t val:   "<<val;		
		cout<<endl<<"\t p2        = "<<p2<<" Pa = "<<p2/rho/g<<" mwc";
		cout<<endl<<"\t p1 = pgas = "<<p1<<" Pa = "<<p1/rho/g<<" mwc";
		cout<<endl<<"\t Q         = "<<Q<<" m3/s";
		cout<<endl<<"\t H(Q)      = "<<H_fun(Q)<<" mwc, (p2-p1)="<<(p2-p1)/rho/g<<" mwc"<<endl;
		//cin.get();
	}
	
}

double SurgeChamber::Get_Q(double dh){

	double QQ;
	if (dh<0)
		QQ=sqrt(fabs(dh)/zeta_open);
	else
		QQ=-sqrt(dh/zeta_closed);

	if (DEBUG){
		cout<<endl<<"Surge chamber "<<name<<" : entering DEBUG mode, Get_Q() :";
		printf("\n\t QQ=%5.3f, H=%5.3f, dh=%5.3f",QQ,H_fun(QQ),dh);
		cin.get();
	}

	return QQ;
}



void SurgeChamber::GetBetaAtFront(double t_target, double & LHS, double & coeff_Q){
	cout<<endl<<"ERROR!!! SurgeChamber::GetBetaAtFront() - This function should not be called..."<<endl;
	exit(-1);

}

void SurgeChamber::Set_BC_Left(string type, double val){
	cout<<endl<<"ERROR!!! SurgeChamber::Set_BC_Left - This function should not be called..."<<endl;
		exit(-1);
}


void SurgeChamber::UpdateInternal(double t_target){

	bool ok=false;
	for (int i=0; i<transient_type.size(); i++){
		if (t>transient_time.at(i)){


			if (transient_type.at(i)=="open"){
				if (!transient_triggered.at(i)){
					cout<<endl<<"!!! Triggering transient event #"<<i<<" of surge tank "<<name;
					cout<<endl<<"!!! type: "<<transient_type.at(i);
					cout<<", time:"<<transient_time.at(i);
					transient_triggered.at(i)=true;
				}
				is_closed=false;
				ok=true;
			}

			
		}
	}
UpdateTime(t_target);
}


void SurgeChamber::UpdateTime(double _t){
double tmp=Vgas_old*pgas_old;
Vgas_old=Vgas;
pgas_old=p1;

Vgas = Vgas_old + (_t-t)*Q;
p1 = tmp/Vgas;
if (t>t_out){
printf("\n %10s t=%5.3fs, Q=%+5.2f m3/s, Vgas=%5.1f m3, pgas=p1=%5.3f vom, p2=%5.3f vom, is_closed:%d",
	name.c_str(),t, Q, Vgas, p1/1000/9.81,p2/1000/9.81,is_closed);
//cout<<endl<<"\t is_closed:"<<is_closed<<", zeta_open="<<zeta_open<<", zeta_closed="<<zeta_closed<<endl;
t_out+=dt_out;
}

t_old=t;
	t=_t;
}

void SurgeChamber::Ini(int){
	//cout<<endl<<"SurgeChamber::Ini";
	if (save_data){
		tmpvec.push_back(t);
		tmpvec.push_back(p1);
		tmpvec.push_back(p2);
		tmpvec.push_back(Q);
		tmpvec.push_back(Vgas);
		tmpvec.push_back(p2);
		data.clear();
		data.reserve(100);
		data.push_back(tmpvec);
	}
	//cout<<" done"<<endl;
};
void SurgeChamber::Ini(){
	Ini(1);
};

void SurgeChamber::Ini(double){
	Ini(1);
};


string SurgeChamber::Info(){

	std::ostringstream oss;
	oss << "\n\n SurgeChamber name     : " << name;
	oss << "\n================================";
	oss << "\n        node_from : " << node_from;
	//oss << "\n          node_to : " << node_to;
	oss << "\n              rho : " << rho<<" kg/m3";
	oss << "\n                t : " << t<<" s";
	oss << "\n      Q (current) : " << Q<<" m3/s = "<<Q*3600<<" m3/h";
	oss << "\n   Vgas (current) : " << Vgas<<" m3";
	oss << "\n   pgas (current) : " << p1<<" bar (=p1)";
	//oss << "\n     p1 (current) : " << p1<<" Pa = "<<p1/rho/g<<" mwc";
	//oss << "\n     p2 (current) : " << p2<<" Pa = "<<p2/rho/g<<" mwc"<<endl;
	oss << "\n           Tmax : "<<Tmax<<" s (estimated time to stop from n_nom, Q_nom)";
	oss << "\n             dt : "<<dt<<" s";
	oss << endl;
	return oss.str();
}

void SurgeChamber::GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void SurgeChamber::GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void SurgeChamber::Save_data(){

	tmpvec.at(0) = t;
	tmpvec.at(1) = p1;
	tmpvec.at(2) = p2;
	tmpvec.at(3) = Q;
	tmpvec.at(4) = p2;
	tmpvec.at(5) = Vgas;

	data.push_back(tmpvec);

}



void SurgeChamber::Write_data(string folder){
	
	if (!save_data) {
		cout << endl << "WARNING! SCP: " << name;
		cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
	}
	else {
		cout << endl << "Saving to " << fname.c_str() << " ... ";
		FILE * pFile;
		pFile = fopen ((folder + "/" + fname).c_str(), "w");
		fprintf(pFile, "t (s); p1 (bar); p2 (bar); Q m3/h; pgas (bar); Vgas (m3)\n");
		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e;\n",
					data.at(i).at(0),
					data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
					data.at(i).at(3)*3600,
					data.at(i).at(4)/1.e5, data.at(i).at(5)
					);
		fclose (pFile);
		cout << " done. ";
	}
}

void SurgeChamber::Add_transient(string _type, double _when, double _val){
	transient_type.push_back(_type);
	transient_time.push_back(_when);
	transient_val.push_back(_val);
	transient_triggered.push_back(false);
}

double SurgeChamber::Get_dprop(string){return 0.;};
