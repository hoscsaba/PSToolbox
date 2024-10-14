
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
#include "Pump.h"

#include <cmath>

Pump::Pump(const string _name,
		const string _n1,
		const string _n2,
		double _rho,
		double _Qnom,
		double _Hnom,
		vector<double> &_coeff_H,
		vector<double> &_coeff_P,
		double _n_nom,
		double _n_act,
		double _theta,
		bool _save_data): PSToolboxBaseEdge("Pump",_name,_n1,_n2), Units() {

	Q_nom=_Qnom;
	H_nom=_Hnom;
	rho=_rho;
	for (unsigned int i=0; i< _coeff_H.size(); i++)
		coeff_H.push_back(_coeff_H.at(i));
	for (unsigned int i=0; i< _coeff_P.size(); i++)
		coeff_P.push_back(_coeff_P.at(i));

	n_nom=_n_nom/60.;
	theta = _theta;
	n_act=_n_act/60.;
	n_act_new=n_act;
	save_data=_save_data;

	double Pmax = P_fun(Q_nom,n_nom);
	double omega0 = n_nom*2.*M_PI;
	double Mmax = Pmax/omega0;
	Tmax = omega0/Mmax/theta;
	t=0.;
	Q=Q_nom;
	dt=0.1;//Tmax/100.;
	fname = name + ".dat";

	is_rigid_element=true;
}

double Pump::H_fun(double QQ, double n){
	double HH=0.;

	if (use_true_H_curve){
		if (QQ<0)
			//HH = coeff_H.at(0)*pow(n/n_nom,2.0);
			HH=coeff_H.at(0)-(20.*coeff_H.at(0))/Q_nom*QQ;
		else
			for (unsigned int i=0; i<coeff_H.size(); i++)
				HH+=coeff_H.at(i)*pow(QQ,i)*pow(n/n_nom,(double)2-i);
	}
	else{
		double H_nom_act=H_nom*(n/n_nom)*(n/n_nom);
		double Q_nom_act=Q_nom*(n/n_nom);
		HH=1.3*H_nom_act-0.3*H_nom_act/Q_nom_act*QQ;
	}
	return HH;
}

void Pump::H_fun_lin(double QQ, double n, double& a0, double& a1)
{
	if (use_true_H_curve){
		a1=0.;
		if (QQ<0)
			a1=-20.*coeff_H.at(0)/Q_nom;
		else
			for (unsigned int i=1; i<coeff_H.size(); i++)
				a1+=coeff_H.at(i)*pow(QQ,i-1)*i*pow(n/n_nom,(double)2-i);
	}
	else{
		double H_nom_act=H_nom*(n/n_nom)*(n/n_nom);
		double Q_nom_act=Q_nom*(n/n_nom);
		a1=-0.3*H_nom_act/Q_nom_act;
	}
	a0=H_fun(QQ,n)-a1*QQ;
	//a0=H_nom-a1*(QQ-Q_nom);


	//			cout<<endl;
	//			cout<<endl<<"\t PUMP H_fun_lin evaluated at Q="<<QQ<<" m3/s = "<<QQ*3600<<" m3/h";
	//			cout<<endl<<"\t\t p1                = "<<p1<<" Pa";
	//			cout<<endl<<"\t\t p2                = "<<p2<<" Pa";
	//			cout<<endl<<"\t\t dp=(p2-p1)/rho/g  = "<<(p2-p1)/rho/g<<" mwc";
	//			cout<<endl<<"\t\t H_fun(QQ)         = "<<H_fun(QQ,n_act)<<" mwc";
	//			cout<<endl<<"\t\t a1                = "<<a1;
	//			cout<<endl<<"\t\t a1*QQ             = "<<a1*QQ;
	//			cout<<endl<<"\t\t H_fun(QQ,n)-a1*QQ = "<<H_fun(QQ,n)-a1*QQ;
	//			cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1<<" + "<<a0<<" + "<<a1<<" * Q (mwc)";
	//			cout<<endl<<"\t\t p2 = p1 + a0 + a1*QQ = "<<p1<<" + "<<a0*rho*g<<" + "<<a1*rho*g<<" * Q (Pa)";
	//			cout<<endl;
	//cin.get();

}

double Pump::P_fun(double QQ, double n){
	//if(QQ <= 0) return coeff_P[0]*pow(n/n_nom,3);

	double PP=0.;
	for (unsigned int i=0; i<coeff_P.size(); i++)
		PP+=coeff_P.at(i)*pow(QQ,i)*pow(n/n_nom,3-i);

	return  PP;
}

void Pump::GetBetaAtFront(double t_target, double & LHS, double & coeff_Q){
	// Equation: p2-p1 = rho*g*H_fun(Q) = rho*g*(a0 + a1*Q)
	// p2-rho*g*a0=p1+rho*g*a1*Q
	double a0,a1;

	H_fun_lin(Q,n_act,a0,a1);
	LHS=p2-rho*g*a0;
	coeff_Q=rho*g*a1;

}

void Pump::GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q){
	// Equation: p2-p1 = rho*g*H_fun(Q) = rho*g*(a0 + a1*Q)
	// p1+rho*g*a0=p2-rho*g*a1*Q
	double a0,a1;
	H_fun_lin(Q,n_act,a0,a1);
	LHS=p1+rho*g*a0;
	coeff_Q=-rho*g*a1;

}

void Pump::GetEdgeEquationCoeffs(double t_target, bool is_front, double & LHS, double & coeff_Q, double & coeff_p1, double & coeff_p2){
	// Equation: LHS = coeff_p1*p1-coeff_p2*p2+coeff_Q*Q
	// p2-p1=rho*g*a0+rho*g*a1*Q
	// -rho*g*a0=p1-p2+rho*g*a1*Q
	double a0,a1;
	//	double Qtmp=Q;
	//	if (Q<-Q_nom){
	//		Q=-Q_nom;
	//		cout<<endl<<endl<<"WARNING! Pump "<<name<<" Q="<<Qtmp;
	//		cout<<", overriden by -Q_nom="<<Q*3600<<"m3/h, H="<<H_fun(Q,n_act);
	//	}

	H_fun_lin(Q,n_act,a0,a1);

	LHS=-rho*g*a0;
	coeff_Q=rho*g*a1;
	coeff_p1=1.;
	coeff_p2=-1.;
}

void Pump::Set_BC_Left(string type, double val){
	p1=val;

	bool is_ok=false;
	if (type.compare("Pressure")==0){
		p1=val;
		double a0,a1;
		H_fun_lin(Q,n_act,a0,a1);
		Q=((p2-rho*g*a0)-val)/a1/rho/g;
		//Q=Get_Q((p2-p1)/rho/g);

		is_ok=true;
	}
	if (type.compare("FlowRate")==0){
		Q=val;
		p2=p1+H_fun(Q,n_act)*rho*g;
		is_ok=true;
	}
	if (!is_ok){
		cout<<endl<<endl<<"!!!!!!!!!!!!!! Pump::Set_BC_Left !!!!!!!!!!!!!!";
		cin.get();
	}

	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"Pump::Set_BC_Left("<<type<<","<<val<<")";
		cout<<endl<<"\t Edge "<<name<<" BC set up:";
		cout<<endl<<"\t p1   = "<<p1<<" Pa (p2="<<p2<<" Pa)";
		cout<<endl<<"\t Q    = "<<Q<<" m3/s";
		cout<<endl<<"\t H(Q) = "<<H_fun(Q,n_act)<<" mwc, (p2-p1)="<<(p2-p1)/rho/g<<" mwc"<<endl;
	}
}

void Pump::Set_BC_Right(string type, double val){
	p2=val;


	bool is_ok=false;

	if (type.compare("Pressure")==0){
		p2=val;
		double a0,a1;
		H_fun_lin(Q,n_act,a0,a1);
		Q=-((p1+rho*g*a0)-val)/a1/rho/g;
		//Q=Get_Q((p2-p1)/rho/g);
		is_ok=true;
	}
	if (type.compare("FlowRate")==0){
		Q=val;
		p2=p1+H_fun(Q,n_act)*rho*g;
		is_ok=true;
	}
	if (!is_ok){
		cout<<endl<<endl<<"!!!!!!!!!!!!!! Pump::Set_BC_Right !!!!!!!!!!!!!!";
		cin.get();
	}

	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"Pump::Set_BC_Right("<<type<<","<<val<<")";
		cout<<endl<<"\t Edge "<<name<<" BC set up:";
		cout<<endl<<"\t type:  "<<type;
		cout<<endl<<"\t val:   "<<val;
		cout<<endl<<"\t p2   = "<<p2<<" Pa = "<<p2/rho/g<<" (p1="<<p1<<" Pa)";
		cout<<endl<<"\t Q    = "<<Q<<" m3/s";
		cout<<endl<<"\t H(Q) = "<<H_fun(Q,n_act)<<" mwc, (p2-p1)="<<(p2-p1)/rho/g<<" mwc"<<endl;
	}

}

double Pump::Get_Q(double dh){
	int iter=0, iter_max=100;
	double QQ=Q, Qnew=Q, a0,a1;
	double err=fabs(dh-H_fun(QQ,n_act)), TOL=0.01;
	DEBUG=false;
	if (DEBUG){
		cout<<endl<<"PUMP "<<name<<" : entering DEBUG mode, Get_Q() iteration:";
		printf("\n\t iter %2d: Q=%5.3f, H=%5.3f, dh=%5.3f",iter,QQ,H_fun(QQ,n_act),dh);
	}
	while ((err>TOL) || (iter<3)){
		iter++;
		H_fun_lin(QQ,n_act,a0,a1);
		Qnew=(dh-a0)/a1;
		Qnew=(QQ+Qnew)/2.;
		err=fabs(dh-H_fun(Qnew,n_act));
		if (DEBUG)
			printf("\n\t iter %2d: Q=%5.3f, H=%5.3f, dh=%5.3f",iter,Qnew,H_fun(Qnew,n_act),dh);
		if (iter==iter_max){
			cout<<endl<<"ERROR!!!! Pump::Get_Q() last step:";
			printf("\n\t iter %2d: Q=%5.3f, H=%5.3f, dh=%5.3f",iter,Qnew,H_fun(Qnew,n_act),dh);
			exit(-1);
		}
		QQ=Qnew;
	}
	if (DEBUG)
		cin.get();

	return QQ;
}


void Pump::UpdateInternal(double t_target){
	//	n_act_new=n_act;
	bool ok=false;
	for (int i=0; i<transient_type.size(); i++){
		if (t>transient_time.at(i)){

			if (transient_type.at(i)=="n_new"){
				// Set n_new only once, allowing later changes.
				if (!transient_triggered.at(i)){
					cout<<endl<<"!!! Triggering transient event #"<<i<<" of pump "<<name;
					cout<<endl<<"!!! type: "<<transient_type.at(i);
					cout<<", time:"<<transient_time.at(i);
					cout<<", new val:"<<transient_val.at(i)<<"rpm"<<endl;

					transient_triggered.at(i)=true;
					n_act_new = transient_val.at(i)/60;
					n_act=n_act_new;
				}
				ok=true;
			}

			if (transient_type.at(i)=="switch_off"){
				if (!transient_triggered.at(i)){
					cout<<endl<<"!!! Triggering transient event #"<<i<<" of pump "<<name;
					cout<<endl<<"!!! type: "<<transient_type.at(i);
					cout<<", time:"<<transient_time.at(i);
					transient_triggered.at(i)=true;
				}
				double M = P_fun(Q,n_act)/(2*M_PI*n_act);
				n_act_new= n_act - (t_target-t)*M/theta;
				if (n_act_new<60.0 / 60.0)
					n_act_new= 1.0 / 60.0 ;
				n_act=n_act_new;
				//use_true_H_curve=false;

				ok=true;
			}

			if (transient_type.at(i)=="use_true_H_curve"){
				if (!transient_triggered.at(i)){
					cout<<endl<<"!!! Triggering transient event #"<<i<<" of pump "<<name;
					cout<<endl<<"!!! type: "<<transient_type.at(i);
					cout<<", time:"<<transient_time.at(i);


					transient_triggered.at(i)=true;
					use_true_H_curve=true;
				}
				ok=true;
			}
		}
	}
	t=t_target;
}


void Pump::UpdateTime(double _t){
	n_act=n_act_new;
	t=_t;
}

void Pump::Ini(int){
	//cout<<endl<<"Pump::Ini";
	if (save_data){
		tmpvec.push_back(t);
		tmpvec.push_back(p1);
		tmpvec.push_back(p2);
		tmpvec.push_back(Q);
		tmpvec.push_back(n_nom*60);
		tmpvec.push_back(n_act*60);
		tmpvec.push_back(H_fun(Q,n_act));
		tmpvec.push_back(P_fun(Q,n_act));
		data.clear();
		data.reserve(100);
		data.push_back(tmpvec);
	}
	//cout<<" done"<<endl;
};
void Pump::Ini(){
	Ini(1);
};

void Pump::Ini(double){
	Ini(1);
};


string Pump::Info(){

	std::ostringstream oss;
	oss << "\n\n Pump name     : " << name;
	oss << "\n================================";
	oss << "\n        node_from : " << node_from;
	oss << "\n          node_to : " << node_to;
	oss << "\n              rho : " << rho<<" kg/m3";
	oss << "\n                t : " << t<<" s";
	oss << "\n      Q (current) : " << Q<<" m3/s = "<<Q*3600<<" m3/h";
	oss << "\n      n (current) : " << n_act*60<<" rpm";
	oss << "\n      H (current) : " << H_fun(Q,n_act)<<" mwc";
	oss << "\n     p1 (current) : " << p1<<" Pa = "<<p1/rho/g<<" mwc";
	oss << "\n     p2 (current) : " << p2<<" Pa = "<<p2/rho/g<<" mwc"<<endl;
	oss << "\n            Q nom : " << Q_nom<<" m3/s";
	oss << "\n H(Q nom) @ n_nom : " << H_fun(Q_nom,n_nom)<<" mwc";
	oss << "\n P(Q nom) @ n_nom : " << P_fun(Q_nom,n_nom)/1000.<<" kW";
	oss << "\n     H(0) @ n_nom : " << H_fun(0,n_nom) << " mwc";
	oss << "\n     P(0) @ n_nom : " << P_fun(0,n_nom)/1000. << " kW";
	oss << "\n           H(Q) : ";
	oss<< coeff_H.at(0)<<" + ";
	oss<< coeff_H.at(1)<<" Q ";
	for (unsigned int i=2; i<coeff_H.size(); i++)
		oss<<" + "<< coeff_H.at(i)<<" Q^"<<i;
	oss << "\n           P(Q) : ";
	oss<< coeff_P.at(0)<<" + ";
	oss<< coeff_P.at(1)<<" Q ";
	for (unsigned int i=2; i<coeff_P.size(); i++)
		oss<<" + "<< coeff_P.at(i)<<" Q^"<<i;
	oss << "\n          n_nom : "<<n_nom*60<<" rpm";
	oss << "\n          theta : "<<theta<<" kg m2";
	oss << "\n           Tmax : "<<Tmax<<" s (estimated time to stop from n_nom, Q_nom)";
	oss << "\n             dt : "<<dt<<" s";
	oss << endl;
	return oss.str();
}

void Pump::GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void Pump::GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&){};

void Pump::Save_data(){

	tmpvec.at(0) = t;
	tmpvec.at(1) = p1;
	tmpvec.at(2) = p2;
	tmpvec.at(3) = Q;
	tmpvec.at(4) = n_nom*60;
	tmpvec.at(5) = n_act*60;
	tmpvec.at(6) = H_fun(Q,n_act);
	tmpvec.at(7) = P_fun(Q,n_act);

	data.push_back(tmpvec);

}



void Pump::Write_data(string folder){
	
	if (!save_data) {
		cout << endl << "WARNING! SCP: " << name;
		cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
	}
	else {
		cout << endl << "Saving to " << fname.c_str() << " ... ";
		FILE * pFile;
		pFile = fopen ((folder + "/" + fname).c_str(), "w");
		fprintf(pFile, "t (s); p1 (bar); p2 (bar); Q m3/h; n_nom (rpm); n_act (rpm); H (mwc); P (kW)\n");
		for (int i = 0; i < data.size(); i++)
			fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e;\n",
					data.at(i).at(0),
					data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
					data.at(i).at(3)*3600,
					data.at(i).at(4), data.at(i).at(5),
					data.at(i).at(6), data.at(i).at(7)/1000.
					);
		fclose (pFile);
		cout << " done. ";
	}
}

void Pump::Add_transient(string _type, double _when, double _val){
	transient_type.push_back(_type);
	transient_time.push_back(_when);
	transient_val.push_back(_val);
	transient_triggered.push_back(false);
}

double Pump::Get_dprop(string){return 0.;};
