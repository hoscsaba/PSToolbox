#define USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>
#include <string>
#include "my_tools.h"
#include "LWP.h"
#include "SCP.h"
#include "Valve.h"
#include "Reservoir.h"
#include "Damper.h"
#include "Connector.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>

// ! Default Constructor
Connector::Connector(bool _DEBUG): demand_k(0), demand_p0(0), demand_t0(0) {
	DEBUG = _DEBUG;
	type = -1;
	demand = 0;
	BC_type = "none";
	BC_value = 0;
	is_front1 = false;
	is_front2 = false;
	is_front3 = false;
	is_front4 = false;
	e1 = nullptr;
	e2 = nullptr;
	e3 = nullptr;
	e4 = nullptr;

	p_prev = 0.;
}

// ! Empty Destructor
Connector::~Connector() = default;

Connector::Connector(
	string _name,
	vector<PSToolboxBaseEdge *> &_e,
	vector<bool> &_is_front,
	vector<int> &_edges_idx,
	double _demand, bool _DEBUG): demand_k(0), demand_p0(0), demand_t0(0) {
	name = std::move(_name);
	edges = _e;
	is_front = _is_front;
	edges_idx = _edges_idx;
	demand = _demand;
	DEBUG = _DEBUG;
	type = 0;

	BC_type = "none";
	BC_value = 0;
	is_front1 = false;
	is_front2 = false;
	is_front3 = false;
	is_front4 = false;
	e1 = nullptr;
	e2 = nullptr;
	e3 = nullptr;
	e4 = nullptr;

	p_prev = 0.;
}

// ! Edges junction for two edges
Connector::Connector(
	string _name,
	PSToolboxBaseEdge *_e1, bool _is_front1,
	PSToolboxBaseEdge *_e2, bool _is_front2,
	double _demand, bool _DEBUG, vector<int> &_edges_idx) {
	name = std::move(_name);
	e1 = _e1;
	e2 = _e2;
	e3 = nullptr;
	e4 = nullptr;
	is_front1 = _is_front1;
	is_front2 = _is_front2;
	is_front3 = false;
	is_front4 = false;
	BC_type = "none";
	BC_value = 0;
	demand = _demand;
	DEBUG = _DEBUG;
	edges_idx = _edges_idx;
	type = 2;

	p_prev = 0.;
}

// ! Edge-to simple BC
Connector::Connector(
	string _name,
	PSToolboxBaseEdge *_e1, bool _is_front1,
	string _BC_type, double _BC_value,
	double _demand, bool _DEBUG) {
	name = _name;
	e1 = _e1;
	e2 = NULL;
	e3 = NULL;
	e4 = NULL;
	is_front1 = _is_front1;
	is_front2 = false;
	is_front3 = false;
	is_front4 = false;
	BC_type = _BC_type;
	BC_value = _BC_value;
	demand = _demand;
	DEBUG = _DEBUG;

	type = 1;
	p_prev = 1.e5;
}

// ! 3 SCP pipes
Connector::Connector(
	string _name,
	PSToolboxBaseEdge *_e1, bool _is_front1,
	PSToolboxBaseEdge *_e2, bool _is_front2,
	PSToolboxBaseEdge *_e3, bool _is_front3,
	double _demand, bool _DEBUG, vector<int> &_edges_idx) {
	name = _name;
	e1 = _e1;
	e2 = _e2;
	e3 = _e3;
	e4 = NULL;
	is_front1 = _is_front1;
	is_front2 = _is_front2;
	is_front3 = _is_front3;
	is_front4 = false;
	BC_type = "none";
	BC_value = 0;
	demand = _demand;
	DEBUG = _DEBUG;
	edges_idx = _edges_idx;

	/*cout << "Internal size: " << edges_idx.size() << endl;
    cout << "External size: " << _edges_idx.size() << endl;*/


	type = 3;

	p_prev = 0.;
}

// ! 4 SCP pipes
Connector::Connector(
	string _name,
	PSToolboxBaseEdge *_e1, bool _is_front1,
	PSToolboxBaseEdge *_e2, bool _is_front2,
	PSToolboxBaseEdge *_e3, bool _is_front3,
	PSToolboxBaseEdge *_e4, bool _is_front4,
	double _demand, bool _DEBUG, vector<int> &_edges_idx) {
	name = _name;
	e1 = _e1;
	e2 = _e2;
	e3 = _e3;
	e4 = _e4;
	is_front1 = _is_front1;
	is_front2 = _is_front2;
	is_front3 = _is_front3;
	is_front4 = _is_front4;
	demand = _demand;
	DEBUG = _DEBUG;
	edges_idx = _edges_idx;

	BC_type = "none";
	BC_value = 0;

	type = 4;
	p_prev = 0.;
}


// This is the main update function
void Connector::Update(double t_target, int update_idx) {
	/*if(name == "30_2")
	{
	 cout << "Con " << name << "updated" << endl;
	}*/
	switch (type) {
		case 0:
			Connector_SCP_NPipes(t_target, update_idx);
			break;

		case 1:
			Connector_SCP_Pipe_Simple_BC(t_target);
			break;

		case 2:
			Connector_SCP_2Pipes(t_target, update_idx);
			break;

		case 3:
			Connector_SCP_3Pipes(t_target, update_idx);
			break;

		case 4:
			Connector_SCP_4Pipes(t_target, update_idx);
			break;

		default:
			cout << endl << "ERROR!!!";
			cout << endl << "Connector::Update() -> unknown type=" << type << endl;
			cin.get();
			break;
	}
}

//! Connect Reservoir to Valve
/*! Connect Reservoir to Valve
  \param t_target, s
  \param Reservoir* reservoir
  \param Valve* valve
  \param p_downstream
  \param &p pressure @ connection
  \param &T temperature @ connection
 */

void Connector::Connector_Reservoir_and_Valve_Inlet(double t_target,
                                                    Reservoir *r1, Valve *v1, bool INLET_PRESSURE_DROP,
                                                    double p_downstream, double &p, double &mp) {
	// Solve the following system for p,T,rho,v,a:
	// (1) mp = MassFlowRate(pt,Tt,pd,Tout,x)
	// (2) rho=rho(p,T)gt
	// (3) a  =a(p,T)
	// (4a) T/Tr=(pr/p)^(kappa-1)/kappa (gas)
	// (4b) pr=p+rho/2*v^2 (liquid)

	double pr = r1->Get_dprop("p");
	double Tr = r1->Get_dprop("T");
	double Abore = v1->Get_dprop("A");

	double rho, v, T, p_new, cp, cV, kappa;

	double err = 1., TOL = 1e-3;
	int iter = 0, MAX_ITER = 10;
	p = r1->Get_dprop("p");
	if (r1->is_Gas)
		T = r1->Get_dprop("T");
	else
		T = 293.;

	if (DEBUG)
		cout << endl << "Connector::Connector_Reservoir_and_Pipe_Front -> entering debug mode" << endl;

	while ((err > TOL) && (iter <= MAX_ITER)) {
		rho = r1->Get_dprop("rho");
		mp = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
		v = mp / rho / Abore;
		if (r1->is_Gas) {
			cp = r1->gas->Get_cp(p, T);
			cV = r1->gas->Get_cV(p, T);
			kappa = cp / cV;
			T = r1->gas->Get_T(rho, p);
			if (INLET_PRESSURE_DROP)
				p_new = pr / pow(T / Tr, kappa / (kappa - 1));
			else
				p_new = pr;
		} else {
			if (INLET_PRESSURE_DROP)
				p_new = pr - rho / 2. * v * v;
			else p_new = pr;
		}
		if (DEBUG)
			printf("\n iter: %2d, p=%5.3e, rho=%5.3e, mp=%5.3e, v=%5.3e, p_new=%5.3e, err=%5.3e",
			       iter, p, rho, mp, v, p_new, err);

		err = fabs(p - p_new) / 1.e5;
		p = p_new;
		if (iter == MAX_ITER) {
			cout << endl << "Connector.Connector_Reservoir_and_Valve_Inlet() ERROR!!!" << endl;
			cin.get();
		}
	}
}

//! Connect LWP end to Valve
/*! Connect LWP end to Valve
  \param t_target, s
  \param LWP* pipe
  \param Valve* valve
  \param p_downstream
  \param &p pressure @ connection
  \param &T temperature @ connection
 */

void Connector::Connector_LWP_Pipe_Back_and_Valve(double t_target,
                                                  LWP *p1, Valve *v1, double p_downstream) {
	// Solve the following system for p,T,rho,v:
	// (1) alpha=p+ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow
	// (3) T=TL
	// (4) rho=rho(p,T)
	// (5) a  =a(T)

	double alpha, Apipe;

	Apipe = p1->Get_dprop("A");
	alpha = p1->GetAlphaPrimitiveAtEnd(t_target);
	double err1 = 1.e5, err2 = 1.e5, mp;
	int iter = 0, MAX_ITER = 100;
	double p = p1->Get_dprop("p_back");
	double T = p1->Get_dprop("T_back");
	double TOL_p = p / 100., TOL_T = T / 100.;
	while ((err1 > TOL_p) || (err2 > TOL_T)) {
		mp = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));

		double a = p1->gas->Get_SonicVel(T, p);
		double f = alpha - (p + mp / Apipe * a);

		double dp = 0.001 * p;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp = v1->Get_MassFlowRate(p + dp, T, p_downstream, 293., v1->Get_dprop("x"));

		a = p1->gas->Get_SonicVel(T, p + dp);
		double f1 = alpha - (p + dp + mp / Apipe * a);
		double df_dp = (f1 - f) / dp;
		double p_new = p - f / df_dp;
		double T_new = p1->gas->Get_T(p_new, p1->gas->Get_rho(p, T));

		T = 0.2 * T + 0.8 * T_new;
		p = 0.2 * p + 0.8 * p_new;

		err1 = fabs(f);
		err2 = fabs(T - T_new);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, p, T, a, mp, err1, err2);
			cin.get();
		}

		if (DEBUG)
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, p_downstream=%5.3e, T=%5.3e, a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, p, p_downstream, T, a, mp, err1, err2);
	}
	//cout << endl << "pipe " << p1->Get_name() << ": mp=" << mp * 1000 << " g/s";
	p1->BCRight(t_target, "StaticPres_Outlet", p, T, true);
}

void Connector::Connector_LWP_Pipe_Back_Valve_Front(double t_target,
                                                    LWP *pipe_u, Valve *v1, LWP *pipe_d) {
	// Solve the following system for mp, v1,v2, p1, p2
	// (1) alpha=p1+ro1*a1*v1 = p1+mp/Apipe1*a1
	// (2) beta =p2-ro2*a2*v2 = p2-mp/Apipe2*a2
	// (3) ro1*v1*Apipe1=valve_mass_flow
	// (4) ro2*v2*Apipe2=valve_mass_flow
	// (5) cp*T1 + v1^2/2 = cp*T2 + v2^2/2

	double Apipe1 = pipe_u->Get_dprop("A");
	double Apipe2 = pipe_d->Get_dprop("A");
	double alpha = pipe_u->GetAlphaPrimitiveAtEnd(t_target);
	double beta = pipe_d->GetBetaPrimitiveAtFront(t_target);

	double pu = pipe_u->Get_dprop("p_back");
	double Tu = pipe_u->Get_dprop("T_back");
	double rhou = pipe_u->Get_dprop("rho_back");

	double pd = pipe_d->Get_dprop("p_front");
	double Td = pipe_d->Get_dprop("T_front");
	double rhod = pipe_d->Get_dprop("rho_front");

	double err1 = 1.e5, err2 = 1.e5, mp;
	int iter = 0, MAX_ITER = 100;

	double TOL_p = 100., TOL_T = 0.1;
	while ((err1 > TOL_p) || (err2 > TOL_T)) {
		mp = v1->Get_MassFlowRate(pu, Tu, pd, Td, v1->Get_dprop("x"));

		double a = pipe_u->gas->Get_SonicVel(Tu, pu);
		double f = alpha - (pu + mp / Apipe1 * a);

		double dp = 0.001 * pu;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp = v1->Get_MassFlowRate(pu + dp, Tu, pd, Td, v1->Get_dprop("x"));

		a = pipe_u->gas->Get_SonicVel(Tu, pu + dp);
		double f1 = alpha - (pu + dp + mp / Apipe1 * a);
		double df_dp = (f1 - f) / dp;
		double p_new = pu - f / df_dp;
		double T_new = pipe_u->gas->Get_T(p_new, pipe_u->gas->Get_rho(pu, Tu));

		Tu = 0.2 * Tu + 0.8 * T_new;
		pu = 0.2 * pu + 0.8 * p_new;

		err1 = fabs(f);
		err2 = fabs(Tu - T_new);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, pu, Tu, a, mp, err1, err2);
			cin.get();
		}

		if (DEBUG)
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, pu=%5.3e, pd=%5.3e, Tu=%5.3e, a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, pu, pd, Tu, a, mp, err1, err2);
	}
	//cout << endl << "pipe " << p1->Get_name() << ": mp=" << mp * 1000 << " g/s";
	pipe_u->BCRight(t_target, "StaticPres_Outlet", pu, Tu, true);
}

void Connector::Update(double t_target, LWP *pipe_u, Damper *damper, LWP *pipe_d) {
	// Solve the following system for mp, v1,v2, p1, p2
	// (1) alpha=p+mp_u/A_u*a_u
	// (2) beta =p+mp_d/A_d*a_d
	// (3) mp_u=mp_d+mp_damper
	// (4) mp_damper = Damper::GetMassFlowRate(p)

	double A_u = pipe_u->Get_dprop("A");
	double A_d = pipe_d->Get_dprop("A");
	double alpha = pipe_u->GetAlphaPrimitiveAtEnd(t_target);
	double beta = pipe_d->GetBetaPrimitiveAtFront(t_target);

	double p_u = pipe_u->Get_dprop("p_back");
	double T_u = pipe_u->Get_dprop("T_back");
	double rho_u = pipe_u->Get_dprop("rho_back");

	double p_d = pipe_d->Get_dprop("p_front");
	double T_d = pipe_d->Get_dprop("T_front");
	double rho_d = pipe_d->Get_dprop("rho_front");

	double err = 1.e5;
	int iter = 0, MAX_ITER = 100;

	double TOL_p = 0.001e5;
	double p = (p_u + p_d) / 2.;
	while (err > TOL_p) {
		double mp_damper = damper->GetMassFlowRate(p);
		double mp_u = (alpha - p) * A_u / pipe_u->gas->Get_SonicVel(T_u, p);
		double mp_d = -(beta - p) * A_d / pipe_d->gas->Get_SonicVel(T_d, p);
		double f = mp_u - (mp_d + mp_damper);

		double dp = 0.001 * p;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp_damper = damper->GetMassFlowRate(p + dp);
		mp_u = (alpha - (p + dp)) * A_u / pipe_u->gas->Get_SonicVel(T_u, p + dp);
		mp_d = -(beta - (p + dp)) * A_d / pipe_d->gas->Get_SonicVel(T_d, p + dp);
		double f1 = mp_u - (mp_d + mp_damper);
		double df_dp = (f1 - f) / dp;
		double p_new = p - f / df_dp;

		p = 0.2 * p + 0.8 * p_new;

		err = fabs(f);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Update(t,LWP,Damper,LWP) -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: p=%5.3e, mp_u=%5.3e, mp_d=%5.3e, mp_damper=%5.3e er1=%+5.3e",
				iter, p, mp_u, mp_d, mp_damper, err);
			cin.get();
		}

		if (DEBUG)
			printf(
				"\n (back) iter #%2d: p=%5.3e, mp_u=%5.3e, mp_d=%5.3e, mp_damper=%5.3e er1=%+5.3e",
				iter, p, mp_u, mp_d, mp_damper, err);
	}
	double mp_u = (alpha - p) * A_u / pipe_u->gas->Get_SonicVel(T_u, p);
	double mp_d = -(beta - p) * A_d / pipe_d->gas->Get_SonicVel(T_d, p);
	if (mp_u > 0.)
		pipe_u->BCRight(t_target, "StaticPres_Outlet", p, T_u, true);
	else
		pipe_u->BCRight(t_target, "StaticPres_and_StaticTemp_Inlet", p, T_d, true);
	if (mp_d < 0.)
		pipe_u->BCLeft(t_target, "StaticPres_Outlet", p, T_d, true);
	else
		pipe_u->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", p, T_u, true);
}

void Connector::Update(double t_target, LWP *pipe_u, Damper *damper, Valve *valve, double p_d, double T_d) {
	// Solve the following system for mp, v1,v2, p1, p2
	// (1) alpha=p+mp_u/A_u*a_u
	// (2) mp_prv=Valve::GetMassFlowRate()
	// (3) mp_u=mp_d+mp_damper
	// (4) mp_damper = Damper::GetMassFlowRate(p)

	double A_u = pipe_u->Get_dprop("A");
	double alpha = pipe_u->GetAlphaPrimitiveAtEnd(t_target);
	double p_u = pipe_u->Get_dprop("p_back");
	double rho_u = pipe_u->Get_dprop("rho_back");
	double T_u = pipe_u->Get_dprop("T_back");

	double err = 1.e5;
	int iter = 0, MAX_ITER = 100;

	double TOL_mp = rho_u * A_u * 0.01;
	double p = p_u;
	while (err > TOL_mp) {
		double mp_damper = damper->GetMassFlowRate(p);
		double mp_u = (alpha - p) * A_u / pipe_u->gas->Get_SonicVel(T_u, p);
		double mp_d = valve->Get_MassFlowRate(p, T_u, p_d, T_d, valve->Get_dprop("x"));
		double f = mp_u - (mp_d + mp_damper);

		double dp = 0.001 * (p - damper->Get_dprop("p"));
		if (fabs(dp) < 1.)
			dp = 1. * (dp < 0. ? -1. : 1.);

		mp_damper = damper->GetMassFlowRate(p + dp);
		mp_u = (alpha - (p + dp)) * A_u / pipe_u->gas->Get_SonicVel(T_u, p + dp);
		mp_d = valve->Get_MassFlowRate(p + dp, T_u, p_d, T_d, valve->Get_dprop("x"));
		double f1 = mp_u - (mp_d + mp_damper);
		double df_dp = (f1 - f) / dp;
		double p_new = p - f / df_dp;

		p = 0.5 * p + 0.5 * p_new;

		err = fabs(f);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Update(t,LWP,Damper,LWP) -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (LWP back - Damper - Valve) iter #%2d: p=%5.3e, dp=%+5.3e, mp_u=%5.3e, mp_d=%5.3e, mp_damper=%5.3e err/TOL=%+5.3e",
				iter, p, p - damper->Get_dprop("p"), mp_u, mp_d, mp_damper, err / TOL_mp);
			cin.get();
		}

		if (DEBUG)
			printf(
				"\n (LWP back - Damper - Valve) iter #%2d: p=%5.3e, dp=%+5.3e, mp_u=%5.3e, mp_d=%5.3e, mp_damper=%5.3e err/TOL=%+5.3e",
				iter, p, p - damper->Get_dprop("p"), mp_u, mp_d, mp_damper, err / TOL_mp);
	}
	double mp_u = (alpha - p) * A_u / pipe_u->gas->Get_SonicVel(T_u, p);
	if (mp_u > 0.)
		pipe_u->BCRight(t_target, "StaticPres_Outlet", p, T_u, true);
	else
		pipe_u->BCRight(t_target, "StaticPres_and_StaticTemp_Inlet", p, T_d, true);
}


/*
void Connector::Connector_LWP_Valve_Pipe_Front(double t_target,
                                               LWP *p1, Valve *v1, double p_upstream) {
	// Solve the following system for p,T,rho,v:
	// (1) alpha=p+ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow
	// (3) T=TL
	// (4) rho=rho(p,T)
	// (5) a  =a(T)

	//double kappa=p1->gas->kappa;
	//double cp   =p1->gas->cp;

	double update_OK = true;
	double rho, v, a, alpha, Apipe;

	Apipe = p1->Get_dprop("A");
	//cout<<endl<<" -------> 1. Connector::Connector_LWP_Pipe_Back_and_Valve() t_target="<<t_target<<endl;
	alpha = p1->GetAlphaPrimitiveAtEnd(t_target);
	//cout<<endl<<" -------> 1. Connector::Connector_LWP_Pipe_Back_and_Valve() done."<<endl;
	double err1 = 0., err2, err = 1.e5, mp, TOL = 1e-5;
	int iter = 0, MAX_ITER = 10;
	double p = p1->Get_dprop("p_back");
	double T = p1->Get_dprop("T_back");
	while ((err > TOL) && (iter <= MAX_ITER)) {
		a = p1->gas->Get_SonicVel(T, p);
		rho = p1->gas->Get_rho(p, T);
		mp = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
		v = mp / rho / Apipe;
		p = alpha - rho * a * v;
		T = p1->gas->Get_T(rho, p);
		//err1 = (alpha - (p + rho * a * v)) / 1.e5;
		err2 = rho * v * Apipe - mp;
		//err = sqrt(err1 * err1 + err2 * err2);
		err = sqrt(err2 * err2);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, p, T, rho, a, v, err1, err2);
			//cin.get();
			update_OK = false;
		}
	}
	//cout<<endl<<" -------> 2. Connector::Connector_LWP_Pipe_Back_and_Valve() t_target="<<t_target<<endl;
	p1->BCRight(t_target, "StaticPres_Outlet", p, T, true);
	//cout<<endl<<" -------> 2. Connector::Connector_LWP_Pipe_Back_and_Valve() done."<<endl;
	//v1->Update(t_target,p,p_downstream);

	//return update_OK;
}*/

void Connector::Connector_LWP_Pipe_Front_and_Valve(double t_target,
                                                   LWP *p1, Valve *v1, double p_upstream, double T_upstream) {
	// Solve the following system for p,T,rho,v:
	// (1) beta=p-ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow
	// (3) T=TL
	// (4) rho=rho(p,T)
	// (5) a  =a(T)

	double p = p1->Get_dprop("p_front");
	double T = p1->Get_dprop("T_front");
	double Apipe = p1->Get_dprop("A");
	double beta = p1->GetBetaPrimitiveAtFront(t_target);

	double err1 = 1.e5, err2 = 1.e5, TOL_p_Pa = fabs(p) / 1000., TOL_T = fabs(T) / 1000., mp;
	int iter = 0, MAX_ITER = 100;
	while ((err1 > TOL_p_Pa) || (err2 > TOL_T)) {
		/*double mp = v1->Get_MassFlowRate(p_upstream, T_upstream, p, T, v1->Get_dprop("x"));
		double a = p1->gas->Get_SonicVel(T, p);
		double rho = p1->gas->Get_rho(p, T);
		double v = mp / rho / Apipe;
		double p_new = beta + rho * a * v;
		double T_new = p1->gas->Get_T(p_new, rho);
		//err = sqrt(err1 * err1 + err2 * err2);
		T = 0.1 * T + 0.9 * T_new;
		p = 0.1 * p + 0.9 * p_new;
		err1 = fabs(beta - (p - rho * a * v));
		err2 = fabs(T - T_new);*/

		mp = v1->Get_MassFlowRate(p_upstream, T_upstream, p, T, v1->Get_dprop("x"));

		double a = p1->gas->Get_SonicVel(T, p);
		double f = beta - (p - mp / Apipe * a);

		double dp = 0.001 * p;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp = v1->Get_MassFlowRate(p_upstream, T_upstream, p + dp, T, v1->Get_dprop("x"));

		a = p1->gas->Get_SonicVel(T, p + dp);

		double f1 = beta - (p + dp - mp / Apipe * a);
		double df_dp = (f1 - f) / dp;
		double p_new = p - f / df_dp;
		double T_new = p1->gas->Get_T(p_new, p1->gas->Get_rho(p, T));

		T = 0.6 * T + 0.4 * T_new;
		p = 0.6 * p + 0.4 * p_new;
		err1 = fabs(f);
		err2 = fabs(T - T_new);

		iter++;

		if (DEBUG)
			printf(
				"\n (front) iter #%2d: beta=%5.3e, p_upstream=%5.3e, p=%5.3e, T_upstream=%5.3e, T=%5.3e,a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, beta, p_upstream, p, T_upstream, T, a, mp, err1, err2);

		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Pipe_Front_and_Valve() -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: beta=%5.3e, p=%5.3e, T=%5.3e, a=%5.3e, mp=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, beta, p, T, a, mp, err1, err2);
			cin.get();
		}
	}
	//cout << endl << "pipe " << p1->Get_name() << ": mp=" << mp * 1000 << " g/s";
	p1->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", p, T, true);
	//cin.get();
}


bool Connector::Connector_LWP_Pipe_Back_and_Valve_with_Absorber(double t_target,
                                                                LWP *p1, Valve_with_Absorber *v1, double p_downstream,
                                                                double &p, double &T) {
	// Solve the following system for p,T,rho,v:
	// (1) alpha=p+ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow
	// (3) T=TL
	// (4) rho=rho(p,T)
	// (5) a  =a(T)

	//double kappa=p1->gas->kappa;
	//double cp   =p1->gas->cp;

	double update_OK = true;
	double rho, v, a, alpha, Apipe;

	Apipe = p1->Get_dprop("A");
	//p1->GetAllPrimitiveAtEnd(t_target, pL, v, T, rho);
	alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	double err1, err2, err = 1.e5, mp, TOL = 1e-5;
	int iter = 0, MAX_ITER = 10;
	p = p1->Get_dprop("p_back");
	while ((err > TOL) && (iter <= MAX_ITER)) {
		a = p1->gas->Get_SonicVel(T, p);
		rho = p1->gas->Get_rho(p, T);
		mp = v1->Get_MassFlowRate(p, T, p_downstream, 293., v1->Get_dprop("x"));
		v = mp / rho / Apipe;
		p = alpha - rho * a * v;
		//T   = p1->gas->Get_T(rho,p);
		err1 = (alpha - (p + rho * a * v)) / 1.e5;
		err2 = rho * v * Apipe - mp;
		err = sqrt(err1 * err1 + err2 * err2);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Pipe_Back_and_Valve() -> MAX_ITER reached, stopping iteration.";
			printf(
				"\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
				iter, alpha, p, T, rho, a, v, err1, err2);
			//cin.get();
			update_OK = false;
		}
	}

	//cout<<endl<<"\n\n exiting Connector...";
	//printf("\n\t p=%5.3e, T=%5.3e, rho=%5.3e, a=%5.3e, v=%5.3e, mp=%+5.3e",p, T, rho, a, v, mp);
	return update_OK;
}


void Connector::Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
                                                       Reservoir *r1, LWP *p1, bool inlet_pressure_drop) {
	//cout << endl << "Entering Connector_LWP_Reservoir_and_Pipe_Front()..." << endl;

	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");

	bool is_inflow = false, is_outflow = false;
	if (!inlet_pressure_drop) {
		is_inflow = p1->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", pt, Tt, false);

		if (DEBUG) {
			// cout << endl << "Connector_LWP_Reservoir_and_Pipe_Front()";
			// cout << endl << "pt=" << pt << ", Tt=" << Tt;
			// cout << endl << "is_inflow = " << is_inflow;
			//cin.get();
		}
	} else {
		is_inflow = Connector_LWP_Reservoir_and_Pipe_Front_inlet(t_target, r1, p1, inlet_pressure_drop, pt, Tt);

		// Old pt, Tt values were overwritten load them again
		pt = r1->Get_dprop("p");
		Tt = r1->Get_dprop("T");
	}

	if (is_inflow)
		is_inflow = p1->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", pt, Tt, true);
	else {
		is_outflow = p1->BCLeft(t_target, "StaticPres_Outlet", pt, 0, false);

		if (is_outflow)
			is_outflow = p1->BCLeft(t_target, "StaticPres_Outlet", pt, 0, true);
		else {
			p1->BCLeft(t_target, "Wall", 0, 0, true);
			cout << endl <<
					"WARNING: Connector_LWP_Reservoir_and_Pipe_Front() -> Neither inflow, nor outflow! Applying wall BC.";
			if (DEBUG)
				cin.get();
		}
	}
	//cout << endl << "Leaving Connector_LWP_Reservoir_and_Pipe_Front()..." << endl;
}


void Connector::Connector_LWP_Reservoir_and_Pipe_Back(double t_target,
                                                      Reservoir *r1, LWP *p1, bool inlet_pressure_drop) {
	const double pt = r1->Get_dprop("p");
	const double Tt = r1->Get_dprop("T");

	//cout << endl << endl << " ****** Tt = " << Tt << " ****** " << endl;

	if (inlet_pressure_drop) {
		cout << endl << endl << "WARNING!!!<<endl";
		cout <<
				"Connector::Connector_LWP_Reservoir_and_Pipe_Back inlet_pressure_drop = true not implemented, proceeding with false.";
		cout << endl << endl;
	}

	// Trying inlet flow BC
	bool is_outflow = false;
	bool is_inflow = p1->BCRight(t_target, "StaticPres_and_StaticTemp_Inlet", pt, Tt, false);

	if (is_inflow)
		is_inflow = p1->BCRight(t_target, "StaticPres_and_StaticTemp_Inlet", pt, Tt, true);
	else {
		// try outflow
		is_outflow = p1->BCRight(t_target, "StaticPres_Outlet", pt, 0., false);

		if (is_outflow)
			is_outflow = p1->BCRight(t_target, "StaticPres_Outlet", pt, 0., true);
		else {
			p1->BCRight(t_target, "Wall", 0, 0, true);
			cout << endl <<
					"WARNING: Connector_LWP_Reservoir_and_Pipe_Back() -> Neither inflow, nor outflow! Applying wall BC.";
			cout << endl << "   Nem of the pipe: " << p1->Get_name() << endl;
			if (DEBUG)
				cin.get();
		}
	}
	if (DEBUG) {
		cout << endl << "Connector_LWP_Reservoir_and_Pipe_Back()";
		cout << endl << "pt=" << pt << ", Tt=" << Tt;
		cout << endl << "p1=" << p1->Get_dprop("p_back");
		//cout << ", v1=" << p1->Get_dprop("v_back") << ", dp=" << p1->Get_dprop("p_back") - pt;
		cout << endl << "is_inflow = " << is_inflow;
		cout << endl << "is_outflow = " << is_outflow << endl;
		//cin.get();
	}
}

bool Connector::Connector_LWP_Reservoir_and_Pipe_Front_inlet(double t_target,
                                                             Reservoir *r1, LWP *p1, bool inlet_pressure_drop,
                                                             double &p, double &T) {
	// Solve the following system for p,T,rho,v:
	// (1)  Tt=T+v^2/2/cp            Isentropic flow from the reservoir to the pipe
	// (2)  p-rho*a*v = beta_front   Charactersistic equation
	// (3)  rho = rho(p,T)
	// (4a) p=pt;
	// (4b) pt/rhot^k = p/rho^k

	if (DEBUG)
		cout << endl << endl << "Connector::Connector_LWP_Reservoir_and_Pipe_Front: entering DEBUG mode" <<
				endl;

	double cp = 1000., kappa_Tv = 1.4;
	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");
	//double rhot = p1->gas->Get_rho(pt, Tt);

	double rho, v, a, beta;
	//p1->GetAllPrimitiveAtFront(t_target, p, v, T, rho);
	beta = p1->GetBetaPrimitiveAtFront(t_target);
	//bool is_ok = false;

	// Assume inflow
	v = p1->Get_dprop("v_front");
	if (v < 0)
		v = 0.;

	double err = 1.e5, TOL = 0.001 /*m/s*/, dv, f, f1, df, vnew;
	int iter = 0, MAX_ITER = 50;
	while (fabs(err) > TOL) {
		f = Connector_LWP_Reservoir_and_Pipe_Front_fun(v, beta, r1, p1, inlet_pressure_drop);

		if (fabs(v) > 0.001)
			dv = v * 0.01;
		else
			dv = 0.001;
		f1 = Connector_LWP_Reservoir_and_Pipe_Front_fun(v + dv, beta, r1, p1, inlet_pressure_drop);

		df = (f1 - f) / dv;
		vnew = v - f / df;

		err = f;

		if (DEBUG)
			printf("\n (front, inflow) inner iter #%2d: pt=%5.3e, beta=%5.3e, v=%5.3e, err=%5.3e", iter, pt,
			       beta, v,
			       err);

		iter++;

		double RELAX = 0.5;
		v = (1. - RELAX) * v + RELAX * vnew;

		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (1)";
			exit(-1);
		}
	}

	double Tnew;
	iter = 0;
	err = 1.e5;
	T = 293;
	while (fabs(err) > 0.1) {
		cp = p1->gas->Get_cp(p, T);
		kappa_Tv = p1->gas->Get_kappa_Tv();
		Tnew = Tt - v * v / 2 / cp;
		if (inlet_pressure_drop)
			p = pt * pow(Tnew / Tt, kappa_Tv / (kappa_Tv - 1.));
		else
			p = pt;
		rho = p1->gas->Get_rho(p, Tnew);
		a = p1->gas->Get_SonicVel(Tnew, p);

		err = T - Tnew;
		T = Tnew;

		if (DEBUG)
			printf(
				"\n (front, inflow) outer iter #%2d: cp=%5.3e, kappa_Tv=%5.3e, p=%5.3e, T=%5.3e, v=%5.3e, err=%5.3e",
				iter, cp, kappa_Tv, p, T, v, err);

		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (2)";
			exit(-1);
		}

		iter++;
	}

	bool is_inflow = false;
	if (v > 0)
		is_inflow = true;
	return is_inflow;
}

void Connector::Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
                                                       Reservoir *r1, LWP *p1, bool inlet_pressure_drop,
                                                       double &p, double &T) {
	// Solve the following system for p,T,rho,v:
	// (1)  Tt=T+v^2/2/cp            Isentropic flow from the reservoir to the pipe
	// (2)  p-rho*a*v = beta_front   Charactersistic equation
	// (3)  rho = rho(p,T)
	// (4a) p=pt;
	// (4b) pt/rhot^k = p/rho^k

	if (DEBUG)
		cout << endl << endl << "Connector::Connector_LWP_Reservoir_and_Pipe_Front: entering DEBUG mode" <<
				endl;

	double cp = 1000., kappa_Tv = 1.4;
	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");
	//double rhot = p1->gas->Get_rho(pt, Tt);

	double rho, v, a, beta;
	//p1->GetAllPrimitiveAtFront(t_target, p, v, T, rho);
	beta = p1->GetBetaPrimitiveAtFront(t_target);
	bool is_ok = false;

	// Assume inflow
	v = 1.;

	double err = 1.e5, TOL = 0.001 /*m/s*/, dv, f, f1, df, vnew;
	int iter = 0, MAX_ITER = 50;
	while (fabs(err) > TOL) {
		f = Connector_LWP_Reservoir_and_Pipe_Front_fun(v, beta, r1, p1, inlet_pressure_drop);

		if (fabs(v) > 0.001)
			dv = v * 0.01;
		else
			dv = 0.001;
		f1 = Connector_LWP_Reservoir_and_Pipe_Front_fun(v + dv, beta, r1, p1, inlet_pressure_drop);

		df = (f1 - f) / dv;
		vnew = v - f / df;

		err = f;

		if (DEBUG)
			printf("\n (front, inflow) inner iter #%2d: pt=%5.3e, beta=%5.3e, v=%5.3e, err=%5.3e", iter, pt,
			       beta, v,
			       err);

		iter++;

		double RELAX = 0.5;
		v = (1. - RELAX) * v + RELAX * vnew;

		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_LWP_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (1)";
			exit(-1);
		}
	}

	double Tnew;
	iter = 0;
	err = 1.e5;
	T = 293;
	while (fabs(err) > 0.1) {
		cp = p1->gas->Get_cp(p, T);
		kappa_Tv = p1->gas->Get_kappa_Tv();
		Tnew = Tt - v * v / 2 / cp;
		if (inlet_pressure_drop)
			p = pt * pow(Tnew / Tt, kappa_Tv / (kappa_Tv - 1.));
		else
			p = pt;
		rho = p1->gas->Get_rho(p, Tnew);
		a = p1->gas->Get_SonicVel(Tnew, p);

		err = T - Tnew;
		T = Tnew;

		if (DEBUG)
			printf(
				"\n (front, inflow) outer iter #%2d: cp=%5.3e, kappa_Tv=%5.3e, p=%5.3e, T=%5.3e, v=%5.3e, err=%5.3e",
				iter, cp, kappa_Tv, p, T, v, err);

		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached. (2)";
			exit(-1);
		}

		iter++;
	}

	if (v < 0) {
		is_ok = false;
		if (DEBUG) {
			cout << endl << endl << " Connector_Reservoir_and_Pipe_Front() -> Assumed inflow, but computed v<0";
			//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err=%5.3e",
			//	iter, pt, rhot, beta, p, T, rho, v, err);
			cout << endl << "Trying outflow";
		}
	} else
		is_ok = true;

	// Need to try again with outflow
	if (!is_ok) {
		//double pM, vM;
		double TM = p1->Get_dprop("T_front");
		double pM = p1->Get_dprop("p_front");
		double rhoM = p1->gas->Get_rho(pM, TM);
		bool is_C0_ok = p1->GetC0AtFront(t_target);
		if (is_C0_ok) {
			rho = rhoM;
			T = TM * pow(rho / rhoM, kappa_Tv - 1.);
			p = p1->gas->Get_p(rho, T);
			a = p1->gas->Get_SonicVel(T, p);
			v = (p - beta) / rho / a;
			if (DEBUG) {
				printf("\n (front) OUTFLOW beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e", beta, p, T, rho,
				       v);
				if (v > 0) {
					cout << endl << "ERROR!!! Assumed outflow but computed inflow." << endl;
					//	exit(-1);
				}
			}
		} else {
			// Still inflow, set velocity to 0.
			T = TM;
			rho = rhoM;
			p = p1->gas->Get_p(rho, T);
			a = p1->gas->Get_SonicVel(T, p);
			v = (p - beta) / rho / a;
			if (DEBUG) {
				printf("\n (front) ERROR !! UNDECIDED beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e", beta,
				       p, T,
				       rho, v);
				//exit(-1);
			}
		}
	}
}

double Connector::Connector_LWP_Reservoir_and_Pipe_Front_fun(double v, double beta, Reservoir *r1, LWP *p1,
                                                             bool inlet_pressure_drop) {
	double pt = r1->Get_dprop("p");
	double Tt = r1->Get_dprop("T");
	//double rhot = p1->gas->Get_rho(pt, Tt);
	double cp = p1->gas->Get_cp(pt, Tt);
	double kappa_Tv = p1->gas->Get_kappa_Tv();

	// The absolute value is very important!
	double T = Tt - v * fabs(v) / 2. / cp;
	double p = pt;
	if (inlet_pressure_drop)
		p = pt * pow(T / Tt, kappa_Tv / (kappa_Tv - 1.));
	double rho = p1->gas->Get_rho(p, T);
	double a = p1->gas->Get_SonicVel(T, p);
	return (p - rho * a * v) - beta;
}

double Connector::signed_sqrt(double x) {
	if (x > 0)
		return sqrt(x);
	else
		return (-sqrt(-x));
}

bool Connector::Connector_SCP_Pipe_Back_and_Valve(double t_target,
                                                  SCP *p1, Valve *v1, double p_downstream, double rho, double a,
                                                  double &p) {
	// Solve the following system for p,T,rho,v:
	// (1) alpha=p+ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow

//cout<<endl<<"p_downstream:"<<p_downstream;
//cout<<endl<<"rho         :"<<rho;
//cout<<endl<<"a           :"<<a;
//cin.get();

	if (DEBUG)
		cout << endl << endl << "Connector::Connector_SCP_Pipe_Back_and_Valve: entering DEBUG mode" << endl;

	bool update_OK = true;

	double f, f1, df, dp, v, Apipe = p1->Get_dprop("A");
	double x = v1->Get_dprop("x");
	double alpha = p1->GetAlphaPrimitiveAtEnd(t_target);

	double err1, err2, err = 1.e5, mp, TOL = 10., pnew;
	int iter = 0, MAX_ITER = 100;
	double RELAX = 0.8;
	while ((fabs(err) > TOL) && (iter <= MAX_ITER)) {
		// 1. sol
		//mp  = v1->Get_MassFlowRate_InCompressible(p, p_downstream, rho, x);
		//v   = mp / rho / Apipe;
		//p   = alpha - rho * a * v;
		//err1 = (alpha - (p + rho * a * v)) / 1.e5;
		//err2 = rho * v * Apipe - (v1->Get_MassFlowRate_InCompressible(p, p_downstream, rho, x));

		// Newton (this is faaar quicker)
		mp = v1->Get_MassFlowRate(p, 293., p_downstream, 293., x);
		v = mp / rho / Apipe;
		f = p - (alpha - rho * a * v);

		dp = 0.01 * p;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp = v1->Get_MassFlowRate(p + dp, 293., p_downstream, 293., x);
		v = mp / rho / Apipe;
		f1 = (p + dp) - (alpha - rho * a * v);
		df = (f1 - f) / dp;

		pnew = p - f / df;
		p = RELAX * pnew + (1 - RELAX) * p;

		err1 = f;
		err2 = 0;
		if (DEBUG)
			printf("\n\t mp=%5.3e, x=%5.3e, v=%5.3e, f=%5.3e, f1=%5.3e, puj=%5.3e", mp, x, v, f, f1, p / 1e5);

		err = sqrt(err1 * err1 + err2 * err2);

		//printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
		//				iter, alpha, p, v, err1, err2);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_SCP_Pipe_Back_and_Valve() -> MAX_ITER reached, exiting";
			printf("\n (back) iter #%2d: alpha=%5.3e, p=%5.3e, v=%5.3e, err1=%+5.3e, err2=%+5.3e",
			       iter, alpha, p, v, err1, err2);
			update_OK = false;
		}
	}

	return update_OK;
}

bool Connector::Connector_SCP_Pipe_Front_and_Valve_Outlet(double t_target,
                                                          SCP *p1, Valve *v1, double p_upstream, double rho,
                                                          double a,
                                                          double &p) {
	// Solve the following system for p,T,rho,v:
	// (1) beta=p-ro*a*v
	// (2) ro*Apipe*v=valve_mass_flow

	if (DEBUG)
		cout << endl << endl << "Connector::Connector_SCP_Pipe_Back_and_Valve: entering DEBUG mode" << endl;

	bool update_OK = true;

	double f, f1, df, dp, v, Apipe = p1->Get_dprop("A");
	double x = v1->Get_dprop("x");
	double beta = p1->GetBetaPrimitiveAtFront(t_target);

	double err = 1.e5, mp, TOL = 10., pnew;
	int iter = 0, MAX_ITER = 100;
	double RELAX = 0.8;
	while ((fabs(err) > TOL) && (iter <= MAX_ITER)) {
		// 1. sol
		//mp  = v1->Get_MassFlowRate_InCompressible(p_upstream, p, rho, x);
		//v   = mp / rho / Apipe;
		//p   = beta + rho * a * v;
		//err1 = (beta - (p - rho * a * v)) / 1.e5;
		//err2 = rho * v * Apipe - (v1->Get_MassFlowRate_InCompressible(p_upstream, p, rho, x));

		// Newton (this is faaar quicker)
		mp = v1->Get_MassFlowRate(p_upstream, 293., p, 293., x);
		v = mp / rho / Apipe;
		f = p - (beta + rho * a * v);

		dp = 0.01 * p;
		if (fabs(dp) < 10.)
			dp = 10.;

		mp = v1->Get_MassFlowRate(p_upstream, 293., p + dp, 293., x);
		v = mp / rho / Apipe;
		f1 = (p + dp) - (beta + rho * a * v);
		df = (f1 - f) / dp;

		pnew = p - f / df;
		p = RELAX * pnew + (1 - RELAX) * p;

		err = fabs(f);
		if (DEBUG)
			printf("\n\t pu=%5.3e, pd=%5.3e, mp=%5.3e, x=%5.3e, v=%5.3e, err=%5.3e", p_upstream, p, mp, x, v,
			       err);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_SCP_Pipe_Front_and_Valve_Outlet() -> MAX_ITER reached, exiting";
			printf("\n (back) iter #%2d: beta=%5.3e, p=%5.3e, v=%5.3e, err=%+5.3e",
			       iter, beta, p, v, err);
			update_OK = false;
		}
	}

	return update_OK;
}


/*!
  Connects SCP front to reservoir
  \param t_target [in] target time value
  \param *r1 [in] pointer to reservoir
  \param *p1 [in] pointer to pipe
  \param rho [in] density, kg/m3
  \param a [in] dsonic velocity, m/s
  \param inlet pressure drop [in] bool, if true, inlet pressure drop issumed, i.e. pr=p+rho/2*v^2
  \param p [out] pressure at the first node of the pipe, result of computations
 */

void Connector::Connector_SCP_Reservoir_and_Pipe_Front(double t_target,
                                                       Reservoir *r1, SCP *p1, double rho, double a,
                                                       bool inlet_pressure_drop, double &p) {
	// Solve the following system for p,T,rho,v:
	// (1) p-rho*a*v = beta_front   Charactersistic equation
	// (2a) p=pt-rho/2*v^2       if inlet_pressure_drop=true
	// (2b) p=pt                 if inlet_pressure_drop=false

	double IPD_mul = 0;
	if (inlet_pressure_drop)
		IPD_mul = 1.;
	double pt = r1->Get_dprop("p");
	double v = 0;
	double beta = p1->GetBetaPrimitiveAtFront(t_target);

	p = pt;
	double err1, err2, err = 1.e5, TOL = 1e-3;
	int iter = 0, MAX_ITER = 50;
	while (fabs(err) > TOL) {
		v = (p - beta) / rho / a; // "Corrector"
		p = pt - IPD_mul * rho / 2 * v * v;
		err1 = pt - (p + IPD_mul * rho * v * v / 2.);
		err2 = (p - rho * a * v - beta) / 1.e5;
		err = sqrt(err1 * err1 + err2 * err2);

		//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);

		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector_Reservoir_and_Pipe_Front() -> MAX_ITER reached.";
			printf("\n (front) iter #%2d: pt=%5.3e, p=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e",
			       iter, pt, p, v, err1, err2);
			cin.get();
		}
	}
	//cin.get();
}

void Connector::Connector_SCP_Pipe_Simple_BC(double t_target) {
	if (is_front1) {
		e1->Set_BC_Left(BC_type, BC_value);
		//cout<<endl<<"Simple_BC - updated front of pipe "<<e1->Get_name();
	} else {
		e1->Set_BC_Right(BC_type, BC_value);
		//cout<<endl<<"Simple_BC - updated back of pipe "<<e1->Get_name();
	}
}


void Connector::Connector_SCP_2Pipes(double t_target, int update_idx) {
	// Solves:
	// p + d1*rho*a*v1 = alpha1
	// p + d2*rho*a*v2 = alpha2
	// p + d3*rho*a*v3 = alpha3
	// d1*A1*v1 + d2*A2*v2 + d3*A3*v3 = mpout/rho
	//
	// Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
	//
	// is_front = false -> d1=+1, alpha=pL+rho*a*vL
	// is_front = true  -> d1=-1, alpha=pR-rho*a*vR

	double alpha1, alpha2;
	// double A1  = e1->Get_dprop("A");
	// double A2  = e2->Get_dprop("A");
	// double rho1 = e1->Get_dprop("rho");
	// double rho2 = e2->Get_dprop("rho");
	// double a1   = e1->Get_dprop("a");
	// double a2   = e2->Get_dprop("a");
	// int d1 = 1, d2 = 1;
	double coeff_Q1, coeff_Q2;
	if (is_front1) {
		e1->GetBetaAtFront(t_target, alpha1, coeff_Q1);
		//alpha1 = e1->GetBetaAtFront(t_target);
		//d1 = -1;
	} else
		e1->GetAlphaAtEnd(t_target, alpha1, coeff_Q1);
	//alpha1 = e1->GetAlphaAtEnd(t_target);

	if (is_front2) {
		e2->GetBetaAtFront(t_target, alpha2, coeff_Q2);
		//alpha2 = e2->GetBetaAtFront(t_target);
		//d2 = -1;
	} else
		e2->GetAlphaAtEnd(t_target, alpha2, coeff_Q2);
	//alpha2 = e2->GetAlphaAtEnd(t_target);

	double p; //, v1,v2,v3, rhoa1, rhoa2,rhoa3;
	//rhoa1=rho1*a1;
	//rhoa2=rho2*a2;
	//rhoa3=rho3*a3;
	//p  = (alpha1*A1/rhoa1 + alpha2*A2/rhoa2 + alpha3*A3/rhoa3 - demand) / (A1/rhoa1 + A2/rhoa2 + A3/rhoa3);
	p = (alpha1 / coeff_Q1 + alpha2 / coeff_Q2 - demand) / (1. / coeff_Q1 + 1. / coeff_Q2);

	if (update_idx == edges_idx.at(0)) {
		if (is_front1)
			e1->Set_BC_Left("Pressure", p);
		else
			e1->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(1)) {
		if (is_front2)
			e2->Set_BC_Left("Pressure", p);
		else
			e2->Set_BC_Right("Pressure", p);
	}
}

void Connector::Connector_SCP_3Pipes(double t_target, int update_idx) {
	// Solves:
	// p + d1*rho*a*v1 = alpha1
	// p + d2*rho*a*v2 = alpha2
	// p + d3*rho*a*v3 = alpha3
	// d1*A1*v1 + d2*A2*v2 + d3*A3*v3 = mpout/rho
	//
	// Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
	//
	// is_front = false -> d1=+1, alpha=pL+rho*a*vL
	// is_front = true  -> d1=-1, alpha=pR-rho*a*vR

	double alpha1, alpha2, alpha3;
	double coeff_Q1, coeff_Q2, coeff_Q3;
	if (is_front1)
		e1->GetBetaAtFront(t_target, alpha1, coeff_Q1);
	else
		e1->GetAlphaAtEnd(t_target, alpha1, coeff_Q1);

	if (is_front2)
		e2->GetBetaAtFront(t_target, alpha2, coeff_Q2);
	else
		e2->GetAlphaAtEnd(t_target, alpha2, coeff_Q2);

	if (is_front3)
		e3->GetBetaAtFront(t_target, alpha3, coeff_Q3);
	else
		e3->GetAlphaAtEnd(t_target, alpha3, coeff_Q3);

	double p = (alpha1 / coeff_Q1 + alpha2 / coeff_Q2 + alpha3 / coeff_Q3 - demand) / (
		           1. / coeff_Q1 + 1. / coeff_Q2 + 1. / coeff_Q3);

	if (update_idx == edges_idx.at(0)) {
		if (is_front1)
			e1->Set_BC_Left("Pressure", p);
		else
			e1->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(1)) {
		if (is_front2)
			e2->Set_BC_Left("Pressure", p);
		else
			e2->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(2)) {
		if (is_front3)
			e3->Set_BC_Left("Pressure", p);
		else
			e3->Set_BC_Right("Pressure", p);
	}
}

void Connector::Connector_SCP_4Pipes(double t_target, int update_idx) {
	double alpha1, alpha2, alpha3, alpha4;
	double coeff_Q1, coeff_Q2, coeff_Q3, coeff_Q4;

	if (is_front1)
		e1->GetBetaAtFront(t_target, alpha1, coeff_Q1);
	else
		e1->GetAlphaAtEnd(t_target, alpha1, coeff_Q1);

	if (is_front2)
		e2->GetBetaAtFront(t_target, alpha2, coeff_Q2);
	else
		e2->GetAlphaAtEnd(t_target, alpha2, coeff_Q2);

	if (is_front3)
		e3->GetBetaAtFront(t_target, alpha3, coeff_Q3);
	else
		e3->GetAlphaAtEnd(t_target, alpha3, coeff_Q3);

	if (is_front4)
		e4->GetBetaAtFront(t_target, alpha4, coeff_Q4);
	else
		e4->GetAlphaAtEnd(t_target, alpha4, coeff_Q4);

	double p;
	p = (alpha1 / coeff_Q1 + alpha2 / coeff_Q2 + alpha3 / coeff_Q3 + alpha4 / coeff_Q4 - demand) / (
		    1. / coeff_Q1 + 1. / coeff_Q2 + 1. / coeff_Q3 + 1. / coeff_Q4);

	if (update_idx == edges_idx.at(0)) {
		if (is_front1)
			e1->Set_BC_Left("Pressure", p);
		else
			e1->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(1)) {
		if (is_front2)
			e2->Set_BC_Left("Pressure", p);
		else
			e2->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(2)) {
		if (is_front3)
			e3->Set_BC_Left("Pressure", p);
		else
			e3->Set_BC_Right("Pressure", p);
	}

	if (update_idx == edges_idx.at(3)) {
		if (is_front4)
			e4->Set_BC_Left("Pressure", p);
		else
			e4->Set_BC_Right("Pressure", p);
	}
}

void Connector::setPressureDependentDemand(double k, double p0, double t0) {
	pressure_dependent_demand = true;
	demand_k = k;
	demand_p0 = p0;
	demand_t0 = t0;
}


void Connector::Connector_SCP_NPipes(double t_target, int update_idx) {
	if (DEBUG) {
		cout << endl << "Connector_SCP_NPipes() - name of connector : " << name;
		cout << endl << "\t t_target = " << t_target;
	}

	vector<double> alpha(edges.size());
	vector<double> coeff_Q(edges.size());
	vector<double> alpha_coeff_Q(edges.size());
	vector<double> delta(edges.size());

	double p = p_prev, p_new;
	double b0 = -demand, b1 = 0., b2 = 0.;
	double err_max = 0.0001 * 1000 * 9.81, err = 10 * err_max;
	int iter = 0, iter_max = 1000000;
	int idx;


	while ((err > err_max) || (iter < 3)) {
		iter++;
		if (iter == iter_max) {
			cout << endl << "!!!!!!!!!!! Connector::Connector_SCP_NPipes" << endl;
			cout << "iter == " << iter << " err = " << err << " Node = " << name << endl;
			while (1) { 1; };
		}


		for (unsigned int i = 0; i < edges_idx.size(); i++) {
			idx = edges_idx.at(i);

			if (is_front.at(i)) {
				edges.at(idx)->GetBetaAtFront(t_target, alpha.at(i), coeff_Q.at(i));
				delta.at(i) = -1.;
			} else {
				edges.at(idx)->GetAlphaAtEnd(t_target, alpha.at(i), coeff_Q.at(i));
				delta.at(i) = 1.;
			}
			if (DEBUG) {
				printf("\n\t\t edge %10s: %+5.3e =  p %+5.3e Q", edges.at(idx)->name.c_str(), alpha.at(i),
				       coeff_Q.at(i));
				printf(", type=%s", edges.at(idx)->edge_type.c_str());
			}
			alpha_coeff_Q.at(i) = 1.;
		}

		//cin.get();
		//b0 frissítése b0 = gyök(p / k)
		if (pressure_dependent_demand && t_target > demand_t0) {
			double delta_p = p - demand_p0;
			double lim = 100;
			if (delta_p > lim) demand = sqrt(delta_p / demand_k);
			else demand = sqrt(lim / demand_k);
			/*else if(delta_p > -lim) demand = delta_p * sqrt(lim / demand_k ) / lim;
			else demand = -sqrt(-delta_p / demand_k );*/

			//cout << "dp = " << delta_p << "\tdemand = " << demand << endl;
			if (isnan(delta_p)) {
				while (1) { 1; };
			}

			if (delta_p < 0) {
				cout << "Pressure dependent demand off!" << endl;
				pressure_dependent_demand = false;
				demand = 0;
			}
		}


		b0 = -demand;
		for (unsigned int i = 0; i < edges_idx.size(); i++)
			b0 *= coeff_Q.at(i);

		for (unsigned int i = 0; i < edges_idx.size(); i++)
			for (unsigned int j = 0; j < edges_idx.size(); j++)
				if (i != j)
					alpha_coeff_Q.at(i) *= coeff_Q.at(j);
		b1 = 0., b2 = 0.;
		for (unsigned int i = 0; i < edges_idx.size(); i++) {
			b1 += alpha_coeff_Q.at(i) * delta.at(i) * alpha.at(i);
			b2 += alpha_coeff_Q.at(i) * delta.at(i);
		}

		p_new = (b0 + b1) / b2;


		double RELAX = 0.0;
		p_new = RELAX * p + (1. - RELAX) * p_new;

		err = fabs(p - p_new);

		if (DEBUG) {
			cout << endl << "\t iter:" << iter << ", p=" << p << " Pa = " << p / 1000 / 9.81 << " mwc";
			cout << ", p_new=" << p_new << ", err=" << fabs(p - p_new);
			double QQ;
			for (unsigned int i = 0; i < edges_idx.size(); i++) {
				QQ = (alpha.at(i) - p) / coeff_Q.at(i);
				cout << endl << "\t\t Q(" << edges_idx.at(i) << ")=" << QQ << " m3/s = " << QQ * 3600 <<
						" m3/h";
			}
			cin.get();
		}

		p = p_new;


		for (unsigned int i = 0; i < edges_idx.size(); i++) {
			idx = edges_idx.at(i);
			if (edges.at(idx)->is_rigid_element) {
				//cout<<endl<<"updating rigid element: "<<edges.at(idx)->name<<", is_front: "<<is_front.at(i);
				if (is_front.at(i))
					edges.at(idx)->Set_BC_Left("Pressure", p);
				else
					edges.at(idx)->Set_BC_Right("Pressure", p);
			}
		}
		//cin.get();
	} //while = iteracio vege

	p_prev = p;

	if (DEBUG)
		cout << endl << "\t Iteration finished. Setting up boundary conditions...";

	for (unsigned int i = 0; i < edges_idx.size(); i++) {
		idx = edges_idx.at(i);
		if (update_idx == idx) {
			//cout<<endl<<"  "<<edges.at(update_idx)->name<<", is_front:"<<is_front.at(i)<<" setting pressure to "<<p;
			if (is_front.at(i))
				edges.at(idx)->Set_BC_Left("Pressure", p);
			else
				edges.at(idx)->Set_BC_Right("Pressure", p);
		}
	}

	if (DEBUG) {
		cout << endl << "===================================" << endl;
		cin.get();
	}
}

void Connector::Connector_LWP_Pipes(double t_target,
                                    LWP *p1, const bool is_front1,
                                    LWP *p2, const bool is_front2,
                                    LWP *p3, const bool is_front3) {
	// Solves:
	// p + d1*rho*a*v1 = alpha1
	// p + d2*rho*a*v2 = alpha2
	// p + d3*rho*a*v3 = alpha3
	// d1*A1*v1 + d2*A2*v2 + d3*A3*v3 = mpout/rho
	//
	// Solution: p=(sum Ai*alphai-mpout*a)/(sum Ai)
	//
	// is_front = false -> d1=+1, alpha=pL+rho*a*vL
	// is_front = true  -> d1=-1, alpha=pR-rho*a*vR

	//double mpout=0.; // Can be changed later to be a nodal demand
	double alpha1, alpha2, alpha3;
	const double A1 = p1->Get_dprop("A");
	const double A2 = p2->Get_dprop("A");
	const double A3 = p3->Get_dprop("A");
	double p = 0., p_new = 0., T1, T2, T3, a1, a2, a3, rho1, rho2, rho3;

	int d1 = 1, d2 = 1, d3 = 1;
	if (is_front1) {
		alpha1 = p1->GetBetaPrimitiveAtFront(t_target);
		T1 = p1->Get_dprop("T_front");
		a1 = p1->Get_dprop("a_front");
		rho1 = p1->Get_dprop("rho_front");
		d1 = -1;
	} else {
		alpha1 = p1->GetAlphaPrimitiveAtEnd(t_target);
		T1 = p1->Get_dprop("T_back");
		a1 = p1->Get_dprop("a_back");
		rho1 = p1->Get_dprop("rho_back");
	}

	if (is_front2) {
		alpha2 = p2->GetBetaPrimitiveAtFront(t_target);
		T2 = p2->Get_dprop("T_front");
		a2 = p2->Get_dprop("a_front");
		rho2 = p2->Get_dprop("rho_front");
		d2 = -1;
	} else {
		alpha2 = p2->GetAlphaPrimitiveAtEnd(t_target);
		T2 = p2->Get_dprop("T_back");
		a2 = p2->Get_dprop("a_back");
		rho2 = p2->Get_dprop("rho_back");
	}

	if (is_front3) {
		alpha3 = p3->GetBetaPrimitiveAtFront(t_target);
		T3 = p3->Get_dprop("T_front");
		a3 = p3->Get_dprop("a_front");
		rho3 = p3->Get_dprop("rho_front");
		d3 = -1;
	} else {
		alpha3 = p3->GetAlphaPrimitiveAtEnd(t_target);
		T3 = p3->Get_dprop("T_back");
		a3 = p3->Get_dprop("a_back");
		rho3 = p3->Get_dprop("rho_back");
	}

	try {
		double err_p = 1.e5, ERR_P_MAX = 10.;
		double K1 = A1 * a2 * a3;
		double K2 = A2 * a1 * a3;
		double K3 = A3 * a1 * a2;
		int iter = 0, ITER_MAX = 20;
		while (err_p > ERR_P_MAX) {
			//p_new  = (A1 * alpha1 + A2 * alpha2 + A3 * alpha3 - mpout * a) / (A1 + A2 + A3);
			p_new = (K1 * alpha1 + K2 * alpha2 + K3 * alpha3) / (K1 + K2 + K3);
			//a = p1->gas->Get_SonicVel(T,p);
			err_p = sqrt((p - p_new) * (p - p_new));
			p = p_new;
			if (DEBUG)
				printf("\n\t iter #%d/%d, p=%5.3e, err_p=%5.3e", iter, ITER_MAX, p_new, err_p);

			iter++;
			if (iter == ITER_MAX)
				throw std::runtime_error("Too many iterations: iter = " + std::to_string(iter));
		}
	} catch (int iter) {
		cout << endl << "!!! ERROR !!!";
		cout << endl << "Connector::LWP_Pipes() -> too many iterations!!!";
		cout << endl << "Exiting..." << endl << endl;
		exit(-1);
	}

	p = p_new;
	const double T = (T1 + T2 + T3) / 3.;
	const double v1 = (alpha1 - p) / (d1 * rho1 * a1);
	const double v2 = (alpha2 - p) / (d2 * rho2 * a2);
	const double v3 = (alpha3 - p) / (d3 * rho3 * a3);

	Set_LWP_BC(p1, t_target, is_front1, p, T, v1);
	Set_LWP_BC(p2, t_target, is_front2, p, T, v2);
	Set_LWP_BC(p3, t_target, is_front3, p, T, v3);
}

void Connector::Connector_LWP_ThrottleValve_LWP(double t_target, LWP *p1, ThrottleValve *tv, LWP *p2) {
	// Solves:
	// p + d1*rho1*a1*v1 = p + d1*a1*mp/A1 = alpha1
	// p + d2*rho2*a2*v2 = p + d2*a2*mp/A2 = alpha2
	//                h1 = h2
	//                mp = flow_rate(Cd,p1,T1,p2,T2)
	// T1 = T_upstream = T(end) @ pipe 1 if mp > 0
	// T2 = T_upstream = T(0)   @ pipe 2 if mp < 0

	const double A1 = p1->Get_dprop("A");
	const double A2 = p2->Get_dprop("A");

	double d1 = 1;
	double alpha1 = p1->GetAlphaPrimitiveAtEnd(t_target);
	double T1 = p1->Get_dprop("T_back");
	double a1 = p1->Get_dprop("a_back");
	double rho1 = p1->Get_dprop("rho_back");
	double pres1 = p1->Get_dprop("p_back");

	double d2 = -1;
	double alpha2 = p2->GetBetaPrimitiveAtFront(t_target);
	double T2 = p2->Get_dprop("T_front");
	double a2 = p2->Get_dprop("a_front");
	double rho2 = p2->Get_dprop("rho_front");
	double pres2 = p2->Get_dprop("p_front");
	double pres1_new = pres1, pres2_new = pres2;

	try {
		double err_p = 1.e5, ERR_P_MAX = 100., mp = 10. * rho1 * A1;
		int iter = 0, ITER_MAX = 100;
		while (err_p > ERR_P_MAX) {
			pres1 = alpha1 - d1 * a1 / A1 * mp;
			pres2 = alpha2 - d2 * a2 / A2 * mp;
			double f = mp - tv->Get_MassFlowRate(pres1, T1, pres2, T2);

			double dmp = mp * 0.01;
			if (fabs(dmp) < 0.001)
				dmp = 0.001;
			double d_pres1 = alpha1 - d1 * a1 / A1 * (mp + dmp);
			double d_pres2 = alpha2 - d2 * a2 / A2 * (mp + dmp);
			double d_f = (mp + dmp) - tv->Get_MassFlowRate(d_pres1, T1, d_pres2, T2);

			double df_dmp = (d_f - f) / dmp;
			double mp_new = mp - f / df_dmp;
			//double mp_new = tv->Get_MassFlowRate(pres1, T1, pres2, T2);

			double RELAX = 0.5;

			mp = (1 - RELAX) * mp + RELAX * mp_new;

			pres1_new = alpha1 - d1 * a1 / A1 * mp;
			pres2_new = alpha2 - d2 * a2 / A2 * mp;

			if (mp > 0)
				T2 = T1;
			else
				T1 = T2;

			err_p = sqrt(pow(pres1 - pres1_new, 2.) + pow(pres2 - pres2_new, 2.));

			pres1 = (1 - RELAX) * pres1 + RELAX * pres1_new;
			pres2 = (1 - RELAX) * pres2 + RELAX * pres2_new;


			a1 = p1->gas->Get_SonicVel(T1, pres1);
			a2 = p2->gas->Get_SonicVel(T2, pres2);

			if (DEBUG)
				printf(
					"\n\t iter #%d/%d, mp=%+5.3e, p1res=%5.3e, p2res=%5.3e, dp=%+5.3e, err_p=%5.3e, T1=%5.3e, T2=%5.3e",
					iter, ITER_MAX, mp, pres1, pres2, pres1_new - pres2_new, err_p, T1, T2);


			iter++;
			if (iter == ITER_MAX)
				throw std::runtime_error("Too many iterations: iter = " + std::to_string(iter));
		}
	} catch (int iter) {
		cout << endl << "!!! ERROR !!!";
		cout << endl << "Connector::LWP_Pipes() -> too many iterations!!!";
		cout << endl << "Exiting..." << endl << endl;
		exit(-1);
	}
	const double v1 = (alpha1 - pres1_new) / (d1 * rho1 * a1);
	const double v2 = (alpha2 - pres2_new) / (d2 * rho2 * a2);

	Set_LWP_BC(p1, t_target, false, pres1_new, T1, v1);
	Set_LWP_BC(p2, t_target, true, pres2_new, T2, v2);
}

void Connector::Set_LWP_BC(LWP *p1, const double t_target, const bool is_front,
                           const double p, const double T, const double v) {
	double is_inflow, is_outflow;
	//cout << endl << "Set_LWP_BC: p = " << p / 1e5 << " bar, v = " << v << " m/s, T=" << T << " K, is_front = " <<
	//		is_front;

	if (is_front) {
		if (v > 0)
			is_inflow = p1->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", p, T, true);
		else
			is_outflow = p1->BCLeft(t_target, "StaticPres_Outlet", p, T, true);
	} else {
		if (v < 0)
			is_inflow = p1->BCRight(t_target, "StaticPres_and_StaticTemp_Inlet", p, T, true);
		else
			is_outflow = p1->BCRight(t_target, "StaticPres_Outlet", p, T, true);
	}
}

bool Connector::connected_to_rigid() {
	for (int i: edges_idx) {
		if (edges[i]->is_rigid_element) {
			return true;
		}
	}
	return false;
}

void Connector::Connector_LWP_Pipes(double t_target,
                                    LWP *p1, LWP *p2) {
	/*
       void Connector::Connector_LWP_Pipes(double t_target,
       LWP *p1, LWP *p2,
       double& pL, double& pR,
       double& TL, double& TR,
       double& rhoL, double& rhoR,
       double& vL, double& vR) {
	 */
	//bool DEBUG=true;

	// Solves:
	// eq. (1) pL + rhoL*aL*vL = alpha_L
	// eq. (2) pR - rhoR*aR*vR = beta_R
	// eq. (3) rhoL*vL*AL   = rhoR*vR*AR (continuity)
	// eq. (4) rhoL*vL^2*AL = rhoR*vR^2*AR+(pR-pL)*max(AL,AR) (impulse)
	// eq. (5) (rhoL*vL*eL + pL*vL)*AL = (rhoR*vR*eR + pR*vR)*AR (energy)
	// eq. (6-7) a = sonic_vel(T,p)
	// eq. (8-9) e = internal_energy(T,p)
	// eq. (10-11) gas law
	// eq. (12) rhoL = rho0 (C0 characteristic)
	//
	// 12 unkowns: p, T, rho, v, a, e @ R(ight) and L(eft)

	bool is_inflow, is_outflow;

	const double AR = p1->Get_dprop("A");
	const double AL = p2->Get_dprop("A");
	double A = AR;
	if (AL > AR)
		A = AL;

	double pL = p1->Get_dprop("p_back");
	double TL = 293.; //p1->Get_dprop("T_back");
	double vL = p1->Get_dprop("v_back");
	double rhoL = p1->gas->Get_rho(pL, TL);
	double aL = p1->gas->Get_SonicVel(TL, pL);
	double eL = p1->gas->Get_e_from_Tp(TL, pL) + vL * vL / 2;

	double pR = p2->Get_dprop("p_front");
	double TR = 293.; //p2->Get_dprop("T_front");
	double vR = p2->Get_dprop("v_front");
	double rhoR = p2->gas->Get_rho(pR, TR);
	double aR = p2->gas->Get_SonicVel(TR, pR);
	double eR = p1->gas->Get_e_from_Tp(TR, pR) + vR * vR / 2;

	double pR_new, pL_new; //,rhoR_new,rhoL_new,vL_new,vR_new;

	double alphaL = p1->GetAlphaPrimitiveAtEnd(t_target);
	double betaR = p2->GetBetaPrimitiveAtFront(t_target);
	if (DEBUG) {
		printf("\n\n Connector::Connector_LWP_Pipes()");
		printf("\n pL=%5.3e, TL=%5.3e, vL=%5.3e, rhoL=%5.3e, aL=%5.3e, eL=%5.3e, alphaL=%5.3e",
		       pL, TL, vL, rhoL, aL, eL, alphaL);
		printf("\n pR=%5.3e, TR=%5.3e, vR=%5.3e, rhoR=%5.3e, aR=%5.3e, eR=%5.3e, betaR =%5.3e",
		       pR, TR, vR, rhoR, aR, eR, betaR);
		cin.get();
	}

	// Preliminary esrimate assuming constant temperature
	double aa = 1. / rhoR / AR - 1. / rhoL / AL;
	double bb = aR + aL;
	double cc = betaR * AR - alphaL * AL;
	double DD = bb * bb - 4 * aa * cc;
	double mp = 0.0;
	if (DD > 0) {
		mp = (-bb + sqrt(DD)) / 2. / aa;
	} else {
		cout << endl << endl << "ERRORRRRRR!" << endl;
		cin.get();
	}

	pL = alphaL - aL / AL * mp;
	pR = betaR + aR / AR * mp;
	rhoL = p1->gas->Get_rho(pL, TL);
	//rhoR = p1->GetC0AtEnd(t_target);
	rhoR = p1->gas->Get_rho(pR, TR);
	vL = mp / rhoL / AL;
	vR = mp / rhoR / AR;

	if (DEBUG) {
		printf("\n\nEstimate:");
		printf("\n mp=%5.3e", mp);
		printf("\n pL=%5.3e, TL=%5.3e, vL=%5.3e, rhoL=%5.3e", pL, TL, vL, rhoL);
		printf("\n pR=%5.3e, TR=%5.3e, vR=%5.3e, rhoR=%5.3e", pR, TR, vR, rhoR);
		cin.get();
	}

	double err_p = 1e5, TOL_p = 10.;
	double err_T = 1e5, TOL_T = .1;
	int iter = 0, MAX_ITER = 1000;

	while ((fabs(err_p) > TOL_p) || (fabs(err_T) > TOL_T)) {
		iter++;
		if (iter == MAX_ITER) {
			cout << endl << "!!!ERROR!!! Connector::Connector_LWP_Pipes() -> MAX_ITER reached. Exiting" << endl;
			//printf("\n (front) iter #%2d: pt=%5.3e, rhot=%5.3e, beta=%5.3e, p=%5.3e, T=%5.3e, rho=%5.3e, v=%5.3e, err1=%5.3e, err2=%5.3e, err3=%5.3e", iter,pt,rhot,beta,p,T,rho,v,err1,err2,err3);
			exit(-1);
		}
		// eq. (1) & (2)
		double tmpL = rhoL * vL * vL * AL;
		double tmpR = rhoR * vR * vR * AR;
		double Atmp = max(AL, AR);
		//if (iter<3)
		//pL_new = alphaL - rhoL * aL * vL;
		pL_new = tmpR / Atmp + pR - tmpL / Atmp;

		pR_new = betaR + rhoR * aR * vR;


		double v_TRESHOLD = 0.01;
		vL = (alphaL - pL_new) / rhoL / aL;
		// eq. (12)
		if (fabs(vL) > v_TRESHOLD) {
			rhoL = p1->GetC0AtEnd(t_target);
			mp = rhoL * vL * AL;
			if (rhoL < 0) {
				cout << endl << "!!! ERROR !!!! rhoL=" << rhoL << " < 0 !!!" << endl;
				cin.get();
			}
			vR = mp / AR / rhoR;
			rhoR = (pL + rhoL * eL - pR) / eR;
		} else {
			pL_new = alphaL - rhoL * aL * vL;
			vL = 0.;
			vR = 0.;
			rhoR = p2->GetC0AtFront(t_target);
			rhoL = rhoR;
		}

		// eq. (6-7)
		aL = p1->gas->Get_SonicVel(TL, pL_new);
		aR = p2->gas->Get_SonicVel(TR, pR_new);
		// eq. (8-9)
		eL = p1->gas->Get_e_from_Tp(TL, pL_new) + vL * vL / 2.;
		eR = p2->gas->Get_e_from_Tp(TR, pR_new) + vR * vR / 2.;
		// eq. (10-11)
		double TL_new = p1->gas->Get_T(rhoL, pL_new);
		double TR_new = p2->gas->Get_T(rhoR, pR_new);

		err_T = sqrt(pow(TL - TL_new, 2.) + pow(TR - TR_new, 2.));
		err_p = sqrt(pow(pL - pL_new, 2.) + pow(pR - pR_new, 2.));
		pL = (pL + pL_new) / 2.;
		pR = (pR + pR_new) / 2.;
		TL = (TL + TL_new) / 2.;
		TR = (TR + TR_new) / 2.;
		mp = vL * rhoL * AL;
		double mpR = vR * rhoR * AR;
		if (DEBUG) {
			printf(
				"\n\t iter=%2d, (L,R) p=%5.3f,%5.3f, T=%5.3f,%5.3f, rho=%5.3f,%5.3f, v=%5.3f,%5.3f, mp=%5.3f, %5.3f, err=%5.3e,%5.3e",
				iter, pL_new / 1.e5, pR_new / 1.e5, TL, TR, rhoL, rhoR, vL, vR, mp, mpR, err_T, err_p);
		}
	}
	if (vR > 0) {
		is_outflow = p1->BCRight(t_target, "StaticPres_Outlet", pL, 0, true);
		is_inflow = p2->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", pR, TR, true);
	} else {
		is_outflow = p2->BCRight(t_target, "StaticPres_Outlet", pR, 0, true);
		is_inflow = p1->BCLeft(t_target, "StaticPres_and_StaticTemp_Inlet", pL, TL, true);
	}
}
