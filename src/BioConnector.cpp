#define USE_MATH_DEFINES

#include <cstdlib>
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
#include "BioConnector.h"
#include "BioPipe.h"
#include "PSToolboxBaseEdge.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>

BioConnector::~BioConnector() = default;

BioConnector::BioConnector(
	string _name,
	vector<PSToolboxBaseEdge *> &_e,
	vector<bool> &_is_front,
	vector<int> &_edges_idx,
	double _demand, bool _DEBUG) {
	name = std::move(_name);
	edges = _e;
	is_front = _is_front;
	edges_idx = _edges_idx;
	demand = _demand;
	DEBUG = _DEBUG;
	type = 0;

	BC_type = "none";
	BC_value = 0;

	// DEBUG=true;
}

BioConnector::BioConnector(
	string _name,
	PSToolboxBaseEdge *_e1, bool _is_front1,
	string _BC_type, double _BC_value,
	double _demand, bool _DEBUG) {
	name = std::move(_name);
	s_edges = _e1;
	s_is_front = _is_front1;
	demand = _demand;
	DEBUG = _DEBUG;

	BC_type = _BC_type;

	if (BC_type == "Velocity")
		type = 1;
	else {
		if (BC_type == "Pressure")
			type = 2;
		else {
			cout << endl << "ERROR!!! BioConnector::BioConnector() BC_type=" << BC_type << " is not recognized" << endl;
			exit(-1);
		}
	}
	BC_value = _BC_value;
	//DEBUG=true;
}

//void BioConnector::SetInletQuality(VectorXd _Cin){
//	Cin=_Cin;
//  }

// This is the main update function
void BioConnector::Update(double t_target, int update_idx) {
	switch (type) {
		case 0:
			BioConnector_N_BioPipes(t_target, update_idx);
			break;

		case 1:
			BioConnector_FreeEnd(t_target);
			break;

		case 2:
			BioConnector_Reservoir(t_target);
			break;

		default:
			cout << endl << "ERROR!!!";
			cout << endl << "BioConnector::Update() -> unknown type=" << type << endl;
			cin.get();
			break;
	}
}


void BioConnector::BioConnector_N_BioPipes(double t_target, int update_idx) {
	//if (name == "52466")
	//	DEBUG = true;
	//else
	DEBUG = false;


	if (DEBUG) {
		cout << endl << "BioConnector_N_BioPipes() - name of connector : " << name;
		cout << endl << "\t t_target = " << t_target;
	}

	vector<bool> is_inflow;
	for (unsigned int i = 0; i < edges_idx.size(); i++) {
		int idx = edges_idx.at(i);

		if ((is_front[i] && (edges[idx]->Get_dprop("v_conv") > 0)) || (
			    !is_front[i] && (edges[idx]->Get_dprop("v_conv") < 0)))
			is_inflow.push_back(false);
		else
			is_inflow.push_back(true);
	}

	if (DEBUG) {
		for (unsigned int i = 0; i < edges_idx.size(); i++) {
			int idx = edges_idx.at(i);
			cout << endl << "\t edge: " << edges[idx]->Get_name() << ", is_front: " << is_front[i] << ", v_conv: " <<
					edges[idx]->Get_dprop("v_conv") << ", is_inflow:" << is_inflow[i];
		}
		//cin.get();
	}

	double sum_mp = 0.0, mp;
	VectorXd C;
	VectorXd sum_C = VectorXd::Zero(edges.at(0)->Get_iprop("num_of_bio_vars"));

	int num_of_inflows = 0;
	for (unsigned int i = 0; i < edges_idx.size(); i++) {
		int idx = edges_idx.at(i);

		if (is_inflow.at(i)) {
			if (DEBUG)
				cout << endl << "Getting inflowing C for edge... " << edges.at(idx)->Get_name();
			if (is_front.at(i))
				C = edges.at(idx)->Get_C_Front();
			else
				C = edges.at(idx)->Get_C_Back();
			mp = edges.at(i)->Get_dprop("mp_back");
			sum_C += C * fabs(mp);
			sum_mp += fabs(mp);
			num_of_inflows++;

			if (DEBUG)
				cout << endl << "\tinflow, edge" << edges.at(idx)->Get_name() << " C=" << C.transpose() << ", sum_C=" <<
						sum_C.transpose() << ", sum_mp=" << sum_mp << endl;
			//cin.get();
		}
	}
	if (num_of_inflows > 0)
		C = sum_C / sum_mp;
	else {
		//cout << endl << "WARNING!!!! at connector " << name << " there is no inflow. Setting C to zero.";
		C = sum_C;
		//	cin.get();
	}

	for (unsigned int i = 0; i < edges_idx.size(); i++) {
		int idx = edges_idx.at(i);
		if (!is_inflow.at(i)) {
			if (is_front.at(i)) {
				edges.at(idx)->Set_inletQualityFront(C);
				if (DEBUG)
					cout << endl << "front of edge " << edges.at(idx)->Get_name() << " is set to " << C;
			} else {
				edges.at(idx)->Set_inletQualityBack(C);
				if (DEBUG)
					cout << endl << "back of edge " << edges.at(idx)->Get_name() << " is set to " << C;
			}
		}
	}

	if (DEBUG)
		cin.get();
}

void BioConnector::BioConnector_Reservoir(double t_target) {
	if (type == 0) {
		cout << endl << "ERROR!!!! BioConnector::BioConnector_Reservoir() was called with type 0." << endl;
		exit(-1);
	} else {
		//cout<<endl<<"Entering BioConnector::BioConnector_Reservoir()"<<endl;
		bool is_inflow = false;

		if ((s_is_front && (s_edges->Get_dprop("v_conv") > 0)) || (!s_is_front && (s_edges->Get_dprop("v_conv") < 0)))
			is_inflow = true;
		//cout<<endl<<"\t is_inflow = "<<is_inflow<<endl;
		if (is_inflow) {
			//cout<<endl<<"Inflow, C="<<Cin.transpose()<<endl;
			//cout<<endl<<"s_edges->name = "<<s_edges->name;
			//cout<<endl<<"s_is_front = "<<s_is_front<<endl;

			if (s_is_front)
				s_edges->Set_inletQualityFront(Cin);
			else
				s_edges->Set_inletQualityBack(Cin);
		}
	}
}

void BioConnector::BioConnector_FreeEnd(double t_target) {
	if (type == 0) {
		cout << endl << "ERROR!!!! BioConnector::BioConnector_FreeEnd() was called with type 0." << endl;
		exit(-1);
	}
	bool is_inflow = false;

	if ((s_is_front && (s_edges->Get_dprop("v_conv") > 0)) || (!s_is_front && (s_edges->Get_dprop("v_conv") < 0)))
		is_inflow = true;
	/*
		if (name == "48854") {
			cout << endl << "entering BioConnector_FreeEnd";
			cout << endl << "s_is_front=" << s_is_front;
			cout << endl << "v_conv    =" << s_edges->Get_dprop("v_conv");
			cout << endl << "is_inflow =" << is_inflow;
			cout << endl << "Cin       =" << Cin;
			cout << endl << "Npts      =" << s_edges->Get_iprop("Npts");
			cout << endl << "C         =" << s_edges->C;
			cin.get();
		}
	*/
	if (is_inflow) {
		if (s_is_front)
			s_edges->Set_inletQualityFront(Cin);
		else
			s_edges->Set_inletQualityBack(Cin);
	}

	/*
	if (name == "48854") {
		cout << endl << "leaving BioConnector_FreeEnd.";
		cin.get();
	}
	*/
}
