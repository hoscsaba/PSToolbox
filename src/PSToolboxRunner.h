#pragma once

#include <vector>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include <Eigen/Dense>


struct RuntimeModifier {
	double time;
	double set_time = 0.0;
	bool done = false;
	bool has_started = false;
	vector<int> connectors;
	vector<double> demands;
	vector<double> old_demands;

	void Add(int idx, double demand) {
		connectors.push_back(idx);
		demands.push_back(demand);
	};
};


class PSToolboxRunner {
public:
	PSToolboxRunner(
		vector<PSToolboxBaseEdge *> &,
		vector<Connector *> &,
		vector<int> &,
		vector<int> &,
		vector<vector<int> > &);

	void Run(double);

	void Set_Save_data(bool newval) { save_data = newval; };

	void Save_data_for_ini();

	void Ini_from_file(string fname, double dt_target);

	void Set_Folder(string newval) { folder = newval; };
	void Set_DEBUG(bool newval) { DEBUG = newval; };
	void Set_Node_mul(int newval) { Node_mul = newval; };

	void Set_Pressure_limit(bool set, double _p_limMin, double _p_limMax, string _fname, string _type) {
		set_limit = set;
		p_limMin = _p_limMin;
		p_limMax = _p_limMax;
		limits_fname = _fname;
		limits_type = _type;
	};

	void Set_Pressure_limit_time(double start, double end) {
		p_start_time = start;
		p_end_time = end;
	};
	void Set_Save_interval(double interval) { save_interval = interval; };
	void Set_Write_interval(double interval) { write_interval = interval; };
	void Add_Runtime_modifier(RuntimeModifier &rmod) { modifiers.push_back(rmod); };

	void ListEdgeInfo() {
		for (int i = 0; i < edges.size(); i++) {
			cout << edges.at(i)->Info();
			cin.get();
		}
	}

	void Ini(double dt_target);

private:
	vector<PSToolboxBaseEdge *> edges;
	vector<Connector *> cons;
	vector<int> con_at_edge_start;
	vector<int> con_at_edge_end;
	bool DEBUG;
	int Node_mul;
	bool ini_done;

	//min-max pressure
	double p_limMin;
	double p_limMax;
	double p_start_time;
	double p_end_time;

	bool set_limit;
	string limits_fname;
	string limits_type;

	//save stuff
	bool save_data;
	double write_interval;
	double save_interval;
	vector<vector<int> > rigid_subsystems;
	vector<vector<int> > flexible_pipes_connected_to_subsystems;
	vector<vector<bool> > flexible_pipes_connected_to_subsystems_is_front;

	void Update_RS(int idx_of_rs, double t_target, bool save);

	void RS_Solve(int rs_idx, double t_target, bool rs_save);

	vector<int> num_of_edges;
	// list of nodes in the rs
	vector<vector<int> > rs_idx_nodes;
	// list of edges in the rs
	vector<vector<int> > rs_idx_edges;
	//  rs_idx_nodes_edges.at(i).at(j) : edges connected to node j in rs i
	vector<vector<vector<int> > > rs_idx_nodes_edges;
	//  rs_idx_nodes_edges.at(i).at(j) : edge direction connected to node j in rs i
	vector<vector<vector<bool> > > rs_idx_nodes_is_front;

	//runtime modifiers
	vector<RuntimeModifier> modifiers;

	//save folder
	string folder;
};
