#pragma once

#include <vector>
#include "PSToolboxBaseEdge.h"
#include "BioConnector.h"
#include <Eigen/Dense>


struct RuntimeBioModifier {
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

class PSToolboxBioRunner {
public:
	PSToolboxBioRunner(
		vector<PSToolboxBaseEdge *> &,
		vector<BioConnector *> &,
		vector<int> &,
		vector<int> &,
		vector<vector<int> > &);

	vector<PSToolboxBaseEdge *> edges;
	vector<BioConnector *> cons;

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
	//void Set_Write_interval(double interval) {write_interval = interval;};
	void Add_Runtime_modifier(RuntimeBioModifier &rmod) { modifiers.push_back(rmod); };

	void ListEdgeInfo() {
		for (int i = 0; i < edges.size(); i++) {
			cout << edges.at(i)->Info();
			cin.get();
		}
	}

	void Ini(double dt_target);

	void Ini(double dt_target, VectorXd Cini);

	void Load_v_conv_from_file(string fname, double mul);

	void Load_v_conv_from_epanet_export(string fname, double mul);

	void Load_nodal_demand_from_epanet_export(string fname, double mul);

	void Set_snapshot_fname(string newval) { snapshot_fname = std::move(newval); }
	void Set_save_snapshot(bool newval) { save_snapshot = newval; }

	void build_statistics();

private:
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

	int snapshot_num;
	string snapshot_fname;
	bool save_snapshot;

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
	vector<RuntimeBioModifier> modifiers;

	//save folder
	string folder;

	string convertSecondsToTimeFormat(double totalSeconds);

	double limitAbs(double value, double maxAbs);

	void Write_snapshot(double);
};
