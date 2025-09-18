using namespace std;

#include <string>
#include "PSToolboxBaseEdge.h"

PSToolboxBaseEdge::PSToolboxBaseEdge(const string _edge_type, const string _name, const string _n1, const string _n2) {
	name = _name;
	edge_type = _edge_type;
	node_from = _n1;
	node_to = _n2;
	dt = 0.;
	edge_type = _edge_type;
	save_data = false;
	DEBUG = false;
	Q = 0.;
	t = 0.;
	dt = 0.;
	g = 9.81;
	p1 = 0.;
	p2 = 0.;
	is_rigid_element = false;
	suppress_all_output = false;
	dt_out = 1.;
};
