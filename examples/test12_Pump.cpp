// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test12_Pump.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "/Users/hoscsaba/program/PSToolbox/SCP.h"
#include "/Users/hoscsaba/program/PSToolbox/Pump.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxBaseEdge.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxPlotter.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxRunner.h"

using namespace std;

int main(int argc, char **argv) {
	bool DEBUG=false;
	vector<PSToolboxBaseEdge *> edges;
	vector<Connector *> cons;

	bool save_data=true;

	//=============== THIS SECTION GOES TO THE DATA READER
	// define edges
	double rho=1000.;
	double Dp=0.1;
	double Q_nom=1*Dp*Dp*M_PI/4.;
	double H_nom=100., H_max=1.2*H_nom;
	vector<double> coeff_H, coeff_P;
	coeff_H.push_back(H_max);
	coeff_H.push_back(0.0);
	coeff_H.push_back((H_nom-H_max)/Q_nom/Q_nom);

	double P_nom=9.81*rho*Q_nom*H_nom;
	double P_0=P_nom*0.2;
	coeff_P.push_back(P_0);
	coeff_P.push_back((P_nom-P_0)/Q_nom);

	edges.push_back(new Pump("pump","n1","n2",rho,Q_nom,coeff_H,coeff_P,/*n_nom*/ 3000.,/*n_act*/ 3000.,0.1,save_data));
	edges.push_back(new SCP("p1","n2","n3",1000,1300,101,0.1,100.,0,0,save_data));

	// define nodes (connectors)
	double demand1 = 0.;
	double demand2 = 0.;
	double demand3 = 0.;
	vector<bool> is_front1;
	is_front1.push_back(false);
	is_front1.push_back(true);
	vector<int> edges_idx1;
	edges_idx1.push_back(0);
	edges_idx1.push_back(1);

	cons.push_back(new Connector("n1",edges.at(0),true,"Pressure",0.e5,demand1,DEBUG));
	cons.push_back(new Connector("n2",edges,is_front1,edges_idx1,demand2,DEBUG));
	cons.push_back(new Connector("n3",edges.at(1),false,"Pressure",0.e5,demand3,DEBUG));

	// We need to add here connectivity info, e.g. by adding edge pointer + start/end info if relevant
	// con_at_edge_start.at(i) stores the idx of the connector connected to the start of the edge
	// con_at_edge_end.at(i)   stores the idx of the connector connected to the end of the edge
	vector<int> con_at_edge_start(edges.size(),-1);
	vector<int> con_at_edge_end(edges.size(),-1);
	con_at_edge_start.at(0)=0; con_at_edge_end.at(0)=1;
	con_at_edge_start.at(1)=1; con_at_edge_end.at(1)=2;
	//con_at_edge_start.at(2)=1; con_at_edge_end.at(2)=3;
	// =============== END OF DATA READER SECTION

	PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);
	r.Set_Save_data(true);
	r.Set_DEBUG(false);
	r.Set_Node_mul(2);
	r.Run(10.);

//	PSToolboxPlotter pl1("p1.dat"); pl1.Plot();
	PSToolboxPlotter pl2("data/p1.dat"); pl2.Plot();
//	PSToolboxPlotter pl3("p3.dat"); pl3.Plot();

}
