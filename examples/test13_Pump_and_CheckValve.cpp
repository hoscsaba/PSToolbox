/*
	-----------------------------------
	PSToolbox Example 13: System with pump and check valve

	Description: Simulates a system with two serially connected pumps and a bypass.
	In the simulation pump1 stops at 60s, then the simulation runs until 180s.

	Last updated: 2024. 05. 06.

	Status: Works without errors, steady state is reached

	-----------------------------------
*/

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include "SCP.cpp"
#include "Pump.cpp"
#include "Tank.cpp"
#include "CheckValve.cpp"
#include "Connector.cpp"
#include "PSToolboxBaseEdge.cpp"
#include "PSToolboxRunner.cpp"

using namespace std;

int main(int argc, char **argv) {
	bool DEBUG=false;
	vector<PSToolboxBaseEdge *> edges;
	vector<Connector *> cons;

	bool save_data=true;

	// define edges
	double rho=1000.;
	double Dp=0.1;
	double Q_nom=2*Dp*Dp*M_PI/4.;
	double H_nom=70., H_max=1.2*H_nom;
	vector<double> coeff_H, coeff_P;
	coeff_H.push_back(H_max);
	coeff_H.push_back(0.0);
	coeff_H.push_back((H_nom-H_max)/Q_nom/Q_nom);
	double zeta=H_nom/1000./Q_nom/Q_nom;

	double P_nom=9.81*rho*Q_nom*H_nom;
	double P_0=P_nom*0.2;
	coeff_P.push_back(P_0);
	coeff_P.push_back((P_nom-P_0)/Q_nom);

	/* 0 */ edges.push_back(new Tank("tank","nx","n0",rho,0.,save_data));
	/* 1 */ edges.push_back(new Pump("pump1","n0","n1",rho,Q_nom,H_nom,coeff_H,coeff_P,/*n_nom*/ 3000.,/*n_act*/ 3000.,0.1,save_data));
	/* 2 */ edges.push_back(new CheckValve("checkvalve_pump1","n1","n2",rho,zeta,10000.*zeta,save_data));
	/* 3 */ edges.push_back(new CheckValve("checkvalve_bypass","n0","n2",rho,zeta,10000.*zeta,save_data));
	/* 4 */ edges.push_back(new SCP("pipe1","n2","n3",rho,1300,1000,Dp,200,0,20,save_data));
	/* 5 */ edges.push_back(new Pump("pump2","n3","n4",rho,Q_nom,H_nom,coeff_H,coeff_P,/*n_nom*/ 3000.,/*n_act*/ 3000.,0.1,save_data));
	/* 6 */ edges.push_back(new CheckValve("checkvalve_pump2","n4","n5",rho,zeta,10000.*zeta,save_data));
	/* 7 */ edges.push_back(new SCP("pipe2","n5","n6",rho,1300,2000,Dp,100,20,50,save_data));
	/* 8 */ edges.push_back(new SCP("pipe3","n5","n7",rho,1300,3000,Dp,100,20,50,save_data));


	edges.at(4)->Set_string_prop("lambda_model","hw");
	edges.at(7)->Set_string_prop("lambda_model","hw");
	edges.at(8)->Set_string_prop("lambda_model","hw");
	//edges.at(1)->Add_transient("use_true_H_curve",1*60.,0.);
	edges.at(1)->Add_transient("switch_off",1*60.,0.);

	// define nodes (connectors)
	double demand_n0 = 0.0, demand_n1 = 0.0, demand_n2 = 0.0, demand_n3 = 0.0, demand_n4 = 0.0, demand_n5 = 0.0;

	vector<bool> is_front_n0; is_front_n0.push_back(false); is_front_n0.push_back(true); is_front_n0.push_back(true);
	vector<int> edges_idx_n0; edges_idx_n0.push_back(0);    edges_idx_n0.push_back(1);   edges_idx_n0.push_back(3);

	vector<bool> is_front_n1; is_front_n1.push_back(false); is_front_n1.push_back(true);
	vector<int> edges_idx_n1; edges_idx_n1.push_back(1);    edges_idx_n1.push_back(2);

	vector<bool> is_front_n2; is_front_n2.push_back(false); is_front_n2.push_back(false); is_front_n2.push_back(true);
	vector<int>  edges_idx_n2 ;edges_idx_n2.push_back(2);   edges_idx_n2.push_back(3); edges_idx_n2.push_back(4);

	vector<bool> is_front_n3; is_front_n3.push_back(false); is_front_n3.push_back(true);
	vector<int> edges_idx_n3; edges_idx_n3.push_back(4);    edges_idx_n3.push_back(5);

	vector<bool> is_front_n4; is_front_n4.push_back(false); is_front_n4.push_back(true);
	vector<int> edges_idx_n4; edges_idx_n4.push_back(5);    edges_idx_n4.push_back(6);

	vector<bool> is_front_n5; is_front_n5.push_back(false); is_front_n5.push_back(true); is_front_n5.push_back(true);
	vector<int> edges_idx_n5; edges_idx_n5.push_back(6);    edges_idx_n5.push_back(7); edges_idx_n5.push_back(8);

	/* 0 */ cons.push_back(new Connector("n0",edges,is_front_n0,edges_idx_n0,demand_n0,false));
	/* 1 */ cons.push_back(new Connector("n1",edges,is_front_n1,edges_idx_n1,demand_n1,false));
	/* 2 */ cons.push_back(new Connector("n2",edges,is_front_n2,edges_idx_n2,demand_n2,false));
	/* 3 */ cons.push_back(new Connector("n3",edges,is_front_n3,edges_idx_n3,demand_n3,false));
	/* 4 */ cons.push_back(new Connector("n4",edges,is_front_n4,edges_idx_n4,demand_n4,false));
	/* 5 */ cons.push_back(new Connector("n5",edges,is_front_n5,edges_idx_n5,demand_n5,false));
	/* 6 */ cons.push_back(new Connector("n6",edges.at(7),false,"Pressure",0.,0.0,false));
	/* 7 */ cons.push_back(new Connector("n7",edges.at(8),false,"NR",0.-rho*1300*2,0.0,false));

	// define rigid subsystems
	vector< vector<int> > rigid_subsystems;
	vector<int> rs1; rs1.push_back(0); rs1.push_back(1); rs1.push_back(2); rs1.push_back(3);
	vector<int> rs2; rs2.push_back(5); rs2.push_back(6);

	rigid_subsystems.push_back(rs1);
	rigid_subsystems.push_back(rs2);

	// We need to add here connectivity info, e.g. by adding edge pointer + start/end info if relevant
	// con_at_edge_start.at(i) stores the idx of the connector connected to the start of the edge
	// con_at_edge_end.at(i)   stores the idx of the connector connected to the end of the edge
	vector<int> con_at_edge_start(edges.size(),-1);
	vector<int> con_at_edge_end(edges.size(),-1);
	con_at_edge_start.at(0)=0; con_at_edge_end.at(0)=0;
	con_at_edge_start.at(1)=0; con_at_edge_end.at(1)=1;
	con_at_edge_start.at(2)=1; con_at_edge_end.at(2)=2;
	con_at_edge_start.at(3)=0; con_at_edge_end.at(3)=2;
	con_at_edge_start.at(4)=2; con_at_edge_end.at(4)=3;
	con_at_edge_start.at(5)=3; con_at_edge_end.at(5)=4;
	con_at_edge_start.at(6)=4; con_at_edge_end.at(6)=5;
	con_at_edge_start.at(7)=5; con_at_edge_end.at(7)=6;
	con_at_edge_start.at(8)=5; con_at_edge_end.at(8)=7;
	// =============== END OF DATA READER SECTION

	PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end, rigid_subsystems);
	r.Set_Save_data(true);
	r.Set_DEBUG(false);
	r.Set_Node_mul(5);
  	r.Set_Save_interval(1.0);
  	r.Set_Write_interval(3*60.0);
	string folder = "data";
  	r.Set_Folder(folder);
	r.Run(3*60.);

/*
	PSToolboxPlotter pl1("data/pipe1.dat","pipe"); pl1.Plot();
	PSToolboxPlotter pl2("data/pipe2.dat","pipe"); pl2.Plot();
	PSToolboxPlotter pl3("data/pipe3.dat","pipe"); pl3.Plot();

	PSToolboxPlotter pl4("data/pump1.dat","pump"); pl4.Plot();
	PSToolboxPlotter pl5("data/pump2.dat","pump"); pl5.Plot();

	PSToolboxPlotter pl6("data/checkvalve_pump1.dat","checkvalve"); pl6.Plot();
	PSToolboxPlotter pl7("data/checkvalve_bypass.dat","checkvalve"); pl7.Plot();*/
}
