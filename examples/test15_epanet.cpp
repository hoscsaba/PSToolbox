/*
 The model includes a simple double pipe, with a pump in the beginning and a reservoir (modeled with constant pressure) at the end.
*/

// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_Csaba.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "SCP.h"
#include "Connector.h"
#include "PSToolboxBaseEdge.h"
//#include "PSToolboxPlotter.h"
#include "PSToolboxRunner.h"
#include "CheckValve.h"
#include "Pump.h"
#include "EpanetReader.h"

using namespace std;

//set propogation velocity for different pipes
void calculatePropagationVelocity(EpanetReader & reader)
{
    //go through all the pipes
    for (int i = 0; i < reader.pipes.size(); i++)
    {
        double r = reader.pipes[i].Roughness;
        double D = reader.pipes[i].Diameter;
        double Ef = 2.18e9; //El. mod. of water
        double Ec, Delta = -1;

				//cases based on diameter are given here
        if(D == 600.0)
        {
            Ec = 170.0e9; 
            Delta = 8.7e-3;
        }
        else if(D == 547.4)
        {
            Ec = 900.0e6; 
            Delta = 41.3e-3;
        }
        else if(D == 503.8)
        {
            Ec = 900.0e6; 
            Delta = 63.1e-3;
        }
        else if(D == 800 || D == 850 || D == 2500)
        {
            Ec = 170.0e9; 
            Delta = 8.7e-3;
        }
        else
        {
            cout << "Pipe data not found: " <<  reader.pipes[i].ID <<endl;
        }

        double Er = 1.0/(1.0/Ef + 0.001*D/(Delta * Ec));

        reader.pipes[i].SpeedOfSound = sqrt(Er/1000.0);
        reader.pipes[i].Delta = Delta;
        
        if(Delta == -1)
        {
            cout << "Pipe " << i << " ID: " << reader.pipes[i].ID << endl;
            cout << "\tRoughness= " << reader.pipes[i].Roughness << endl;
            cout << "\t       D = " << reader.pipes[i].Diameter << endl;
            cout << "\t   Delta = " << reader.pipes[i].Delta << endl;
            cout << "\t       a = " << reader.pipes[i].SpeedOfSound << endl << endl;
            while(true) 1;
        }
        cout << reader.pipes[i].ID << "," << reader.pipes[i].Roughness  << "," <<  reader.pipes[i].Diameter << "," <<  reader.pipes[i].Delta << "," <<  reader.pipes[i].SpeedOfSound << endl;
    }
    
}

int main(int argc, char **argv) 
{
  /*
  ================= Simulation settings =================
  */
  string folder = "epanet_test";

  /*
  ================= Manege libraries =================
  */
  string path = folder + "/*";
  string removecom = "rm " + path;
  int o = system(removecom.c_str());
  string mkdir = "mkdir " + folder;
  o = system(mkdir.c_str());

  /*
  ================ read EPANET file ================
  */
  EpanetReader reader;
  string location = "epanet_test.inp";
  reader.readFromFile(location);
  calculatePropagationVelocity(reader);
  reader.convertToRunner2();
  
  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;

  /*
  ================ Add pump to nodes ================
  */
  int idx_sziv1 = reader.findPipeByID("sziv1");
  int idx_sziv2 = reader.findPipeByID("sziv2");
  int idx_sziv3 = reader.findPipeByID("sziv3");
  int idx_szel1 = reader.findPipeByID("szel1");
  int idx_szel2 = reader.findPipeByID("szel2");
  int idx_szel3 = reader.findPipeByID("szel3");

  double dp = 0.6e5;
  double Q_szelep = 500 / 3600.0; //500CMH mellett legyen 0.2bar pressure loss
  double k_open = dp / Q_szelep / Q_szelep;
  double zeta_open = k_open / 1000 / 9.81;
  double zeta_closed = 1.0e6*zeta_open;
  cout << "zeta open = " << zeta_open << " zeta close = " << zeta_closed << endl;

	//parameters
  double n_nom = 2880;
  double n_act = 48*60;
  double Q_nom = 460 / 3600;
  double H_nom = 60.0;
  double theta = 3.0; 

	//charactestic curves
	vector<double> coeff_H, coeff_P;
	coeff_H.push_back(70.0);
	coeff_H.push_back(-100.0);

	coeff_P.push_back(40.0 * 1000.0);
	coeff_P.push_back(400.0 * 1000.0);

	// -- pump and valve has to be added manually as a rigid subsystem --
  edges[idx_sziv1] = new Pump("sziv1","kozos2","n_sziv1",1000,Q_nom,H_nom,coeff_H,coeff_P,n_nom,n_act,theta,true);
  edges[idx_szel1] = new CheckValve("szel1","n_sziv1","kozos3",1000,zeta_open,zeta_closed,true);

  edges[idx_sziv1]->Add_transient("switch_off",900.,0.); //t_off=900, not simulated now


  //nodes
  int idx_kozos2 = reader.findNodeByID("kozos2");
  int idx_kozos3 = reader.findNodeByID("kozos3");
  int idx_nsziv1 = reader.findNodeByID("n_sziv1");

  //node: kozos2
  vector<int> edges_idx1{reader.findPipeByID("pk2"), idx_sziv1};
  vector<bool> edges_front1 {false, true};
  cons[idx_kozos2] = new Connector("kozos2",edges,edges_front1,edges_idx1,0,0);

  //node: kozos3
  vector<int> edges_idx2{idx_szel1, reader.findPipeByID("pk3")};
  vector<bool> edges_front2  {false, true};
  cons[idx_kozos3] = new Connector("kozos3",edges,edges_front2,edges_idx2,0,0);

  //node: n_sziv1
  vector<int> edges_idx3{idx_sziv1, idx_szel1};
  vector<bool> edges_front3{false, true};
  cons[idx_nsziv1] = new Connector("n_sziv1",edges,edges_front3,edges_idx3,0,0);

  //rigid subsystem def.
  vector< vector<int> > rigid_subsystems;
	vector<int> rs1; 
  rs1.push_back(idx_sziv1);
  rs1.push_back(idx_szel1);
  rigid_subsystems.push_back(rs1);

  /*
  ================== Running ==================
  */
  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end, rigid_subsystems);

  r.Set_DEBUG(false);
  r.Set_Folder(folder);
  r.Set_Save_data(true);
  r.Set_Save_interval(1.0);
  r.Set_Write_interval(300.0);

	//save pressure warnings above and below a threshold
  string fname = folder + "/pressure_warning";
  string type = "Maximum";
  r.Set_Pressure_limit(true,1.5e5,6.0e5,fname,type);
  r.Set_Pressure_limit_time(0,600);

  r.Set_Node_mul(1);
  r.Run(600.);
}



















