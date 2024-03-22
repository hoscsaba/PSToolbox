#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include "SCP.h"
#include "Connector.h"
#include "EpanetReader.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "Nyiregyhaza.h"
#include "Pump.h"
#include "CheckValve.h"

using namespace std;

int main(int argc, char **argv) 
{
  /*
  ================= Szimuláció beállítás =================
  */
  double Qipar = 278.0;
  string ipar = "Be2Ki";
  string kornyezet = "Max";
  string folder = "data_" + ipar + "_" + kornyezet + "_Q_" + to_string(int(Qipar)) + "_szivattyu";

  /*
  ================= Könyvtárak kezelése =================
  */
  string path = folder + "/*";
  string removecom = "rm " + path;
  int o = system(removecom.c_str());
  string mkdir = "mkdir " + folder;
  o = system(mkdir.c_str());

  /*
  ================ EPANET fájl beolvasása ================
  */
  EpanetReader reader;
  string location = "ipar_v4.inp";
  reader.readFromFile(location);
  calculatePropagationVelocity(reader);
  modifyDemand(reader,kornyezet);
  reader.convertToRunner2();
  
  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;

  /*
  ================ Szivattyú mókolás ================
  */ 
  bool save_data = true; 
  int idx_edge = edges.size();
  int idx_cons = cons.size();

  //nodeok
  string n_sziv_e = "29"; //szivattyú eleje
  string n_sziv_v = "30"; //szivattyú vége
  string n_tank_e = "tank_e"; //tank eleje
  string n_tank_v = "tank_v"; //tank vége, nem kell!

  //szivattyú paraméterek
  double n_nom = 1475;
  double n_act = 1475;
  double Q_nom = 632.09 / 3600.0;
  double H_nom = 16.65;
  double P_nom = 33.68*1000;
  double theta = 0.15; 

  //Pressure Reduce Valve
  double H_total = 161.0 - 113.0;
  double p_total = H_total * 9.81 * 1000;

  //számolt szivattyú dolgok
  double H_max=1.2*H_nom;
	double P_0=P_nom*0.2;
	vector<double> coeff_H, coeff_P;
	coeff_H.push_back(H_max);
	coeff_H.push_back(0.0);
	coeff_H.push_back((H_nom-H_max)/Q_nom/Q_nom);
	coeff_P.push_back(P_0);
	coeff_P.push_back((P_nom-P_0)/Q_nom);

  //3 szivattyú, 2 CV, 1 tank -> merev alrendszer
  double zeta_min = 1.0e-4;
  double zeta_max = 1.0e5;
  edges.push_back(new Pump("pump1",n_sziv_e,n_sziv_v,1000,Q_nom,H_nom,coeff_H,coeff_P,n_nom,n_act,theta,save_data)); //+0
  edges.push_back(new Pump("pump2",n_sziv_e,n_sziv_v,1000,Q_nom,H_nom,coeff_H,coeff_P,n_nom,n_act,theta,save_data)); //+1
  edges.push_back(new Pump("pump3",n_sziv_e,n_sziv_v,1000,Q_nom,H_nom,coeff_H,coeff_P,n_nom,n_act,theta,save_data)); //+2
  edges.push_back(new CheckValve("pump_bypass",n_sziv_e,n_sziv_v,1000,zeta_min,100*zeta_max,save_data)); //+3
  edges.push_back(new CheckValve("valve_pres",n_sziv_v,n_tank_e,1000,zeta_max,zeta_max,save_data)); //+4
  edges.push_back(new Tank("tank",n_tank_e,n_tank_v,1000,p_total,save_data)); //+5


  //szivattyú eleje connector
  int idx_of_sze = reader.findNodeByID(n_sziv_e);
  int idx_of_76 = reader.findPipeByID("76");
  cout << "Szivattyu eleje: " << idx_of_sze << "\t76-os cso indexe: " << idx_of_76 << endl;

  vector<int> edges_idx_sze {idx_edge + 0, idx_edge + 1, idx_edge + 2, idx_edge + 3, idx_of_76 };
  vector<bool> edges_front_sze  {true,         true        , true        , true        , false};
  cons[idx_of_sze] = new Connector(n_sziv_e,edges,edges_front_sze,edges_idx_sze,0.0,false);


  //szivattyú vége connector 
  int idx_of_szv = reader.findNodeByID(n_sziv_v);
  int idx_of_77 = reader.findPipeByID("77");
  cout << "Szivattyu  vege: " << idx_of_szv << "\t77-os cso indexe: " << idx_of_77 << endl;

  vector<int> edges_idx_szv {idx_edge + 0, idx_edge + 1, idx_edge + 2, idx_edge + 3, idx_of_77, idx_edge + 4 };
  vector<bool> edges_front_szv  {false,        false       , false       , false        ,true , true};
  cons[idx_of_szv] = new Connector(n_sziv_v,edges,edges_front_szv,edges_idx_szv,0.0,false);

  //tank connectorai
  vector<int> edges_idx_te {idx_edge + 4, idx_edge + 5};
  vector<bool> edges_front_te  {false,        true};
  cons.push_back(new Connector(n_tank_e,edges,edges_front_te,edges_idx_te,0.0,false));

  //rigid subsystem def.
  vector< vector<int> > rigid_subsystems;
	vector<int> rs1; 
  for (int i = 0; i < 6; i++)
  {
    rs1.push_back(idx_edge + i);
  }
  rigid_subsystems.push_back(rs1);

  //pump 1
  con_at_edge_start.push_back(reader.findNodeByID(n_sziv_e)); 
  con_at_edge_end.push_back(reader.findNodeByID(n_sziv_v)); 
  //pump 2
  con_at_edge_start.push_back(reader.findNodeByID(n_sziv_e)); 
  con_at_edge_end.push_back(reader.findNodeByID(n_sziv_v)); 
  //pump 3
  con_at_edge_start.push_back(reader.findNodeByID(n_sziv_e)); 
  con_at_edge_end.push_back(reader.findNodeByID(n_sziv_v)); 
  //check valve - bypass
  con_at_edge_start.push_back(reader.findNodeByID(n_sziv_e)); 
  con_at_edge_end.push_back(reader.findNodeByID(n_sziv_v)); 
  //check valve - pres.
  con_at_edge_start.push_back(reader.findNodeByID(n_sziv_v)); 
  con_at_edge_end.push_back(idx_cons); 
  //tank
  con_at_edge_start.push_back(idx_cons); 
  con_at_edge_end.push_back(idx_cons); 
  
  

  /*
  ================ DEBUG ============================
  */
    cout << "Size of edges:\t" << edges.size() << endl;
 for (int i = 0; i < edges.size(); i++)
  {
      cout << "i = " << i << "\t" << edges[i]->name << endl;
  }


  cout << " Size of cons:\t" << cons.size() << endl;
  for (int i = 0; i < cons.size(); i++)
  {
      cout << "i = " << i << "\t name: " << cons[i]->name <<"\tBC: " << cons[i]->BC_type << "\t size of edgeidx: "<< cons[i]->edges_idx.size() << endl;
      cout << "\tedges:";
      for (int j = 0; j < cons[i]->edges_idx.size(); j++)
      {
        cout << cons[i]->edges_idx[j] << "( "<< cons[i]->is_front[j]  <<"), "; 
      }
      cout << endl;
      
  }


  /*
  ================ Demand módosítás ================
  */
  RuntimeModifier iparBE, iparKI, korny;

  // -------------- IPAR -----------------
  iparKI.Add(reader.findNodeByID("n_DIP"),0.0 / 3600.0);
  iparBE.Add(reader.findNodeByID("n_DIP"),Qipar / 3600.0);
  if(ipar == "Be2Ki")
  {
    iparBE.time = 0.0;
    iparKI.time = 900.0;
  }
  else if(ipar == "Ki2Be")
  {
    iparBE.time = 900.0;
    iparKI.time = 0.0;
  }
  else
  {
    cout << "Ismeretlep IPAR allapotvaltozas" << endl;
    while(1) 1;
  }

  // ----------- KÖRNYEZET ---------------
  korny.time = 0.0;
  double R2_old = 145.46;
    //------------------- MAX -----------------------
  if(kornyezet == "Max" && ipar == "Be2Ki" && Qipar == 417.0)
  {
    korny.Add(reader.findNodeByID("r2"),(145.14 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"), -7.27 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),-16.52 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 91.93 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"), 17.56 / 3600.0);
  }
  else if(kornyezet == "Max" && ipar == "Be2Ki" && Qipar == 278.0)
  {
    korny.Add(reader.findNodeByID("r2"),(145.59 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"),  4.56 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"), -8.66 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 95.11 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"), 52.91 / 3600.0);
  }
  else if(kornyezet == "Max" && ipar == "Be2Ki" && Qipar == 208.0)
  {
    korny.Add(reader.findNodeByID("r2"),(145.80 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"),  8.53 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"), -4.38 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 91.44 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"), 74.04 / 3600.0);
  }
  else if(kornyezet == "Max" && ipar == "Ki2Be")
  {
    korny.Add(reader.findNodeByID("r2"),(146.37 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"), 12.37 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),  7.66 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 80.98 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"),   136 / 3600.0);
  } 
  
  
  //------------------- MIN -----------------------
  else if(kornyezet == "Min" && ipar == "Be2Ki" && Qipar == 417.0)
  {
    korny.Add(reader.findNodeByID("r2"),(146.42 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"), -5.58 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),-17.82 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 54.86 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"),-147.57/ 3600.0);
  }
  else if(kornyezet == "Min" && ipar == "Be2Ki" && Qipar == 278.0)
  {
    korny.Add(reader.findNodeByID("r2"),(146.69 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"),  4.81 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),-11.54 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 49.65 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"),-108.77/ 3600.0);
  }
  else if(kornyezet == "Min" && ipar == "Be2Ki" && Qipar == 208.0)
  {
    korny.Add(reader.findNodeByID("r2"),(146.80 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"),  7.34 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),  8.14 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 46.74 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"),-89.72 / 3600.0);
  }
  else if(kornyezet == "Min" && ipar == "Ki2Be")
  {
    korny.Add(reader.findNodeByID("r2"),(147.19 - R2_old)*1000*9.81);
    korny.Add(reader.findNodeByID("n1208"), 10.02 / 3600.0);
    korny.Add(reader.findNodeByID("n1568"),  3.54 / 3600.0);
    korny.Add(reader.findNodeByID("n1181"), 36.04 / 3600.0);
    korny.Add(reader.findNodeByID("n3040"), -22.86/ 3600.0);
  }
  else
  {
    cout << "Ismeretlep KORNYEZET beallitas" << endl;
    while(1) 1;
  }

  /*
  ================== Futtatás ==================
  */
  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end, rigid_subsystems);

  r.Set_DEBUG(false);
  r.Set_Folder(folder);
  r.Set_Save_data(true);
  r.Set_Save_interval(1.0);
  r.Set_Write_interval(50.0);
  r.Add_Runtime_modifier(iparBE);
  r.Add_Runtime_modifier(iparKI);
  r.Add_Runtime_modifier(korny);

  string fname = folder + "/pressure_warning";
  string type = "Maximum";
  r.Set_Pressure_limit(true,1.0e5,8.0e5,fname,type);
  r.Set_Pressure_limit_time(500,1800);

  r.Set_Node_mul(1);
  r.Run(1800.);
}
