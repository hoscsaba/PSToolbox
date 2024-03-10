// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_Csaba.cpp

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

using namespace std;

int main(int argc, char **argv) {

  //remove data directory
  string path = "data/*";
  string removecom = "rm " + path;
  system(removecom.c_str());


  EpanetReader reader;
  string location = "ipar_v3_new.inp";
  reader.readFromFile(location);
  calculatePropagationVelocity(reader); 
  //calculateLambda(reader); 
  cout << "T1: " << reader.junctions[10].Demand;
  modifyDemand(reader);
  cout << "T1: " << reader.junctions[10].Demand;
  reader.convertToRunner2();
  
  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;




  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);


  string fname = "data/pressure_warning";
  string type = "All";
  r.Set_Pressure_limit(true,1.0e5,10.0e5,fname,type);

  r.Set_DEBUG(false);
  r.Set_Node_mul(1);
  r.Set_Save_data(true);
  r.Run(1800.);


}
