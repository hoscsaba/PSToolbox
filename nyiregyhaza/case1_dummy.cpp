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
  string location = "verification.inp";
  reader.readFromFile(location);
  //calculatePropagationVelocity(reader); 
  //calculateLambda(reader); 
  for (int i = 0; i < reader.pipes.size(); i++)
    {
    reader.pipes[i].SpeedOfSound = 103.0;
  }
  


  reader.convertToRunner2();
  
  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;




  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);
  r.Set_Save_data(true);
  r.Set_DEBUG(false);
  r.Set_Node_mul(2);

  //setup pressure limit warnings
  string fname = "data/pressure_warning";
  string type = "Maximum";
  r.Set_Pressure_limit(true,5.0e5,12.0e5,fname,type);
  
  r.Run(500.);


}
