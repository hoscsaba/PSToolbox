// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_Csaba.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include "SCP.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "my_tools.h"
#include "Nyiregyhaza.h"

using namespace std;

int main(int argc, char **argv) {

  //remove data directory
  string path = "data/*";
  string removecom = "rm " + path;
  system(removecom.c_str());

  EpanetReader reader;
  string location = "dummy_v02.inp";
  reader.readFromFile(location);
  //calculatePropagationVelocity(reader); 
  //calculateLambda(reader); 
  for (int i = 0; i < 4; i++)
    {
    reader.pipes[i].Lambda = 0.03;
    reader.pipes[i].SpeedOfSound = 100.0;
  }
  reader.pipes[1].SpeedOfSound = 121.0;
  reader.pipes[2].SpeedOfSound = 96.0;
  


  reader.convertToRunner2();
  //reader.PrintData();
  
  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;




  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);
  r.Set_Save_data(true);
  r.Set_DEBUG(false);
  r.Set_Node_mul(1);
  r.Run(100.);


}
