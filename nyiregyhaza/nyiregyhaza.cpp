// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_Csaba.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "SCP.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "my_tools.h"
#include "Nyiregyhaza.h"

using namespace std;

int main(int argc, char **argv) {


  EpanetReader reader;
  string location = "ipar_v3.inp";
  reader.readFromFile(location);
  calculatePropagationVelocity(reader); 
  calculateLambda(reader); 
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
  r.Set_Node_mul(3);
  r.Run(100.);


}
