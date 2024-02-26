// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_dani.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "SCP.h"
#include "Connector.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "EpanetReader.h"
#include "Nyiregyhaza.h"

using namespace std;

int main(int argc, char **argv) {

  EpanetReader reader;
  //string location = "dummy_v1.inp";
  string location = "ipar_v2_without_pump_1.inp";
  reader.readFromFile(location);
  calculatePropagationVelocity(reader);   
  reader.convertToRunner2();

  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;

  PSToolboxRunner r(edges,cons,con_at_edge_end, con_at_edge_end);
  r.Set_Save_data(true);
  r.Run(10.);


}
