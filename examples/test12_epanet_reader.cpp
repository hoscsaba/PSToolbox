// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test12_epanet_reader.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "/Users/hoscsaba/program/PSToolbox/SCP.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxBaseEdge.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxPlotter.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxRunner.h"
#include "/Users/hoscsaba/program/PSToolbox/EpanetReader.h"

using namespace std;

int main(int argc, char **argv) {
  bool DEBUG=false;
  bool save_data=true;

 EpanetReader reader;
  string location = "dummy_v1.inp";
  //beolvasás
  reader.readFromFile(location);
  //peti adatai alapján Roughness -> Delta és a
  calculatePropagationVelocity(reader); 
  //ez a függvény hozza létre az edges, cons stb... dolgokat
  reader.convertToRunner2();

  //adatok kinyerése a reader osztályból
  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;


  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);
  r.Set_Save_data(true);
  r.Set_DEBUG(false);
  r.Set_Node_mul(5);
  r.Run(10.);

  PSToolboxPlotter pl1("p1.dat"); pl1.Plot();
  PSToolboxPlotter pl2("p2.dat"); pl2.Plot();
  PSToolboxPlotter pl3("p3.dat"); pl3.Plot();

}
