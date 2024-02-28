// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner_Csaba.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "/Users/hoscsaba/program/PSToolbox/SCP.h"
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
  edges.push_back(new SCP("p1","n1","n2",1000,1300,100,0.1,0.02,0,0,save_data));
  edges.push_back(new SCP("p2","n2","n3",1000,1300,101,0.1,0.02,0,0,save_data));
  edges.push_back(new SCP("p3","n3","n4",1000,1300, 99,0.1,0.02,0,0,save_data));

  // define nodes (connectors)
  double demand1 = 0.;
  double demand2 = 0.;
  double demand3 = 0.;
  double demand4 = 0.;
  vector<int> id1;
  id1.push_back(0);
  id1.push_back(1);
  id1.push_back(2);
 // vector<int> id2;
 // id2.push_back(1);
 // id2.push_back(2);
  cons.push_back(new Connector("n1",edges.at(0),true,"Pressure",1.e5,demand1,DEBUG));
  cons.push_back(new Connector("n2",edges.at(0),false,edges.at(1),true,edges.at(2),true,demand2,DEBUG,id1));
// cons.push_back(new Connector("n3",edges.at(0),false,edges.at(1),true,demand2,DEBUG,id1));
  cons.push_back(new Connector("n4",edges.at(1),false,"Pressure",3.e5,demand3,DEBUG));
  cons.push_back(new Connector("n4",edges.at(2),false,"Pressure",1.e5,demand3,DEBUG));

  // We need to add here connectivity info, e.g. by adding edge pointer + start/end info if relevant
  // con_at_edge_start.at(i) stores the idx of the connector connected to the start of the edge 
  // con_at_edge_end.at(i)   stores the idx of the connector connected to the end of the edge 
  vector<int> con_at_edge_start(edges.size(),-1);
  vector<int> con_at_edge_end(edges.size(),-1);
  con_at_edge_start.at(0)=0; con_at_edge_end.at(0)=1;
  con_at_edge_start.at(1)=1; con_at_edge_end.at(1)=2;
  con_at_edge_start.at(2)=1; con_at_edge_end.at(2)=3;
  // =============== END OF DATA READER SECTION

  PSToolboxRunner r(edges,cons,con_at_edge_start, con_at_edge_end);
  r.Set_Save_data(true);
  r.Set_DEBUG(DEBUG);
  r.Set_Node_mul(2);
  r.Run(10.);

  PSToolboxPlotter pl1("p1.dat"); pl1.Plot();
  PSToolboxPlotter pl2("p2.dat"); pl2.Plot();
  PSToolboxPlotter pl3("p3.dat"); pl3.Plot();

}
