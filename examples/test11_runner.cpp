// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner.cpp

#include <stdio.h>
#include "SCP.h"
#include "Connector.h"
#include "Reservoir.h"
#include "Valve.h"
#include "Valve_with_Absorber.h"
#include "PSToolboxBaseEdge.h"
#include "PSToolboxRunner.h"
#include "my_tools.h"
#include "EpanetReader.h"

using namespace std;

int main(int argc, char **argv) {
  bool DEBUG=false;

  //=============== THIS SECTION GOES TO THE DATA READER
<<<<<<< Updated upstream
  /*/ define edges
  edges.push_back(new SCP("p1","n1","n2",1000,1300,100,0.1,0.02,0,0,false));
  edges.push_back(new SCP("p2","n2","n3",1000,1100,100,0.1,0.02,0,0,false));
=======
  // define edges
  edges.push_back(new SCP("p1","n1","n2",1000,1100,100,0.1,0.02,0,0,save_data));
  edges.push_back(new SCP("p2","n2","n3",1000,1300, 70,0.1,0.02,0,0,save_data));
  edges.push_back(new SCP("p3","n2","n4",1000, 900,100,0.1,0.02,0,0,save_data));
>>>>>>> Stashed changes

  // define nodes (connectors)
  double demand1 = 0.;
  double demand2 = 0.;
  double demand3 = 0.;
<<<<<<< Updated upstream
  cons.push_back(new Connector(edges.at(0),true,"Pressure",1.e5,demand1,DEBUG));
  cons.push_back(new Connector(edges.at(0),false,edges.at(1),false,demand2,DEBUG));
  cons.push_back(new Connector(edges.at(1),true,"Velocity",0.,demand3,DEBUG));
=======
  double demand4 = 0.;
  vector<int> id1;
  id1.push_back(0);
  id1.push_back(1);
  id1.push_back(2);
  //  vector<int> id1{1};
  cons.push_back(new Connector("n1",edges.at(0),true,"Pressure",3.e5,demand1,DEBUG));
  cons.push_back(new Connector("n2",edges.at(0),false,edges.at(1),true,edges.at(2),true,demand2,DEBUG,id1));
  cons.push_back(new Connector("n3",edges.at(1),false,"Pressure",1.e5,demand3,DEBUG));
  cons.push_back(new Connector("n4",edges.at(2),false,"Pressure",2.e5,demand4,DEBUG));
>>>>>>> Stashed changes

  // We need to add here connectivity info, e.g. by adding edge pointer + start/end info if relevant
  // con_at_edge_start.at(i) stores the idx of the connector connected to the start of the edge 
  // con_at_edges_end.at(i)   stores the idx of the connector connected to the end of the edge 
  vector<int> con_at_edge_start(edges.size(),-1);
  vector<int> con_at_edge_end(edges.size(),-1);
  con_at_edge_start.at(0)=0; con_at_edge_end.at(0)=1;
  con_at_edge_start.at(1)=1; con_at_edge_end.at(1)=2;
<<<<<<< Updated upstream
  */
  EpanetReader reader;
  string location = "dummy_v1.inp";
  reader.readFromFile(location);
  reader.convertToRunner2();

  //végén ID1 
  //connector(ID,magassag,..)


  vector<PSToolboxBaseEdge *> edges = reader.edges;
  vector<Connector *> cons = reader.cons;
  vector<int> con_at_edge_start = reader.con_at_edge_start;
  vector<int> con_at_edge_end = reader.con_at_edge_end;
 
  return 0;

=======
  con_at_edge_start.at(2)=1; con_at_edge_end.at(2)=3;
>>>>>>> Stashed changes
  // =============== END OF DATA READER SECTION

  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
    edges.at(i)->Ini();

  // Simulation
  double t_global=0.,t_max=1.;
  double t_next;
  int update_idx;
  vector<bool> update_edges(edges.size());

<<<<<<< Updated upstream
  while (t_global<t_max){
=======
  //  std::ofstream outputFile("output.txt");

  // Check if the file is opened successfully
  //if (!outputFile.is_open()) {
  //  std::cerr << "Error opening the file!" << std::endl;
  //  return 1; // Return an error code
  // }

  while (t_global<t_max){
    //outputFile<<t_global<<endl;
    if (t_global>t_out){
      cout<<endl<<round(t_global/t_max*100)<<"%";
      t_out+=dt_out;
    }
>>>>>>> Stashed changes
    // Find the edge that needs update
    // TODO: if two edges are very close to each other, we need to update them at once
    for (unsigned int i=0; i<edges.size(); i++){
      fill(update_edges.begin(),update_edges.end(),false);
      t_next=t_global+1.e5;
      for (unsigned int i=0; i<edges.size(); i++)
        if (edges.at(i)->Get_tnext()<t_next){
          update_idx=i;
          t_next=edges.at(i)->Get_tnext();
        }
    }
    update_edges.at(update_idx)=true;

<<<<<<< Updated upstream
    // Perform update. 
    for (unsigned int i=0; i<edges.size(); i++)
      if (update_edges.at(i)){
        // This updates internal points only, without any info on the BCs
        edges.at(i)->Step();
        // Take care of the BCs
        cons.at(con_at_edge_start.at(update_idx))->Update(t_next);
        cons.at(con_at_edge_end.at(update_idx))->Update(t_next);
      }
    t_global=t_next;
    cout<<endl<<t_global/t_max<<" -> update: "<<update_idx;
  }
=======
    edges.at(update_idx)->Step();
    cons.at(con_at_edge_start.at(update_idx))->Update(t_next,update_idx);
    cons.at(con_at_edge_end.at(update_idx))->Update(t_next,update_idx);
    edges.at(update_idx)->UpdateTime(t_next);

    t_global=t_next;

    if (save_data)
      for (unsigned int i=0; i<edges.size(); i++)
        edges.at(i)->Save_data();
  }

  //outputFile.close();

  if (save_data)
    for (unsigned int i=0; i<edges.size(); i++)
      edges.at(i)->Write_data();

  PSToolboxPlotter pl1("p1.dat"); pl1.Plot();
  PSToolboxPlotter pl2("p2.dat"); pl2.Plot();
  PSToolboxPlotter pl3("p3.dat"); pl3.Plot();
>>>>>>> Stashed changes
}
