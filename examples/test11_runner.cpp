// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lpython2.7 -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long test11_runner.cpp

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include "/Users/hoscsaba/program/PSToolbox/SCP.h"
#include "/Users/hoscsaba/program/PSToolbox/Connector.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxBaseEdge.h"
#include "/Users/hoscsaba/program/PSToolbox/PSToolboxPlotter.h"

using namespace std;

int main(int argc, char **argv) {
  bool DEBUG=false;
  vector<PSToolboxBaseEdge *> edges;
  vector<Connector *> cons;

  bool save_data=true;

  //=============== THIS SECTION GOES TO THE DATA READER
  // define edges
  edges.push_back(new SCP("p1","n1","n2",1000,1300,100,0.1,0.0,0,0,save_data));
  edges.push_back(new SCP("p2","n2","n3",1000,1300,30,0.1,0.0,0,0,save_data));
  //edges.push_back(new SCP("p3","n2","n4",1000,900,100,0.1,0.02,0,0,save_data));

  // define nodes (connectors)
  double demand1 = 0.;
  double demand2 = 0.;
  double demand3 = 0.;
  double demand4 = 0.;
  double H1 = 0.;
  double H2 = 0.;
  double H3 = 0.;
//  vector<int> id1{0};
  vector<int> id1;
  id1.push_back(0);
  id1.push_back(1);
//  vector<int> id1{1};
  cons.push_back(new Connector("n1",H1,edges.at(0),true,"Pressure",3.e5,demand1,DEBUG));
  cons.push_back(new Connector("n2",H2,edges.at(0),false,edges.at(1),true,demand2,DEBUG,id1));
  //cons.push_back(new Connector("n2",edges.at(0),false,edges.at(1),true,edges.at(2),true,demand2,DEBUG));
  cons.push_back(new Connector("n3",H3,edges.at(1),false,"Pressure",1.e5,demand3,DEBUG));
  //  cons.push_back(new Connector("n4",edges.at(2),false,"Velocity",0.,demand4,DEBUG));

  // We need to add here connectivity info, e.g. by adding edge pointer + start/end info if relevant
  // con_at_edge_start.at(i) stores the idx of the connector connected to the start of the edge 
  // con_at_edge_end.at(i)   stores the idx of the connector connected to the end of the edge 
  vector<int> con_at_edge_start(edges.size(),-1);
  vector<int> con_at_edge_end(edges.size(),-1);
  con_at_edge_start.at(0)=0; con_at_edge_end.at(0)=1;
  con_at_edge_start.at(1)=1; con_at_edge_end.at(1)=2;
  //con_at_edge_start.at(2)=1; con_at_edge_end.at(2)=3;
  // =============== END OF DATA READER SECTION

  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
    edges.at(i)->Ini();

  // Simulation
  double t_global=0.,t_max=1., dt_out=t_max/10., t_out=-1.e-10;
  double t_next;
  int update_idx;
  vector<bool> update_edges(edges.size());


  std::ofstream outputFile("output.txt");

  // Check if the file is opened successfully
  if (!outputFile.is_open()) {
    std::cerr << "Error opening the file!" << std::endl;
    return 1; // Return an error code
  }

  while (t_global<t_max){
    outputFile<<t_global<<endl;
    if (t_global>t_out){
      cout<<endl<<round(t_global/t_max*100)<<"%";
      t_out+=dt_out;
    }
    // Find the edge that needs update
    // TODO: if two edges are very close to each other, we need to update them at once

    fill(update_edges.begin(),update_edges.end(),false);

    update_idx=0;
    t_next=edges.at(update_idx)->Get_tnext();
    for (unsigned int i=1; i<edges.size(); i++){
      //        cout<<endl<<" pipe "<<i<<" ("<<edges.at(i)->Get_name()<<"): t="<<edges.at(i)->Get_t();
      //          cout<<", tnext="<<edges.at(i)->Get_tnext()<<", current: "<<t_next;
      if (edges.at(i)->Get_tnext()<t_next){
        update_idx=i;
        t_next=edges.at(i)->Get_tnext();
      }
    }

    update_edges.at(update_idx)=true;

    edges.at(update_idx)->Step();
    cons.at(con_at_edge_start.at(update_idx))->Update(t_next,update_idx);
    cons.at(con_at_edge_end.at(update_idx))->Update(t_next,update_idx);
    edges.at(update_idx)->UpdateTime(t_next);
/*
    if (update_idx==1){
      cout<<edges.at(1)->Info();
      cin.get();
    }
*/
    //}
    t_global=t_next;

    if (save_data)
      for (unsigned int i=0; i<edges.size(); i++)
        edges.at(i)->Save_data();
}

outputFile.close();

if (save_data)
  for (unsigned int i=0; i<edges.size(); i++)
  edges.at(i)->Write_data();

  //cout<<edges.at(0)->Info();
  //cout<<edges.at(1)->Info();

  PSToolboxPlotter pl1("p1.dat"); pl1.Plot();
  PSToolboxPlotter pl2("p2.dat"); pl2.Plot();
  //PSToolboxPlotter pl3("p3.dat"); pl3.Plot();



  }
