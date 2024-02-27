#include "PSToolboxRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

PSToolboxRunner::PSToolboxRunner(
    vector<PSToolboxBaseEdge *> &_e, 
    vector<Connector *> &_c,
    vector<int> _con_at_edge_start,
    vector<int> _con_at_edge_end){
  edges = _e;
  cons = _c;
  con_at_edge_start = _con_at_edge_start;
  con_at_edge_end = _con_at_edge_end;
  save_data=false;
  DEBUG=false;
  Node_mul=1;
};


void PSToolboxRunner::Run(double t_max){

  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
    edges.at(i)->Ini(Node_mul);

  // Simulation
  double t_global=0., dt_out=t_max/10., t_out=-1.e-10;
  double t_next;
  int update_idx;
  vector<bool> update_edges(edges.size());
  int STEP=0;

  double T_TOL=1.e-3;
  for (unsigned int i=0; i<edges.size(); i++)
    if (edges.at(i)->Get_dt()<T_TOL)
      T_TOL=edges.at(i)->Get_dt();
  T_TOL/=100.;


  while (t_global<t_max){
    STEP++;
    if (DEBUG){
      cout<<endl<<"STEP     : "<<STEP;
      cout<<endl<<"t_global : "<<t_global;
    }
    if (t_global>t_out){
      cout<<endl<<round(t_global/t_max*100)<<"%";
      t_out+=dt_out;
    }

    // Find smallest next target time
    fill(update_edges.begin(),update_edges.end(),false);

    update_idx=0;
    t_next=edges.at(update_idx)->Get_tnext();
    for (unsigned int i=0; i<edges.size(); i++){
      if (DEBUG){
        cout<<endl<<" pipe "<<i<<" ("<<edges.at(i)->Get_name()<<"): t="<<edges.at(i)->Get_t();
        cout<<", tnext="<<edges.at(i)->Get_tnext()<<", current: "<<t_next;
      }
      if (edges.at(i)->Get_tnext()<t_next){
        update_idx=i;
        t_next=edges.at(i)->Get_tnext();

        if (DEBUG)
          cout<<"  <-- t_target ";
      }
    }
    update_edges.at(update_idx)=true;

    // Now find all elements within T_GLOBAL_TOL
    for (unsigned int i=0; i<edges.size(); i++)
      if (fabs((edges.at(i)->Get_tnext())-t_next)<T_TOL)
        update_edges.at(i)=true;
    if (DEBUG){
      cout<<endl<<endl<<"The following edges will be updated: ";
      for (unsigned int i=0; i<edges.size(); i++)
        if (update_edges.at(i))
          cout<<" "<<i;
    }

    for (unsigned int i=0; i<edges.size(); i++){ 
      if (update_edges.at(i)){
        if (DEBUG)
          cout<<endl<<"\t Updating edge "<<i<<" ... ";
        edges.at(i)->UpdateInternal();

        cons.at(con_at_edge_start.at(i))->Update(t_next,i);

        cons.at(con_at_edge_end.at(i))->Update(t_next,i);

        edges.at(i)->UpdateTime(t_next);
        if (DEBUG)
          cout<<"done.";
        if (save_data)
          edges.at(i)->Save_data();
      }
    }
    if (DEBUG)
      cin.get();

    t_global=t_next;
  }
  if (save_data)
    for (unsigned int i=0; i<edges.size(); i++)
      edges.at(i)->Write_data();

}
