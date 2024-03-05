#include "PSToolboxRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

PSToolboxRunner::PSToolboxRunner(
    vector<PSToolboxBaseEdge *> &_e, 
    vector<Connector *> &_c,
    vector<int> &_con_at_edge_start,
    vector<int> &_con_at_edge_end){
  edges = _e;
  cons = _c;
  con_at_edge_start = _con_at_edge_start;
  con_at_edge_end = _con_at_edge_end;
  save_data=false;
  DEBUG=false;
  Node_mul=1;
  set_limit = false;
  p_lim = 0.0;
};


void PSToolboxRunner::Run(double t_max){

  //create vectors for critical values
  vector<double> p_crit;
  vector<double> x_crit;
  vector<double> t_crit;
  vector<string> ID_crit;


  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
    edges.at(i)->Ini(Node_mul);

  // Simulation
  double t_global=0., dt_out=t_max/100., t_out=-1.e-10;
  double t_next;
  int update_idx;
  vector<bool> update_edges(edges.size());
  int STEP=0;

  double T_TOL=1.e-10;
  /*for (unsigned int i=0; i<edges.size(); i++)
    if (edges.at(i)->Get_dt()<T_TOL)
      T_TOL=edges.at(i)->Get_dt();
  T_TOL/=100.;*/


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

    // Now find all elements within T_TOL
    for (unsigned int i=0; i<edges.size(); i++)
      if (fabs((edges.at(i)->Get_tnext())-t_next)<T_TOL)
        update_edges.at(i)=true;
    if (DEBUG){
      cout<<endl<<endl<<"The following edges will be updated: ";
      for (unsigned int i=0; i<edges.size(); i++)
        if (update_edges.at(i))
          cout<<" "<<i;
    }

    for (unsigned int i=0; i<edges.size(); i++)
    { 
      if (update_edges.at(i))
      {
        if (DEBUG)
          cout<<endl<<"\t Updating edge "<<i<<" ... ";
        edges.at(i)->UpdateInternal();

        cons.at(con_at_edge_start.at(i))->Update(t_next,i);

        cons.at(con_at_edge_end.at(i))->Update(t_next,i);

        if (DEBUG)
        {  
          cout<<"done.";
          cout << edges.at(i)->Info();
        }

        edges.at(i)->UpdateTime(t_next);

        if (save_data)
        {
          if(edges.at(i)->Get_name() == "105" || edges.at(i)->Get_name() == "3" || edges.at(i)->Get_name() == "74" || edges.at(i)->Get_name() == "p1100" || edges.at(i)->Get_name() == "24" )
            edges.at(i)->Save_data();
        }
          

        if(set_limit)
        {
          edges.at(i)->GetLargePressureValues(p_lim,x_crit,p_crit,t_crit,ID_crit);
        }
      }
    }
    if (DEBUG)
      cin.get();

    t_global=t_next;
  }
  if (save_data)
    for (unsigned int i=0; i<edges.size(); i++)
    {
      if(edges.at(i)->Get_name() == "105" || edges.at(i)->Get_name() == "3" || edges.at(i)->Get_name() == "74" || edges.at(i)->Get_name() == "p1100" || edges.at(i)->Get_name() == "24" )
        edges.at(i)->Write_data();
    }
      

  if(set_limit)
  {
    ofstream ofs = ofstream(limits_fname);
    if(limits_type == "All")
    {
      for (int i = 0; i < ID_crit.size(); i++)
      {
          ofs << "Pressure limit exceeded in pipe " << ID_crit.at(i);
          ofs << "\t t = "<< t_crit.at(i) << " s";
          ofs << "\t x = "<< x_crit.at(i) << " m";
          ofs << "\t p = "<< 1e-5*p_crit.at(i) << " bar" << endl;
      } 
    }
    if(limits_type == "Maximum")
    {
      for(int i = 0; i < edges.size(); i++) //go through all the pipes
      {
        int count = 0;
        double p_max = p_lim;
        int idx = -1;
        for (int j = 0; j < ID_crit.size(); j++) //go through dataset
        {
          if(ID_crit.at(j) == edges.at(i)->Get_name() && p_crit.at(j) >= p_max)
          {
            p_max = p_crit.at(j);
            idx = j;
            count++;
          }
        }
        if(idx >= 0)
        {
          ofs << "Pressure limit exceeded in pipe " << ID_crit.at(idx) << ", " << count << "\t times. Max: ";
          ofs << "\t t = "<< t_crit.at(idx)<< " s";
          ofs << "\t x = "<< x_crit.at(idx)<< " m";
          ofs << "\t p = "<< 1e-5*p_crit.at(idx) << " bar"<< endl;
        }
        else
        {
          ofs << "Pressure limit in pipe " << edges.at(i)->Get_name()  << "\tis NOT exceeded!" << endl;
        }
        
      }
    }
    
    ofs.flush();
    ofs.close();
    
  }

}
