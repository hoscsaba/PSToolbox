#include "PSToolboxRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include <fstream>

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
  p_limMin = 0.0;
  p_limMax = 0.0;
  write_interval = 1000.0;
  save_interval = 1.0;
  p_start_time = 0.0;
  p_end_time = 1.0e10;
};


void PSToolboxRunner::Run(double t_max){

  //create vectors for critical values
  vector<double> p_crit;
  vector<double> x_crit;
  vector<double> t_crit;
  vector<double> last_save;
  vector<string> ID_crit;


  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
  {
    edges.at(i)->Ini(Node_mul);
    last_save.push_back(0.0);
  }

  // Simulation
  double t_global=0., dt_out=t_max/100., t_out=-1.e-10;
  double t_next;
  double t_save = write_interval;
  int update_idx;
  vector<bool> update_edges(edges.size());
  int STEP=0;

  double T_TOL=1.e-10;
  /*for (unsigned int i=0; i<edges.size(); i++)
    if (edges.at(i)->Get_dt()<T_TOL)
      T_TOL=edges.at(i)->Get_dt();
  T_TOL/=100.;*/


  while (t_global<t_max)
  {
    /*
    ======================== TIMESTEP PREPARATIONS ===============================
    */
    STEP++;
    if (DEBUG)
    {
      cout<<endl<<"STEP     : "<<STEP;
      cout<<endl<<"t_global : "<<t_global;
    }
    if (t_global>t_out)
    {
      cout<<endl<<round(t_global/t_max*100)<<"%";
      t_out+=dt_out;
    }

    // Find smallest next target time
    fill(update_edges.begin(),update_edges.end(),false);

    update_idx=0;
    t_next=edges.at(update_idx)->Get_tnext();
    for (unsigned int i=0; i<edges.size(); i++)
    {
      if (DEBUG)
      {
        cout<<endl<<" pipe "<<i<<" ("<<edges.at(i)->Get_name()<<"): t="<<edges.at(i)->Get_t();
        cout<<", tnext="<<edges.at(i)->Get_tnext()<<", current: "<<t_next;
      }
      if (edges.at(i)->Get_tnext()<t_next)
      {
        update_idx=i;
        t_next=edges.at(i)->Get_tnext();

        if (DEBUG)
          cout<<"  <-- t_target ";
      }
    }
    update_edges.at(update_idx)=true;

    // Now find all elements within T_TOL
    for (unsigned int i=0; i<edges.size(); i++)
    {
      if (fabs((edges.at(i)->Get_tnext())-t_next)<T_TOL)
        update_edges.at(i)=true;
    }
    if (DEBUG)
    {
      cout<<endl<<endl<<"The following edges will be updated: ";
      for (unsigned int i=0; i<edges.size(); i++)
      {
        if (update_edges.at(i))
          cout<<" "<<i;
      }
    }


    /*
    ======================== UPDATE EDGES ===============================
    */
    for (unsigned int i=0; i<edges.size(); i++)
    { 
      if (update_edges.at(i))
      {
        if (DEBUG)
        {
          cout<<endl<<"\t Updating edge "<<i<<" ... ";
        }

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
          if(last_save[i] + save_interval < t_next)
          {
            edges.at(i)->Save_data();
            last_save[i] = int(t_next/save_interval)*save_interval;
          }
        }
          
        if(set_limit && t_next > p_start_time && t_next < p_end_time)
        {
          edges.at(i)->GetLargePressureValues(p_limMax,x_crit,p_crit,t_crit,ID_crit);
          edges.at(i)->GetSmallPressureValues(p_limMin,x_crit,p_crit,t_crit,ID_crit);
        }
      }//end of if - update edges
    }//end of for - update edges

    if (DEBUG)
    	cin.get();
    t_global=t_next;

    /*
    ======================== Save at end of timestep ===============================
    */

    if (save_data && t_global > t_save)
    {
      t_save = t_next + write_interval;
      for (unsigned int i=0; i<edges.size(); i++)
      {
          edges.at(i)->Write_data();
      }
    }

    if(t_global > 900) //turn off/on demands
    {
      //cons.at(172)->demand = 0.0808889;
      //cons.at(174)->demand = 0;
    }

  }//end of while 


  /*
  ======================== Save at end of simulation ===============================
  */
  if (save_data)
  {
    for (unsigned int i=0; i<edges.size(); i++)
    {
        edges.at(i)->Write_data();
    }
  }
      

  if(set_limit)
  {
    ofstream ofsMax = ofstream(limits_fname + "_max.log");
    ofstream ofsMin = ofstream(limits_fname + "_min.log");
    if(limits_type == "All")
    {
      for (int i = 0; i < ID_crit.size(); i++)
      {
        if (p_limMax < p_crit.at(i))
        {
          ofsMax << "Pressure limit exceeded in pipe " << ID_crit.at(i);
          ofsMax << "\t t = "<< t_crit.at(i) << " s";
          ofsMax << "\t x = "<< x_crit.at(i) << " m";
          ofsMax << "\t p = "<< 1e-5*p_crit.at(i) << " bar" << endl;
        }
        else
        {
          ofsMin << "Pressure limit exceeded in pipe " << ID_crit.at(i);
          ofsMin << "\t t = "<< t_crit.at(i) << " s";
          ofsMin << "\t x = "<< x_crit.at(i) << " m";
          ofsMin << "\t p = "<< 1e-5*p_crit.at(i) << " bar" << endl;
        }
      } 
    }
    if(limits_type == "Maximum")
    {
      for(int i = 0; i < edges.size(); i++) //go through all the pipes
      {
        int countMax = 0;
        int countMin = 0;
        double p_max = p_limMax;
        double p_min = p_limMin;
        int idxMax = -1;
        int idxMin = -1;
        for (int j = 0; j < ID_crit.size(); j++) //go through dataset
        {
          if(ID_crit.at(j) == edges.at(i)->Get_name() && p_crit.at(j) >= p_max)
          {
            p_max = p_crit.at(j);
            idxMax = j;
            countMax++;
          }
          if(ID_crit.at(j) == edges.at(i)->Get_name() && p_crit.at(j) <= p_min)
          {
            p_min = p_crit.at(j);
            idxMin = j;
            countMin++;
          }
        }
        if(idxMax >= 0)
        {
          ofsMax << "Pressure limit exceeded in pipe " << edges.at(i)->Get_name() << ", " << countMax << "\t times. Max: ";
          ofsMax << "\t t = "<< t_crit.at(idxMax)<< " s";
          ofsMax << "\t x = "<< x_crit.at(idxMax)<< " m";
          ofsMax << "\t p = "<< 1e-5*p_crit.at(idxMax) << " bar"<< endl;
        }
        else
        {
          ofsMax << "Pressure limit in pipe " << edges.at(i)->Get_name()  << "\tis NOT exceeded!" << endl;
        }

        if(idxMin >= 0)
        {
          ofsMin << "Pressure limit exceeded in pipe " << edges.at(i)->Get_name() << ", " << countMin << "\t times. Min: ";
          ofsMin << "\t t = "<< t_crit.at(idxMin)<< " s";
          ofsMin << "\t x = "<< x_crit.at(idxMin)<< " m";
          ofsMin << "\t p = "<< 1e-5*p_crit.at(idxMin) << " bar"<< endl;
        }
        else
        {
          ofsMin << "Pressure limit in pipe " << edges.at(i)->Get_name()  << "\tis NOT exceeded!" << endl;
        }
        
      }
    }
    
    ofsMax.flush();
    ofsMax.close();
    ofsMin.flush();
    ofsMin.close();
    
  }//end of pressure min-max writing

}//end of run function
