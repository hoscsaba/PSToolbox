#pragma once

#include <vector>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

class PSToolboxRunner
{
  public:
    PSToolboxRunner(vector<PSToolboxBaseEdge *>&, vector<Connector *>&,vector<int>&,vector<int>&);
    void Run(double);
    void Set_Save_data(bool newval){save_data=newval;};
    void Set_DEBUG(bool newval){DEBUG=newval;};
    void Set_Node_mul(int newval){Node_mul=newval;};
    void Set_Pressure_limit(bool set, double _p_limMin, double _p_limMax, string _fname, string _type) {set_limit = set; p_limMin = _p_limMin; p_limMax = _p_limMax; limits_fname=_fname; limits_type=_type;};
    void Set_Pressure_limit_time(double start, double end) {p_start_time = start; p_end_time = end;};
    void Set_Save_interval(double interval) {save_interval = interval;};
    void Set_Write_interval(double interval) {write_interval = interval;};
  private:
    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> cons;
    vector<int> con_at_edge_start;
    vector<int> con_at_edge_end;
    bool DEBUG;
    int Node_mul;

    //min-max pressure
    double p_limMin;
    double p_limMax;
    double p_start_time;
    double p_end_time;
    
    bool set_limit;
    string limits_fname;
    string limits_type;

    //save stuff
    bool save_data;
    double write_interval;
    double save_interval;


};
