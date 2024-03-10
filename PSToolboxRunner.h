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
    void Set_Pressure_limit(bool set, double _p_lim, string _fname, string _type) {set_limit = set; p_lim = _p_lim; limits_fname=_fname; limits_type=_type;};
  private:
    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> cons;
    vector<int> con_at_edge_start;
    vector<int> con_at_edge_end;
    bool save_data;
    bool DEBUG;
    int Node_mul;
    double p_lim;
    bool set_limit;
    string limits_fname;
    string limits_type;
};
