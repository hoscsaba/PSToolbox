#pragma once

#include "Reservoir.h"
#include "LWP.h"
#include "SCP.h"
#include "Valve.h"
#include "ThrottleValve.h"
#include "Damper.h"
#include "Valve_with_Absorber.h"
#include "PSToolboxBaseEdge.h"

//! Class for connecting elements
/*!
  The class contains data for Connectors
  */
class Connector {
public:
    Connector(bool DEBUG);

    //type=0
    Connector(
        string name,
        vector<PSToolboxBaseEdge *> &_e,
        vector<bool> &_is_front,
        vector<int> &_edges_idx,
        double _demand, bool _DEBUG);

    // type=1
    Connector(
        string name,
        PSToolboxBaseEdge *e1, bool is_front,
        string BC_type, double BC_value,
        double demand, bool DEBUG);

    // type=2
    Connector(
        string name,
        PSToolboxBaseEdge *e1, bool is_front1,
        PSToolboxBaseEdge *e2, bool is_front2,
        double demand, bool DEBUG, vector<int> &edges_idx);

    // type=3
    Connector(
        string name,
        PSToolboxBaseEdge *e1, bool is_front1,
        PSToolboxBaseEdge *e2, bool is_front2,
        PSToolboxBaseEdge *e3, bool is_front3,
        double demand, bool DEBUG, vector<int> &edges_idx);

    // type=4
    Connector(
        string name,
        PSToolboxBaseEdge *e1, bool is_front1,
        PSToolboxBaseEdge *e2, bool is_front2,
        PSToolboxBaseEdge *e3, bool is_front3,
        PSToolboxBaseEdge *e4, bool is_front4,
        double demand, bool DEBUG, vector<int> &edges_idx);


    ~Connector();

    void Update(double t_target, int update_idx);

    void Connector_Reservoir_and_Valve_Inlet(double t_target,
                                             Reservoir *r, Valve *v, bool INLET_PRESSURE_DROP, double pd, double &p,
                                             double &mp);

    void Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
                                                Reservoir *r, LWP *p, bool inlet_pressure_drop, double &pf, double &Tf);

    bool Connector_LWP_Reservoir_and_Pipe_Front_inlet(double t_target,
                                                      Reservoir *r, LWP *p, bool inlet_pressure_drop, double &pf,
                                                      double &Tf);

    void Connector_LWP_Reservoir_and_Pipe_Front(double t_target,
                                                Reservoir *r, LWP *p, bool inlet_pressure_drop);

    void Connector_LWP_Reservoir_and_Pipe_Back(double t_target,
                                               Reservoir *r, LWP *p, bool inlet_pressure_drop);

    bool Connector_LWP_Pipe_Back_and_Valve(double t_target,
                                           LWP *p, Valve *v, double p_downstream, double &pb, double &Tb);

    void Connector_LWP_Pipe_Back_and_Valve(double t_target,
                                           LWP *p, Valve *v, double p_downstream);

    void Connector_LWP_Pipe_Front_and_Valve(double t_target,
                                            LWP *p, Valve *v, double p_upstream, double T_upstream);

    void Connector_LWP_Pipe_Back_Valve_Front(double t_target,
                                             LWP *p1, Valve *v1, LWP *p2);

    bool Connector_LWP_Pipe_Back_and_Valve_with_Absorber(double t_target,
                                                         LWP *p, Valve_with_Absorber *v, double p_downstream,
                                                         double &pb, double &Tb);

    void Connector_SCP_Reservoir_and_Pipe_Front(double t_target,
                                                Reservoir *r, SCP *p, double rho, double a, bool inlet_pressure_drop,
                                                double &pf);

    bool Connector_SCP_Pipe_Back_and_Valve(double t_target,
                                           SCP *p, Valve *v, double rho, double a, double p_downstream, double &pb);

    bool Connector_SCP_Pipe_Front_and_Valve_Outlet(double t_target,
                                                   SCP *p, Valve *v, double p_upstream, double rho, double a,
                                                   double &pb);

    void Connector_SCP_Pipes(double t_target);

    void Connector_SCP_Pipe_Simple_BC(double t_target);

    void Connector_SCP_Pipes(double t_target,
                             PSToolboxBaseEdge *p1, bool is_front1,
                             PSToolboxBaseEdge *p2, bool is_front2,
                             double mpout,
                             double &p, double &v1, double &v2);

    void Connector_SCP_2Pipes(double t_target, int update_idx);

    void Connector_SCP_3Pipes(double t_target, int update_idx);

    void Connector_SCP_4Pipes(double t_target, int update_idx);

    void Connector_SCP_NPipes(double t_target, int update_idx);

    void Connector_LWP_Pipes(double t_target, LWP *p1, LWP *p2);

    void Connector_LWP_Pipes(double t_target,
                             LWP *p1, bool is_front1,
                             LWP *p2, bool is_front2,
                             LWP *p3, bool is_front3);

    void Connector_LWP_ThrottleValve_LWP(double t_target,
                                         LWP *p1, ThrottleValve *tv, LWP *p2);

    void Update(double t_target, LWP *pipe_u, Damper *d, LWP *pipe_d);

    void Update(double t_target, LWP *pipe_u, Damper *damper, Valve *valve, double p_d, double T_d);

    //private:
    string name;
    bool DEBUG;
    PSToolboxBaseEdge *e1;
    PSToolboxBaseEdge *e2;
    PSToolboxBaseEdge *e3;
    PSToolboxBaseEdge *e4;
    bool is_front1, is_front2, is_front3, is_front4;
    string BC_type;
    double BC_value;
    double demand; // kg/s
    int type;

    vector<PSToolboxBaseEdge *> edges;
    vector<bool> is_front;
    vector<int> edges_idx;

    static void Set_LWP_BC(LWP *p1, const double t_target, const bool is_front,
                           const double p, const double T, const double v);

    double Connector_LWP_Reservoir_and_Pipe_Front_fun(double rho, double beta,
                                                      Reservoir *r1, LWP *p1, bool inlet_pressure_drop);

    double signed_sqrt(double x);

    bool connected_to_rigid();

    void setPressureDependentDemand(double k, double p0, double t0);

    double p_prev;
    double demand_k{};
    double demand_p0{};
    double demand_t0{};
    bool pressure_dependent_demand = false;
};
