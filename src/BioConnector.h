#pragma once

#include "Reservoir.h"
#include "LWP.h"
#include "SCP.h"
#include "Valve.h"
#include "ThrottleValve.h"
#include "Damper.h"
#include <cstdlib>
#include "BioPipe.h"
//#include "PSToolboxBaseEdge.h"

class BioConnector {
public:
    // BioConnector(bool DEBUG);

    //type=0
    BioConnector(
        string name,
        vector<PSToolboxBaseEdge *> &_e,
        vector<bool> &_is_front,
        vector<int> &_edges_idx,
        double _demand, bool _DEBUG);

    // type=1
    BioConnector(
        string name,
        PSToolboxBaseEdge *e1, bool is_front,
        string BC_type, double BC_value,
        double demand, bool DEBUG);

    ~BioConnector();

    void Update(double t_target, int update_idx);

    void BioConnector_N_BioPipes(double t_target, int update_idx);

    void BioConnector_Reservoir(double t_target);

    void BioConnector_FreeEnd(double t_target);

    string name;
    bool DEBUG;
    VectorXd Cin; // If type=1 (free end) or 2 (reservoir), set the ingoing quality
    string BC_type;
    double BC_value;
    double demand; // kg/s
    int type;

    PSToolboxBaseEdge *s_edges; // s_ = single
    bool s_is_front;
    int s_edges_idx;

    vector<PSToolboxBaseEdge *> edges;
    vector<bool> is_front;
    vector<int> edges_idx;

    double signed_sqrt(double x);

    bool connected_to_rigid();

    void Set_InletQuality(VectorXd C_inlet) {
        Cin = C_inlet;
    }
};
