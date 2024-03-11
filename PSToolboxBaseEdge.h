#ifndef PSToolboxBaseEdge_H
#define PSToolboxBaseEdge_H

#include <iostream>
#include <vector>

using namespace std;

class PSToolboxBaseEdge
{
public:
    string name; 
    string edge_type;
    string node_from, node_to;
    bool save_data;
    double t, dt;
    double Q;
    double g;
    bool DEBUG;

    PSToolboxBaseEdge(const string edge_type, const string name,  const string n1, const string n2);
    virtual ~PSToolboxBaseEdge(){};

    string Get_name(){return name;};
    string Get_type(){return edge_type;};

    double Get_t(){return t;};
    double Get_dt(){return dt;};
    double Get_tnext(){return t+dt;};
    void Set_DEBUG(bool newval){DEBUG=newval;};

    virtual void Ini(int)=0;
    virtual void Ini()=0;
    virtual string Info()=0;
    virtual void Set_BC_Left(string type, double val)=0;
    virtual void Set_BC_Right(string type, double val)=0;
    virtual void GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&)=0;
//<<<<<<< HEAD
//=======
    virtual void GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&)=0;
//>>>>>>> d971427eb223c4d6b1e9771c74598bc3624272c8
    virtual void Save_data()=0;
    virtual void Write_data()=0;
    virtual void GetBetaAtFront(double t_target, double & LHS, double & coeff_Q)=0;
    virtual void GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q)=0;
    virtual double Get_dprop(string)=0;
    virtual void UpdateInternal()=0;
    virtual void UpdateTime(double)=0;

//<<<<<<< HEAD

private:


//=======
//>>>>>>> d971427eb223c4d6b1e9771c74598bc3624272c8
};
#endif
