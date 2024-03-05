#ifndef PSToolboxBaseEdge_H
#define PSToolboxBaseEdge_H
#include <iostream>
#include <vector>
using namespace std;

class PSToolboxBaseEdge
{
  protected:
    string name; 
    unsigned int type;
    bool save_data();
    double t=0., dt;
  public:
    PSToolboxBaseEdge(const string name);
    ~PSToolboxBaseEdge(){};

    string Get_name(){return name;};
    unsigned int Get_type(){return type;};

    double Get_t(){return t;};
    double Get_dt(){return dt;};
    double Get_tnext(){return t+dt;};
    virtual void Ini(int)=0;
    virtual void Ini()=0;

    virtual string Info()=0;

    virtual void Set_BC_Left(string type, double val)=0;
    virtual void Set_BC_Right(string type, double val)=0;
    virtual void Save_data()=0;
    virtual void Write_data()=0;
    //  virtual double GetBetaAtFront(double)=0;
    //virtual double GetAlphaAtEnd(double)=0;
    virtual void GetBetaAtFront(double, double &, double &)=0;
    virtual void GetAlphaAtEnd(double, double &, double &)=0;
    virtual double Get_dprop(string)=0;
    virtual void UpdateInternal()=0;
    virtual void UpdateTime(double)=0;
};
#endif
