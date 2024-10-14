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
    double Q, p1,p2;
    double g;
    double n_act, n_act_new;
    bool DEBUG;
bool is_rigid_element;
double dt_out;
	vector<string> transient_type;
	vector<double> transient_time;
	vector<double> transient_val;
	vector<bool> transient_triggered;


    PSToolboxBaseEdge(const string edge_type, const string name,  const string n1, const string n2);
    virtual ~PSToolboxBaseEdge(){};

    string Get_name(){return name;};
    string Get_type(){return edge_type;};

    void Set_t(double tnew){t=tnew;};
    double Get_t(){return t;};
    double Get_dt(){return dt;};
    double Get_tnext(){return t+dt;};
    void Set_DEBUG(bool newval){DEBUG=newval;};
    void Set_Q(double newval){
		Q=newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_Q -> Q set to "<<Q*3600<<" m3/h";
	};
	void Set_p1(double newval){
		p1=newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_p1 -> p1 set to "<<p1/1000/9.81<<" mwc";
	};
	void Set_p2(double newval){
		p2=newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_p2 -> p2 set to "<<p2/1000/9.81<<" mwc";
	};
	double Get_Q(){return Q;};
	double Get_p1(){return p1;};
	double Get_p2(){return p2;};


    virtual void Ini(int)=0;
    virtual void Ini()=0;
     virtual void Ini(double)=0;
    virtual string Info()=0;
    virtual void Set_BC_Left(string type, double val)=0;
    virtual void Set_BC_Right(string type, double val)=0;
    virtual void GetLargePressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&)=0;
    virtual void GetSmallPressureValues(double, vector<double>&, vector<double>&,vector<double>&, vector<string>&)=0;
    virtual void Save_data()=0;
    virtual void Write_data(string folder)=0;
    virtual void GetBetaAtFront(double t_target, double & LHS, double & coeff_Q)=0;
    virtual void GetAlphaAtEnd(double t_target, double & LHS, double &coeff_Q)=0;
    virtual void GetEdgeEquationCoeffs(double t_target, bool is_front, double & LHS, double & coeff_Q,double & coeff_p1,double & coeff_p2)=0;

    virtual double Get_dprop(string)=0;
    virtual void UpdateInternal(double t_target)=0;
    virtual void UpdateTime(double)=0;
    virtual void Set_string_prop(string,string)=0;
	virtual void Add_transient(string,double,double)=0;

   void Set_dt_out(double val){dt_out=val;};

private:

};
#endif
