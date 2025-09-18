#ifndef PSToolboxBaseEdge_H
#define PSToolboxBaseEdge_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class PSToolboxBaseEdge {
public:
	string name;
	string edge_type;
	string node_from, node_to;
	bool save_data;
	double t, dt;
	double Q, p1, p2;
	double g;
	double n_act, n_act_new;
	bool DEBUG;
	bool is_rigid_element;
	double dt_out;
	vector<string> transient_type;
	vector<double> transient_time;
	vector<double> transient_val;
	vector<bool> transient_triggered;

	bool suppress_all_output;

	MatrixXd C;

	PSToolboxBaseEdge(const string edge_type, const string name, const string n1, const string n2);

	virtual ~PSToolboxBaseEdge() {
	};

	string Get_name() { return name; };
	string Get_type() { return edge_type; };

	void Set_t(double tnew) { t = tnew; };
	double Get_t() { return t; };
	double Get_dt() { return dt; };
	double Get_tnext() { return t + dt; };
	void Set_DEBUG(bool newval) { DEBUG = newval; };

	void Set_Q(double newval) {
		Q = newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_Q -> Q set to "<<Q*3600<<" m3/h";
	};

	void Set_p1(double newval) {
		p1 = newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_p1 -> p1 set to "<<p1/1000/9.81<<" mwc";
	};

	void Set_p2(double newval) {
		p2 = newval;
		//cout<<endl<<"\t edge "<<name<<": PSToolboxBaseEdge::Set_p2 -> p2 set to "<<p2/1000/9.81<<" mwc";
	};
	double Get_Q() { return Q; };
	double Get_p1() { return p1; };
	double Get_p2() { return p2; };

	void set_supress_all_output(bool newval) { suppress_all_output = newval; };

	virtual void Ini(int) {
	};

	virtual void Ini() {
	};

	virtual void Ini(double) {
	};

	virtual void Ini(double vini, double pstart, double dt_target) {
	};

	virtual string Info() { return "null"; };

	virtual void Set_BC_Left(string type, double val) {
	};

	virtual void Set_BC_Right(string type, double val) {
	};

	virtual void GetLargePressureValues(double, vector<double> &, vector<double> &, vector<double> &,
	                                    vector<string> &) {
	};

	virtual void GetSmallPressureValues(double, vector<double> &, vector<double> &, vector<double> &,
	                                    vector<string> &) {
	};

	virtual void Save_data() =0;

	virtual void Write_data(string folder) =0;

	virtual void GetBetaAtFront(double t_target, double &LHS, double &coeff_Q) {
	};

	virtual void GetAlphaAtEnd(double t_target, double &LHS, double &coeff_Q) {
	};

	virtual void GetEdgeEquationCoeffs(double t_target, bool is_front, double &LHS, double &coeff_Q, double &coeff_p1,
	                                   double &coeff_p2) {
	};

	virtual double Get_dprop(string) { return 0.0; };

	virtual void Set_dprop(string, double) {
	};


	virtual int Get_iprop(string prop_string) {
		cout << endl << "Get_iprop() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
		return -1;
	}

	virtual void UpdateInternal(double t_target) =0;

	virtual void UpdateTime(double) =0;

	virtual void Set_string_prop(string, string) =0;

	virtual void Add_transient(string, double, double) {
	};

	void Set_dt_out(double val) { dt_out = val; };

	virtual MatrixXd Get_C() {
		cout << endl << "Get_C() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
		return MatrixXd();
	}


	virtual VectorXd Get_InnerMeanC() {
		cout << endl << "Get_InnerMeanC() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
		return VectorXd();
	}

	virtual void Set_inletQualityFront(VectorXd Cin) {
		cout << endl << "Set_inletQualityFront() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
	}

	virtual void Set_inletQualityBack(VectorXd Cin) {
		cout << endl << "Set_inletQualityBack() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
	}

	virtual void Set_max_vals(VectorXd _max_vals) {
		cout << endl << "Set_max_vals() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
	}


	virtual VectorXd Get_C_Front() {
		cout << endl << "Get_C_Front() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
		return (Eigen::VectorXd(1) << 1.0).finished();
	}

	virtual VectorXd Get_C_Back() {
		cout << endl << "Get_C_Back() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
		return (Eigen::VectorXd(1) << 1.0).finished();
	}

	virtual void Set_v_conv(double v) {
		cout << endl << "Set_v_conv() called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
	}

	virtual void Ini(double dt_target, VectorXd Cini) {
		cout << endl << "Ini(double dt_target, VectorXd Cini) called from PSToolboxBaseEdge.";
		cout << endl << "This should not happen." << endl;
		cin.get();
	}

private:
};
#endif
