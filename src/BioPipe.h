#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "PSToolboxBaseEdge.h"
#include "Units.h"

using namespace std;
using namespace Eigen;

class BioPipe : public PSToolboxBaseEdge, Units {
private:
    //Data
    double L, D, A, g, v_conv, dx;
    string bio_type;
    int Npts_min, Npts_max;
    double kb = 0.5 / 24. / 3600.;
    VectorXd x;
    VectorXd v_conv_vec;
    MatrixXd Cnew;
    vector<string> name_of_bio_vars;
    int num_of_bio_vars;
    bool is_v_conv_set;

    bool ini_done; //!< whether the pipe has been initialised or not
    void BCLeft(string type, double val, double &pstart, double &vstart);

    void BCRight(string type, double val, double &pend, double &vend);

    //double Source(int i);
    MatrixXd Source();


    //Saving and plotting data
    vector<vector<double> > data;
    vector<double> tmpvec;
    vector<double> tstatus; //time data for the status
    vector<vector<double> > vstatus, pstatus;
    bool save_data; //!< whether to save data or not
    //bool do_plot_runtime;
    //bool ylim_isset;
    //double lim_pbottom, lim_pup, lim_vbottom, lim_vup;

public:
    BioPipe(string _nev,
            string _cspe_nev,
            string _cspv_nev,
            double _L,
            double _D,
            string _bio_type,
            bool save_data);

    MatrixXd C;

    int Npts;

    ~BioPipe();

    string Info() override;

    void Ini() override;

    void Ini(double dt_target) override;

    void Ini(double dt_target, VectorXd Cini) override;

    double Get_dprop(string prop_string) override;

    int Get_iprop(string) override;

    vector<double> Get_dvprop(string prop_string);

    void Set_dprop(string prop_string, double val);

    void Step(
        string BC_start_type, double BC_start_val,
        string BC_end_type, double BC_end_val);

    void UpdateInternal(double) override;

    void UpdateTime(double) override;

    void Save_data() override;

    void Write_data(string folder) override;

    string fname;

    void Set_string_prop(string, string) override;

    VectorXd Get_C_Front() override {
        return C.col(0);
    }

    VectorXd Get_C_Back() override {
        return C.col(Npts - 1);
    }

    double Get_mp() {
        return (v_conv * 1000. * A);
    }

    void Set_inletQualityFront(VectorXd Cin) override {
        C.col(0) = Cin;
    }

    void Set_inletQualityBack(VectorXd Cin) override {
        C.col(Npts - 1) = Cin;
    }

    void Set_v_conv(double v) override {
        v_conv = v;
        is_v_conv_set = true;
        if (fabs(v_conv) < 0.001)
            v_conv = 0.001;
    }
};
