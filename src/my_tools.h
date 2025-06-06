#pragma once

#include <vector>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <algorithm>

using namespace Eigen;
using namespace std;

int getIndexOf(std::vector<int> v, int num);

double interpolate_val(vector<double> x, vector<double> y, double x0);

vector<double> interp1(vector<double> x, vector<double> y, vector<double> xi);

bool NR(
    VectorXd(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, const VectorXd & pars, const bool show_iter, const double err_max, const int iter_max,



const bool stop_on_failure
);

MatrixXd num_jac(
    VectorXd(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, const VectorXd & pars);

double signum(double x);

double max(vector<double> x);

double min(vector<double> x);

int max_idx(vector<double> x);

int min_idx(vector<double> x);

double mean(vector<double> x);

double stdev(vector<double> x);

void SteepestDescent(
    double(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, const VectorXd & pars, const bool show_iter);

VectorXd num_grad(
    double(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, const VectorXd & pars);

void NelderMead2D(
    double(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, double radius, const VectorXd & pars, const bool show_iter);

void NelderMead(
    double(*fun)(VectorXd & x, const VectorXd & pars),
    VectorXd & x, double radius, const VectorXd & pars, const double x_TOL, const double fun_TOL, const bool show_iter);

void swap(VectorXd x, VectorXd y);

double gen_random_double(double min, double max);

VectorXd linspace_VectorXd(double x_min, double x_max, int Nsteps, bool is_linear);

vector<double> linspace(double x_min, double x_max, int Nsteps, bool is_linear);

double signum(double x);

//Add methods for easier reading of data

struct Datapair {
    std::string name;
    double value;
};

class Input {
public:
    void add_datapair(Datapair _dp);

    double extract_value(std::string name);

    void parse_file(std::string filename);

private:
    Datapair dataarray[100];
    //int elements{ 0 };
    int elements;
};
