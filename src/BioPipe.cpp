#define _USE_MATH_DEFINES

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdio>
#include "my_tools.h"
#include "PSToolboxBaseEdge.h"
#include "BioPipe.h"
#include <cmath>

//! Constructor for the slightly compressible pipe class
BioPipe::BioPipe(const string _name, //!< [in] Name of the slightls compressible pipe
                 const string _cspe_name, //!< [in] Name of the previous node
                 const string _cspv_name, //!< [in] Name of the next node
                 const double _L, //!< [in] Length of the pipe section
                 const double _D,
                 const string _bio_type, //!< [in] Diameter of the pipe
                 const bool _save_data /**< [in] whether to save data*/
) : PSToolboxBaseEdge("BioPipe", _name, _cspe_name, _cspv_name), Units() {
  save_data = _save_data;
  L = _L;
  D = _D;
  rh = D / 4.;
  A = D * D * M_PI / 4.; //Cross sectional area of the pipe
  ini_done = false;
  fname = name + ".dat"; //name of the file for saving data about the pipe
  g = 9.81;

  bio_type = _bio_type;

  // Dummy initialization
  Npts = 1;
  dt = 0.;
  v_conv = 0.;

  is_rigid_element = false;

  Npts_min = 10;
  Npts_max = 200;
  is_v_conv_set = false;
}


//! Desctructor of the BioPipe class
BioPipe::~BioPipe() = default;

string BioPipe::Info() {
  bool show_pts = true;
  if (!ini_done) {
    cout << endl << "ERROR!!!";
    cout << endl << "BioPipe::Info() called on " << name << " but the element is uninitialized!" << endl;
    cin.get();
  }

  std::ostringstream oss;
  oss << "\n\n Pipe name     : " << name;
  oss << "\n================================";
  oss << "\n      node_from : " << node_from;
  oss << "\n        node_to : " << node_to;
  oss << "\n              L : " << L << " m = " << L * m_to_inch << " in";
  oss << "\n              D : " << D << " m= " << D * m_to_inch << " in";
  oss << "\n             dt : " << dt << " s";
  oss << "\n       bio_type : " << bio_type;
  oss << "\n           Npts : " << Npts;
  oss << "\n         v_conv : " << v_conv << " m/s" << endl;

  //oss << "\n        node #  : ";
  //  for (int i = 0; i < Npts; i++)
  //    oss << setw(4) << setfill(' ') << i << " ";

  oss << std::showpos << std::scientific << std::setprecision(2) << std::setw(4);
  for (int i = 0; i < C.rows(); i++) {
    oss << "\n " << name_of_bio_vars.at(i) << " :";
    for (int j = 0; j < C.cols(); j++)
      oss << C(i, j) << "  ";
  }
  oss << endl;

  return oss.str();
}

void BioPipe::Ini() {
  bio_type = "water_age";
  Ini(1.0, (VectorXd(1) << 1.0).finished());
}

void BioPipe::Ini(double dt_target) {
  Ini(1.0, (VectorXd(num_of_bio_vars) << 0.0).finished());
}


/*! \brief Initialize the BioPipe pipe
  Initializes the BioPipe pipe with ....
  */
void BioPipe::Ini(double dt_target, VectorXd Cini) {
  t = 0.;
  if (!is_v_conv_set) {
    cout << endl << "ERROR!!!";
    cout << endl << "BioPipe::Ini() called on " << name << " but v_conv is not set!" << endl;
    cin.get();
  }

  Npts = round(L / fabs(v_conv) / dt_target); // CFL condition reorganized
  printf("\n BioPipe %10s:  L=%6.1f m, D=%6.1f m, v_conv=%5.3f m/s, Tpipe=%5.3f s, dt_target=%5.3e s, Npts=%3d ",
         name.c_str(), L, D, v_conv, L / fabs(v_conv), dt_target, Npts);
  if (Npts < Npts_min) {
    Npts = Npts_min;
    printf(" -> %d", Npts);
  }
  if (Npts > Npts_max) {
    Npts = Npts_max;
    printf(" -> %d", Npts);
  }

  // CFL criteria for selecting dt
  dx = L / (Npts - 1);
  dt_conv = dx / fabs(v_conv);
  dt = dt_conv;

  if (dt > dt_target)
    dt = dt_target;

  if (bio_type == "tracer") {
    num_of_bio_vars = 1;
    name_of_bio_vars.push_back("    tracer (-)");
    v_conv_vec = VectorXd::Zero(num_of_bio_vars);
    v_conv_vec(0) = v_conv;
  } else {
    if (bio_type == "water_age") {
      num_of_bio_vars = 1;
      name_of_bio_vars.push_back("water_age (hours)");
      v_conv_vec = VectorXd::Zero(num_of_bio_vars);
      v_conv_vec(0) = v_conv;
    } else {
      if (bio_type == "chlorine") {
        num_of_bio_vars = 1;
        name_of_bio_vars.push_back("C (mg/l)");
        v_conv_vec = VectorXd::Zero(num_of_bio_vars);
        v_conv_vec(0) = v_conv;
      } else {
        if (bio_type == "microbiology") {
          num_of_bio_vars = 4;
          name_of_bio_vars.push_back("  S_b (mg/l)  ");
          name_of_bio_vars.push_back("  X_b (mg/l) ");
          name_of_bio_vars.push_back("  S_w (mg/l)  ");
          name_of_bio_vars.push_back("  X_w (mg/m2) ");
          v_conv_vec = VectorXd::Zero(num_of_bio_vars);
          v_conv_vec(0) = v_conv;
          v_conv_vec(1) = v_conv;
          v_conv_vec(2) = 0.0;
          v_conv_vec(3) = 0.0;
        } else {
          cout << endl << "ERROR! BioPipe::BioPipe()";
          cout << endl << "\tUnknown bio_type: " << bio_type;
          cout << endl << "\tPossible choices: tracer, water_age, chlorine, microbiology" << endl;
        }
      }
    }
  }
  x = VectorXd::Zero(Npts);
  for (int i = 0; i < Npts; i++)
    x(i) = i * L / (Npts - 1);

  C = MatrixXd::Zero(num_of_bio_vars, Npts);
  for (int j = 0; j < Npts; ++j)
    C.col(j) = Cini;
  Cnew = C;

  ini_done = true;

  if (save_data) {
    tmpvec.push_back(t);
    for (int i = 0; i < num_of_bio_vars; i++)
      tmpvec.push_back(C(i, 0));
    for (int i = 0; i < num_of_bio_vars; i++)
      tmpvec.push_back(C(i, Npts - 1));
    data.clear();
    data.reserve(100);
    data.push_back(tmpvec);
  }
}

int BioPipe::Get_iprop(string prop_string) {
  int out;
  if (prop_string == "Npts")
    out = Npts;
  else if (prop_string == "num_of_bio_vars")
    out = num_of_bio_vars;
  else {
    cout << endl
        << "ERROR! BioPipe::Get_iprop(prop_string), unknown input: prop_string=" << prop_string << endl
        << endl;
    cout << endl << "Name of pipe: " << name << endl;
    exit(-1);
  }
  return out;
}

/*! \brief Return a property of the BioPipe.
  \param prop_string: A string describing the needed property ("L" | "L_feet" | "D" | "D_inch" | "A" | "dt" | "p_front" | "p_back" | "v_front" | "v_back" | "mp_front" | "mp_back" | "frek" | "tnext" | "lambda" | "phi" | "mu" | "rho" | "a")
  \return The value that was looked up.
  */
double BioPipe::Get_dprop(string prop_string) {
  double out = 0.0;
  if (prop_string == "L")
    out = L;
  else if (prop_string == "L_feet")
    out = L * m_to_ft;
  else if (prop_string == "D")
    out = D;
  else if (prop_string == "D_inch")
    out = D * m_to_inch;
  else if (prop_string == "A")
    out = A;
  else if (prop_string == "dt")
    out = dt;
  else if (prop_string == "tnext")
    out = t + dt;
  else if (prop_string == "v_conv")
    out = v_conv;
  else if (prop_string == "mp_back")
    out = v_conv * A * 1000.;
  else if (prop_string == "mp_front")
    out = v_conv * A * 1000.;
  else if (prop_string == "chlorine_kb")
    out = kb;
  else if (prop_string == "bio_mumax")
    out = mumax;
  else if (prop_string == "bio_Y")
    out = Y;
  else if (prop_string == "bio_a")
    out = a;
  else if (prop_string == "bio_kmort")
    out = kmort;
  else if (prop_string == "bio_kb")
    out = kb_bio;
  else if (prop_string == "bio_kfs")
    out = kfs;
  else {
    cout << endl
        << "ERROR! BioPipe::Get_dprop(prop_string), unknown input: prop_string=" << prop_string << endl
        << endl;
    cout << endl << "Name of pipe: " << name << endl;
    exit(-1);
  }
  return out;
}

/*! \brief Set a parameter
  Allows to reset the diameter or length the pipe. Updating the length will rescale the system.
  \param prop_string [in] The paramtere to change. Accepted values are "L" and "D"
  \param val [in] The new value of the set parameter, in SI units.
  \sa Get_dprop
  */
void BioPipe::Set_dprop(string prop_string, double val) {
  if (prop_string == "D") {
    D = val;
  } else if (prop_string == "L") {
    L = val;
    dt = L / (Npts - 1) / v_conv;
  } else if (prop_string == "chlorine_kb") {
    kb = val;
  } else if (prop_string == "bio_mumax")
    mumax = val;
  else if (prop_string == "bio_Y")
    Y = val;
  else if (prop_string == "bio_a")
    a = val;
  else if (prop_string == "bio_kmort")
    kmort = val;
  else if (prop_string == "bio_kb")
    kb_bio = val;
  else if (prop_string == "bio_kfs")
    kfs = val;
  else {
    cout << endl
        << "HIBA! BioPipe::Set_dprop(prop_string), ismeretlen bemenet: prop_string=" << prop_string << endl
        << endl;
  }
}

void BioPipe::Set_string_prop(string prop_string, string val) {
  if (prop_string == "bio_type") {
    bio_type = val;
  } else {
    cout << endl
        << "ERROR! BioPipe::Set_string_prop, unknown input:" << prop_string;
    cout << endl << "Valid inputs: bio_type (water_age, chlorine, microbiology)" << endl;
  }
}


void BioPipe::UpdateInternal(double t_target) {
  DEBUG = false;
  double ERR_MIN;
  double ERR_MAX;
  if (bio_type == "water_age") {
    ERR_MIN = 1.;
    ERR_MAX = 10.;
  } else {
    ERR_MIN = 1e-6;
    ERR_MAX = 1e-5;
  }
  double dt_ = t_target - t;
  double t_local = t;
  //double dt_ = dt;

  const int MAX_RETRIES = 25;
  const double MIN_DT = 1e-14; // safety floor to avoid underflow
  double err;

  //if (name == "P-45")
  // printf("\n %s", name.c_str());
  // printf("\n\t  t=%5.3e (pipe internal) | t_local=%5.3e, dt_=%5.3e, t_local+dt=%5.3e, t_target=%5.3e ***", t, t_local,
  //        dt_, t_local + dt_, t_target);
  if (dt_ < 1.e-10) {
    printf("\n %s: extremely small time step detected: %g", name.c_str(), dt_);
    cin.get();
  }

  while (t_local < t_target) {
    int retries = 0;
    while (true) {
      // This leads to problems, might give very small (1e-15) time step values, which remain there for the rest of the simulations.
      // Do not use.
      //if (t_local + dt_ > t_target) {
      //  dt_ = t_target - t_local;
      //      if (name == "P-45")
      //  printf("\n\t                   last step: | t_local=%5.3e, dt_=%5.3e, t_local+dt=%5.3e, t_target=%5.3e",
      //         t_local, dt_, t_local + dt_, t_target);
      //}

      // Slopes
      MatrixXd k1 = RHS(C);
      MatrixXd Y_Eu = C + dt_ * k1; // Euler predictor (not strictly needed for Heun update)
      MatrixXd k2 = RHS(Y_Eu);

      // Heun update (candidate)
      MatrixXd Ccandidate = C + dt_ * (k1 + k2) / 2.0;

      // Error indicator (curvature of RHS over the step)
      err = (k2 - k1).norm();

      // Accept or reject
      if (err <= ERR_MAX) {
        // Accept the step
        Cnew = std::move(Ccandidate);
        dt = dt_; // keep the (possibly reduced) step size
        if (dt > dt_conv)
          dt = dt_conv;
        break;
      }
      // if (name == "P_45")
      // printf("\n\t rejecting step | t_local=%5.3e, dt_=%5.3e, t_local+dt=%5.3e, t_target=%5.3e, err=%5.3e", t,
      //        t_local, dt_,
      //        t_local + dt_, t_target, err);

      // Reject the step: shrink dt and retry
      dt_ *= 0.5;
      ++retries;

      // Safety stops to avoid infinite loop or absurdly small dt
      if (retries >= MAX_RETRIES || dt_ < MIN_DT) {
        // Force accept to prevent stalling; you can also throw if preferred.
        Cnew = std::move(Ccandidate);
        dt = dt_;
        cout << endl << "Pipe " << name << " integrator reached MAX_RETRIES (" << MAX_RETRIES << ") or MIN_DT (" <<
            MIN_DT
            << "), ";
        cout << "but error is still " << err << " > ERR_MAX=" << ERR_MAX <<
            ". Accepting step but this might cause a numerical problem. " << endl;
        break;
      }
    }
    // if (name == "P-45")
    // printf("\n\t  t=%5.3e (pipe internal) | t_local=%5.3e, dt_=%5.3e, t_local+dt=%5.3e, t_target=%5.3e, err=%5.3e", t,
    //         t_local,
    //         dt_, t_local + dt_, t_target, err);

    t_local += dt_;
    // Update values
    double min_vals = 1.e-8;
    for (int i = 0; i < Npts; i++)
      for (int j = 0; j < num_of_bio_vars; j++) {
        if (Cnew(j, i) < min_vals)
          C(j, i) = min_vals;
        else {
          if (Cnew(j, i) > max_vals(j))
            C(j, i) = max_vals(j);
          else
            C(j, i) = Cnew(j, i);
        }
      }
  }
  // if (name == "P-45") {
  // cout << endl << "\t closing time step.";
  // cin.get();
  // }

  if (err < ERR_MIN) {
    dt = 1.2 * dt_;
    if (DEBUG)
      cout << endl << "Pipe " << name << ": err=" << err << " -> dt increased: " << dt << endl;
  }
  t = t_target;
}

MatrixXd BioPipe::RHS(MatrixXd C) {
  MatrixXd S = Source();
  MatrixXd F = MatrixXd::Zero(num_of_bio_vars, Npts);

  if (bio_type == "microbiology") {
    if (v_conv > 0.) {
      for (int i = 1; i < Npts; i++)
        for (int j = 0; j < num_of_bio_vars; j++)
          F(j, i) = -v_conv_vec(j) * (C(j, i) - C(j, i - 1)) / dx + S(j, i);
      //C_new(j, i) = 0.*C(j, i) + (-v_conv_vec(j) * (C(j, i) - C(j, i - 1)) / dx + S(j, i)) * dt;
    } else {
      for (int i = 0; i < Npts - 1; i++)
        for (int j = 0; j < num_of_bio_vars; j++)
          F(j, i) = -v_conv_vec(j) * (C(j, i + 1) - C(j, i)) / dx + S(j, i);
      //C_new(j, i) = C(j, i) + (-v_conv_vec(j) * (C(j, i + 1) - C(j, i)) / dx + S(j, i)) * dt;
    }
  } else {
    if (v_conv_vec(0) > 0.)
      for (int i = 1; i < Npts; i++)
        F.col(i) = -v_conv_vec(0) * (C.col(i) - C.col(i - 1)) / dx + S.col(i);
    else
      for (int i = 0; i < Npts - 1; i++)
        F.col(i) = -v_conv_vec(0) * (C.col(i + 1) - C.col(i)) / dx + S.col(i);
  }


  return F;
}

void BioPipe::UpdateTime(double _t) {
  t = _t;

  /*
    if (name == "5"){
      cout << endl<<endl<<" Edge 5: "<<C;
      cin.get();
      }
    if (name == "9"){
      cout << endl<<endl<<" Edge 9 : "<<C;
  cin.get();}

   */
}

void BioPipe::Step(
  string BC_start_type, double BC_start_val,
  string BC_end_type, double BC_end_val) {
  double dummy_t_target = 0.;

  cout << endl << endl << "??????????? Should not be using this function!!!" << endl;
  cin.get();

  UpdateInternal(dummy_t_target);

  //BCLeft(BC_start_type, BC_start_val, pnew(0), vnew(0));

  //BCRight(BC_end_type, BC_end_val, pnew(Npts - 1), vnew(Npts - 1));

  //	t += dt;
}

void BioPipe::Save_data() {
  tmpvec.at(0) = t;
  for (int i = 0; i < num_of_bio_vars; i++)
    tmpvec.at(tmpvec.size() + i) = C(i, 0);
  for (int i = 0; i < num_of_bio_vars; i++)
    tmpvec.at(tmpvec.size() + i) = C(i, Npts - 1);

  data.push_back(tmpvec);
}

MatrixXd BioPipe::Source() {
  MatrixXd S = MatrixXd::Zero(num_of_bio_vars, Npts);

  if (bio_type == "tracer") {
    for (unsigned int i = 0; i < Npts; i++)
      for (unsigned int j = 0; j < num_of_bio_vars; j++)
        S(j, i) = 0.0;
  } else {
    if (bio_type == "water_age") {
      for (unsigned int i = 0; i < Npts; i++)
        for (unsigned int j = 0; j < num_of_bio_vars; j++)
          S(j, i) = 1.0;
    } else {
      if (bio_type == "chlorine") {
        for (unsigned int i = 0; i < Npts; i++)
          for (unsigned int j = 0; j < num_of_bio_vars; j++)
            S(j, i) = -kb * C(j, i);
      } else {
        if (bio_type == "microbiology") {
          for (unsigned int i = 0; i < Npts; i++) {
            double Sb = C(0, i);
            double Xb = C(1, i);
            double Sw = C(2, i);
            double Xw = C(3, i);

            // S_b (mg/l)
            // X_b (mg/l)
            // S_w (mg/l)
            // X_w (mg/m2)


            S(0, i) = -1 / Y * mumax * Sb / (kb_bio + Sb) * Xb + a * kmort * Xb - kfs / rh * (Sb - Sw);
            S(1, i) = mumax * Sb / (kb_bio + Sb) * Xb - kmort * Xb - kfs / rh * (Xb - Xw / rh) * 0.;
            S(2, i) = -1 / Y * mumax * Sw / (kb_bio + Sw) * Xw + a * kmort * Xw / rh + kfs / rh * (Sb - Sw);
            S(3, i) = mumax * Sw / (kb_bio + Sw) * Xw - kmort * Xw / rh + kfs / rh * (Xb - Xw / rh) * 0.;

            // if (name == "800014807") {
            //   if (i == 2) {
            //     cout << endl << endl << "name: " << name;
            //     printf("\nNpoints  = %d", Npts);
            //     printf("\ndt       = %g", dt);
            //     printf("\ndtconv   = %g", L / (Npts - 1) / fabs(v_conv));
            //     printf("\nmumax    = %5.3e", mumax);
            //     printf("\nkb_bio   = %5.3e", kb_bio);
            //     printf("\nSb       = %5.3e", Sb);
            //     printf("\nXb       = %5.3e", Xb);
            //     printf("\nSource b = %5.3e", S(1, i));
            //     printf("\nSw       = %5.3e", Sw);
            //     printf("\nXw       = %5.3e", Xw);
            //     printf("\nSource w = %5.3e", S(3, i));
            //     printf("\nkfs.     = %5.3e", kfs);
            //     printf("\nrh.      = %5.3e", rh);
            //     printf("\nt_diff   = %5.3e", rh / kfs / 3600.);
            //     cin.get();
            //   }
            // }
          }
        } else {
          cout << endl << "ERROR!BioPipe::Source()";
          cout << endl << "\tUnknown bio_type: " << bio_type;
          cout << endl << "\tPossible choices: water_age, chlorine, microbiology" << endl;
        }
      }
    }
  }

  return S;
}

/*! \brief Exports saved data
  Exports data saved from previous settings (such as steps and initialization.) Only works if the
  save_data flag was set to true, otherwise an error message is displayed. Data is saved to the previously determined
  savefile name. No inputs/outputs.
  */
void BioPipe::Write_data(string folder) {
  //char fname [50];
  //sprintf (fname, "%s.dat", name.c_str());

  if (!save_data) {
    cout << endl << "WARNING! BioPipe: " << name;
    cout << endl << " --> save_data = false was set in the constructor, cannot save anything..." << endl << endl;
  } else {
    cout << endl << "Saving to " << fname.c_str() << " ... ";

    FILE *pFile;
    pFile = fopen((folder + "/" + fname).c_str(), "w");
    fprintf(pFile, "t (s); p(0) (bar); p(L) (bar); v(0) m/s; v(L) (m/s); Q(0) (m3/h), Q(L) (m3/h), L, D, lambda\n");
    for (int i = 0; i < data.size(); i++)
      fprintf(pFile, "%8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e; %8.6e\n",
              data.at(i).at(0),
              data.at(i).at(1) / 1.e5, data.at(i).at(2) / 1.e5,
              data.at(i).at(3), data.at(i).at(4),
              data.at(i).at(5), data.at(i).at(6),
              L, D, D);
    fclose(pFile);
    cout << " done. ";
  }
}

vector<double> BioPipe::Get_dvprop(string prop_string) {
  int Ntime = data.size();
  //	int Nvars = data.at(0).size();
  //cout<<endl<<"Ntime="<<Ntime<<", Nvars="<<Nvars<<endl;
  vector<double> out(Ntime);
  if (prop_string == "t")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(0);
  else if (prop_string == "p_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1);
  else if (prop_string == "p_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(2);
  else if (prop_string == "p_front_bar")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(1) / 1.e5;
  else if (prop_string == "p_back_bar")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(2) / 1.e5;
  else if (prop_string == "v_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(3);
  else if (prop_string == "v_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(4);
  else if (prop_string == "mp_front")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(5);
  else if (prop_string == "mp_back")
    for (unsigned int i = 0; i < Ntime; i++)
      out.at(i) = data.at(i).at(6);
  else {
    cout << endl
        << "ERROR! BioPipe::Get_dvprop(prop_string), unknown input: prop_string=" << prop_string << endl
        << endl;
    cout << endl << "Name of valve: " << name << endl;
    cin.get();
  }
  return out;
}
