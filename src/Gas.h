#pragma once
class Gas
{
protected:
    double p;
    double rho;
    double T;
    double e;

public:
    Gas();
    virtual ~Gas();

    virtual double Get_p(double rho, double T) = 0;
    virtual double Get_rho(double p, double T) = 0;
    virtual double Get_T(double p, double rho) = 0;
    virtual double Get_e_from_Tp(double T,double p) = 0;
    virtual double Get_SonicVel(double T, double p) = 0;
    virtual double Get_T_from_ep(double e, double p) = 0;
    virtual double Get_T_from_erho(double e, double rho) = 0;
		virtual double Get_Prandtl_from_Tp(double T, double p)=0;
		virtual double Get_ThermalConductivity_from_Tp(double T, double p)=0;
		virtual double Get_DynamicViscosity_from_Tp(double T, double p)=0;
    virtual double Get_kappa_pv() = 0;
    virtual double Get_kappa_Tv() = 0;
    virtual double Get_R() = 0;
    virtual double Get_kappa_Tp() = 0;
    virtual double Get_cp(double p, double T) = 0;
    virtual double Get_cV(double p, double T) = 0;
    virtual double Get_eta_crit() = 0;
    virtual double Get_MassFlux(double pu, double rhou, double pd, double rhod) = 0;
};
