#pragma once
#include "Gas.h"

class IdealGas : public Gas {
public:
    IdealGas(
        double kappa, ///< [in] Adiabatic exponent*/
        double R ///< [in] specific gas constant, J/kg/K */
    );

    explicit IdealGas(IdealGas *g);

    IdealGas();

    ~IdealGas() override;

    double Get_p(double rho, double T) override {
        return rho * R * T;
    }

    double Get_rho(double p, double T) override {
        return p / R / T;
    }

    double Get_T(double p, double rho) override {
        return p / rho / R;
    }

    double Get_e_from_Tp(double T, double p) override {
        return T * cV;
    }

    double Get_SonicVel(double T, double p) override {
        return sqrt(kappa * R * T);
    }

    double Get_T_from_ep(double e, double p) override {
        return e / cV;
    }

    double Get_T_from_erho(double e, double rho) override {
        return e / cV;
    }

    double Get_Prandtl_from_Tp(double T, double p) override;

    double Get_ThermalConductivity_from_Tp(double T, double p) override;

    double Get_DynamicViscosity_from_Tp(double T, double p) override;

    double Get_kappa_pv() override;

    double Get_kappa_Tv() override;

    double Get_kappa_Tp() override;

    double Get_R() override;

    double Get_pcrit();

    double Get_cp(double p, double T) override;

    double Get_cV(double p, double T) override;

    double Get_MassFlux(double pu, double Tu, double pd, double Td) override;

    double Get_eta_crit() override;

    double kappa;
    double R;
    double cp;
    double cV;
};
