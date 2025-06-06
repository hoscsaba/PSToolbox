// g++ -I/usr/local/include/eigen3 -L/Users/hoscsaba/program/PSToolbox -lPSToolbox -lmy_tools -pedantic -O3 -Wall -Wno-c++11-long-long 2pipes_with_ideal_gas.cpp
#include <iostream>
#include <stdio.h>
#include "../src/LWP.h"
#include "../src/IdealGas.h"
#include "../src/Connector.h"

using namespace std;

int main(int argc, char **argv) {

    double p0 = 1.0e5, dp = 0.1*1.e5, T0 = 293;
    double L1 = 10., D1 = 0.1, lambda1 = 0.02;
    double L2 = 10., D2 = 0.02, lambda2 = 0.02;
    double L3 = 10., D3 = 0.05, lambda3 = 0.02;

    IdealGas* gas = new IdealGas(1.4, 287.);
    LWP *p1 = new LWP("p1", "csp1", "csp2", L1, D1, lambda1, 0, 0, true, false, gas);
    LWP *p2 = new LWP("p2", "csp2", "csp3", L2, D2, lambda2, 0, 0, true, false, gas);
    LWP *p3 = new LWP("p3", "csp3", "csp4", L3, D3, lambda3, 0, 0, true, false, gas);
    Connector c(true);

    p1->IniUniform(0., p0+dp, T0, 2);
    p2->IniUniform(0., p0, T0, 2);
    p3->IniUniform(0., p0, T0, 2);

    double pi, pL, pR, po;
    double vi, vL, vR, vo;
    double TL=293.;
    double t = 0., tmax = 2., dt = 1., dt2 = dt, dt3 = dt;

    while (t < tmax) {
        dt  = p1->Get_dt() * 0.95;
        dt2 = p2->Get_dt() * 0.95;
        dt3 = p3->Get_dt() * 0.95;
        if (dt > dt2)
            dt = dt2;
        if (dt > dt3)
            dt = dt3;
        t += dt;

        c.Connector_LWP_Pipes(t, p1, false, p2, true, p3, true);
        bool x1 = p1->BCLeft(t,"StaticPres_and_StaticTemp_Inlet",p0+dp,TL,true);
        bool x2 = p2->BCRight(t,"StaticPres_Outlet",p0,0.,true);
        bool x3 = p3->BCRight(t,"StaticPres_Outlet",p0,0.,true);
        p1->Step(dt);
        p2->Step(dt);
        p3->Step(dt);

        pi=p1->Get_dprop("p_front");
        pL=p1->Get_dprop("p_back");
        pR=p2->Get_dprop("p_front");
        po=p2->Get_dprop("p_back");

        vi=p1->Get_dprop("v_front");
        vL=p1->Get_dprop("v_back");
        vR=p2->Get_dprop("v_front");
        vo=p2->Get_dprop("v_back");

        printf("\n t=%5.3f s, p: %5.3f -> %5.3f -> %5.3f -> %5.3f bar",
                t, pi/1.e5, pL/1.e5, pR / 1.e5, po/1.e5);
    }

    p1->Save_data();
    p2->Save_data();
    p3->Save_data();

}
