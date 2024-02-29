#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "PSToolboxBaseEdge.h"
#include "EpanetReader.h"
#include "Connector.h"
#include "SCP.h"

using namespace std;


double getPN10D150Delta(double D)
{
    const int N = 9;
    double Dn[N] = {90., 110., 140., 180., 200., 225., 280., 315., 400.};
    double Delta[N] = {5.4e-3, 6.6e-3, 8.3e-3, 10.7e-3, 11.9e-3, 13.4e-3, 16.6e-3, 18.7e-3, 23.7e-3};

    int idx = -1;
    for(int i = 0; i < N; i++)
    {
        if(D < Dn[i])
        {
            idx = i;
            break;
        }
    }
    if(idx == 0) return Delta[0]; //below first
    if(idx == -1) return Delta[N-1]; //above last

    return (Delta[idx-1]*(Dn[idx]-D) + Delta[idx]*(D - Dn[idx-1])) / (Dn[idx] - Dn[idx-1]);
}

double getPN10D101Delta(double D)
{
    const int N = 6;
    double Dn[N] = {80., 100., 150., 300., 400., 600.};
    double Delta[N] = {0.5e-3*(96.-80.), 0.5e-3*(122.-100.),0.5e-3*(177.3-150.),0.5e-3*(345.4-300.),0.5e-3*(453.1-400.),0.5e-3*(667.-600.)};

    int idx = -1;
    for(int i = 0; i < N; i++)
    {
        if(D < Dn[i])
        {
            idx = i;
            break;
        }
    }
    if(idx == 0) return Delta[0]; //below first
    if(idx == -1) return Delta[N-1]; //above last

    return (Delta[idx-1]*(Dn[idx]-D) + Delta[idx]*(D - Dn[idx-1])) / (Dn[idx] - Dn[idx-1]);
}

void calculatePropagationVelocity(EpanetReader & reader)
{
    //go through all the pipes
    for (int i = 0; i < reader.pipes.size(); i++)
    {
        double r = reader.pipes[i].Roughness;
        double D = reader.pipes[i].Diameter;
        double Ef = 2.18e9; //El. mod. of water
        double Ec, Delta = -1;

        if(r == 100) //KPE PN10 pipes
        {
            Ec = 300.0e6; //El. mod. of pipe
            if(D == 150) Delta = 9.5e-3;
            if(D == 173) Delta = 10.7e-3;
        }
        if(r == 101) //ACPN10 pipes
        {
            Ec = 19613.3e6;
            Delta = getPN10D101Delta(D);
        }
        if(r == 150 && (reader.pipes[i].ID == "23" || reader.pipes[i].ID == "24")) //KPE PN16 pipes
        {
            Ec = 300.0e6; //El. mod. of pipe
            Delta = 50.8e-3;
        }
        if(r == 150) //KPE PN16 pipes
        {
            Ec = 1400.0e6; //El. mod. of pipe
            Delta = getPN10D150Delta(D);
        }

        if(Delta == -1)
        {
            cout << "Pipe data not found: " <<  reader.pipes[i].ID <<endl;
            while(true) 1;
        }

        double Er = 1.0/(1.0/Ef + 0.001*D/(Delta * Ec));

        reader.pipes[i].SpeedOfSound = sqrt(Er/1000.0);
        reader.pipes[i].Delta = Delta;
        
        if(Delta == -1)
        {
            cout << "Pipe " << i << " ID: " << reader.pipes[i].ID << endl;
            cout << "\tRoughness= " << reader.pipes[i].Roughness << endl;
            cout << "\t       D = " << reader.pipes[i].Diameter << endl;
            cout << "\t   Delta = " << reader.pipes[i].Delta << endl;
            cout << "\t       a = " << reader.pipes[i].SpeedOfSound << endl << endl;
        }
    }
    
}


void calculateLambda(EpanetReader & reader)
{
    for (int i = 0; i < reader.pipes.size(); i++)
    {
        double v = 1.0;
        double D = reader.pipes[i].Diameter*0.001;
        double nu = 1e-6;
        
        double Re = v*D / nu;

        double lambda;
        if(Re < 2300.0)
        {
            lambda = 64.0 / Re;
        }
        else
        {
            lambda = 0.316/pow(Re,0.25);
        }


        //cout << "i = " << i << " D = "<< D << " Re = " << Re << endl; 

        reader.pipes[i].Lambda = lambda;
    }
    
}