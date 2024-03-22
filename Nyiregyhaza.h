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
    double Dn[N] = {79., 98., 141., 176., 200., 220., 278., 300., 450};
    double Delta[N] = {5.4e-3, 6.6e-3, 8.3e-3, 12.0e-3, 13.4e-3, 14.8e-3, 18.7e-3, 21.1e-3, 26.5e-3};

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
    //double Delta[N] = {0.5e-3*(96.-80.), 0.5e-3*(122.-100.),0.5e-3*(177.3-150.),0.5e-3*(345.4-300.),0.5e-3*(453.1-400.),0.5e-3*(667.-600.)};
    double Delta[N] = {9.5e-3, 11.0e-3, 13.5e-3, 24.e-3, 32.e-3, 46.e-3};

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
            Ec = 900.0e6; //El. mod. of pipe
            if(D == 150) Delta = 10.6e-3;
            if(D == 173) Delta = 11.9e-3;
        }
        else if(r == 101) //ACPN10 pipes
        {
            Ec = 19613.3e6;
            Delta = getPN10D101Delta(D);
        }
        else if(r == 150 && (reader.pipes[i].ID == "23" || reader.pipes[i].ID == "24")) //KPE PN16 pipes
        {
            Ec = 900.0e6; //El. mod. of pipe
            Delta = 63.1e-3;
        }
        else if(r == 150) //KPE PN16 pipes
        {
            Ec = 3200.0e6; //El. mod. of pipe
            Delta = getPN10D150Delta(D);
        }
        else
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


void modifyDemand(EpanetReader & reader, string setting)
{
    for (int i = 0; i < reader.junctions.size(); i++)
    {
        if(reader.junctions[i].Pattern == "11")
        {   
            if(setting == "Max")
            {
                reader.junctions[i].Demand *= 0.0611;
            }
            else if(setting == "Min")
            {
                reader.junctions[i].Demand *= 0.0138;
            }
            else
            {
                cout << "Setting does not exist"<< endl;
                while(1) 1;
            }
        }
        else if(reader.junctions[i].Pattern == "Fix")
        {
            reader.junctions[i].Demand *= 1.0;
        }
        else if(reader.junctions[i].Pattern == "ipar")
        {
            reader.junctions[i].Demand *= 0.0416;
        }
        else
        {
            cout << "Node " << reader.junctions[i].ID << " unknown demand pattern: " << reader.junctions[i].Pattern << endl;
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