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

void EpanetReader::readFromFile(const std::string &filename)
{
    ifstream file(filename);
    string line;
    string section;

    while (std::getline(file, line))
    {
        // cout << "len=" << line.length() << "\t";
        if (line.find("[JUNCTIONS]") != string::npos)
        {
            section = "JUNCTIONS";
        }
        else if (line.find("[RESERVOIRS]") != string::npos)
        {
            section = "RESERVOIRS";
        }
        else if (line.find("[TANKS]") != string::npos)
        {
            section = "TANKS";
        }
        else if (line.find("[PIPES]") != string::npos)
        {
            section = "PIPES";
        }
        else if (line.find("[PUMPS]") != string::npos)
        {
            section = "PUMPS";
        }
        else if (line.find("[VALVES]") != string::npos)
        {
            section = "VALVES";
            cout << "Valve found!\n";
        }
        else if (line.length() <= 1 || line[0] == ';')
        {
            cout << "Line skipped"
                    << "\n";
            continue; // Skip empty lines or comments
        }
        else
        {
            // cout << "Section = " << section<< "\n";
            if (section == "JUNCTIONS")
            {
                JunctionReader junction = parseJunction(line);
                cout << "Junction ID = " << junction.ID << "\telev = " << junction.Elev << "\tDemand = " << junction.Demand << "\t pattern = " << junction.Pattern << "\n";

                junctions.push_back(junction);
            }
            else if (section == "RESERVOIRS")
            {
                ReservoirReader reservoir = parseReservoir(line);
                cout << "Reservoir ID = " << reservoir.ID << "\thead = " << reservoir.Head << "\tPattern = " << reservoir.Pattern << "\n";

                reservoirs.push_back(reservoir);
            }
            else if (section == "TANKS")
            {
                TankReader tank = parseTank(line);
                cout << "Tank ID = " << tank.ID << "\tCurve = " << tank.VolCurve << "\t Overflow: " << tank.Overflow << "\n";

                tanks.push_back(tank);
            }
            else if (section == "PIPES")
            {
                PipeReader pipe = parsePipe(line);
                if (pipe.ID == "ClosedPipe - exclude")
                {
                    cout << "Closed PipeReader\n";
                }
                else
                {
                    cout << "Pipe ID = " << pipe.ID << "\tNode1 = " << pipe.Node1 << "\t Node2: " << pipe.Node2 << "\n";

                    pipes.push_back(pipe);
                }
            }
            else if (section == "PUMPS")
            {
                PumpReader pump = parsePump(line);
                cout << "Pump ID = " << pump.ID << "\tNode1 = " << pump.Node1 << "\t Node2: " << pump.Node2 << "\tPars: " << pump.Parameters << "\n";

                pumps.push_back(pump);
            }

            else if (section == "VALVES")
            {
                ValveReader valve = parseValve(line);
                cout << "Valve ID = " << valve.ID << "\tNode1 = " << valve.Node1 << "\t Node2: " << valve.Node2 << "\n";

                valves.push_back(valve);
            }
        }
    }
}

void EpanetReader::convertToRunner2()
{
    // allow the uniform treatment of junctionsans reservoirs
    unifyJunctions();
    cout << "Junctions unified!\n";

    for (int i = 0; i < pipes.size(); i++) // go through all the pipes
    {
        double a = pipes[i].SpeedOfSound; //m/s
        double L = pipes[i].Length; //m
        double D = 0.001 * pipes[i].Diameter; //m
        //double lambda = pipes[i].Lambda; //hs model
        double lambda = pipes[i].Roughness; //hw model
        int idx_e = findNodeByID(pipes[i].Node1);
        int idx_v = findNodeByID(pipes[i].Node2);
        double he = junctions[idx_e].Elev; //m
        double hv = junctions[idx_v].Elev; //m

        if(junctions[idx_e].type == 1) //reservoir at the start of pipe
        {
            he = hv;
        }
        if(junctions[idx_v].type == 1) //reservoir at the end of pipe
        {
            hv = he;
        }

        edges.push_back(new SCP(pipes[i].ID, pipes[i].Node1, pipes[i].Node2, 1000, a, L, D, lambda, he, hv, true));
        cout << "Pipe Call SCP (" <<  pipes[i].ID << ",\t" << pipes[i].Node1 << ",\t" << pipes[i].Node2 << ",\t1000," << a << ",\tL=" << L << ",\tD=" << D << ",\tlambda=" << lambda << ",\the=" << he << ",\thv=" << hv << ",false)" << endl; 
    }
    cout << "SCP pipes ready!\n";

    // go through all the pipes and find their respective nodes
    for (int i = 0; i < pipes.size(); i++)
    {
        // find the node index of pipe ends
        int idx1 = findNodeByID(pipes[i].Node1); // idx is start (1)
        int idx2 = findNodeByID(pipes[i].Node2); // idx is end (0)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(false);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(true);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);
    }
    cout << "Pipe end nodes found!\n";

    // create the connectivities
    for (int i = 0; i < junctions.size(); i++)
    {
        int size = junctions[i].idxPipe.size();
            cout << "i = " << i << "\tsize =" << size <<"\t";

        if(size == 0 || size > 4)
        {
            cout << "Undefined number of connections! Stopping code" << endl;
            while(true) 1;
        }
        if (size == 1)
        {
            // Connect pipe-BC
            int idx1 = junctions[i].idxPipe[0];
            bool end1 = junctions[i].end[0];
            string name = junctions[i].ID;

            // check type
            if (junctions[i].type == 0) //"free" end
            {
                double demand = 0.06*junctions[i].Demand / 3600.0;
                double D = pipes[idx1].Diameter*0.001;
                double A = 0.25*D*D*3.14159265;
                double v = demand / A; //????
                if(!end1) v = -v;
                cons.push_back(new Connector(name,edges[idx1], !end1, "Velocity", v, demand, true));
                cout << "Node " << name << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ",Velocity," << v <<"," << demand << ")\n";
            }
            if (junctions[i].type == 1) // reservoir
            {
                double head = junctions[i].Head;
                int idx_other = findNodeByID(getOtherNodeOfPipe(idx1,junctions[i].ID));
                double elev = junctions[idx_other].Elev;
                double p = (head-elev)*1000.0*9.81; //????
                double demand = 0;
                cons.push_back(new Connector(name,edges[idx1], !end1, "Pressure", p, demand, true));
                cout << "Node " << name << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ",Pressure,"<< p << "," << demand << ")\n";
            }
        }

        if (size == 2)
        {
            // Connect pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            double demand = 0.06*junctions[i].Demand / 3600.0;
            string name = junctions[i].ID;
            vector<int> idx_edge;
            idx_edge.push_back(idx1);
            idx_edge.push_back(idx2);
            cons.push_back(new Connector(name, edges[idx1], !end1, edges[idx2], !end2, demand, true, idx_edge));
            cout << "Node " << name << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << demand << ")\n";
        }
        if (size == 3)
        {
            // Connect pipe-pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            int idx3 = junctions[i].idxPipe[2];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            bool end3 = junctions[i].end[2];
            double demand = 0.06*junctions[i].Demand / 3600.0;
            string name = junctions[i].ID;
            vector<int> idx_edge;
            idx_edge.push_back(idx1);
            idx_edge.push_back(idx2);
            idx_edge.push_back(idx3);
            cons.push_back(new Connector(name,edges[idx1],!end1,edges[idx2],!end2,edges[idx3],!end3,demand,true,idx_edge));
            cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << pipes[idx3].ID << "," << end3  << "," << demand << ")\n";
        }
        if (size == 4)
        {
            // Connect pipe-pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            int idx3 = junctions[i].idxPipe[2];
            int idx4 = junctions[i].idxPipe[3];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            bool end3 = junctions[i].end[2];
            bool end4 = junctions[i].end[3];
            double demand = 0.06*junctions[i].Demand / 3600.0;
            string name = junctions[i].ID;
            vector<int> idx_edge;
            idx_edge.push_back(idx1);
            idx_edge.push_back(idx2);
            idx_edge.push_back(idx3);
            idx_edge.push_back(idx4);
            cons.push_back(new Connector(name,edges[idx1],!end1,edges[idx2],!end2,edges[idx3],!end3,edges[idx4],!end4,demand,false,idx_edge));
            cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << pipes[idx3].ID << "," << end3  << "," << pipes[idx4].ID << "," << end4  << "," << demand << ")\n";
        }
    }
}

void EpanetReader::PrintData()
{
    cout << endl << endl << "------------------CPP FILE--------------------" << endl << endl;
    cout << "Size of edges:\t" << edges.size() << endl;
    for (int i = 0; i < edges.size(); i++)
    {
        cout << "i = " << i << "\t" << edges[i]->name << endl;
    }


    cout << " Size of cons:\t" << cons.size() << endl;
    for (int i = 0; i < cons.size(); i++)
    {
        cout << "i = " << i << "\t BC: " << cons[i]->BC_type << "\t size of edgeidx: "<< cons[i]->edges_idx.size() << endl;
    }

    cout << "Size of start:\t" << con_at_edge_start.size() << endl;
    for (int i = 0; i < con_at_edge_start.size(); i++)
    {
        cout << "i = " << i << "\t" << con_at_edge_start[i] << endl;
    }

    cout << "  Size of end:\t" << con_at_edge_end.size() << endl;
    for (int i = 0; i < con_at_edge_end.size(); i++)
    {
        cout << "i = " << i << "\t" << con_at_edge_end[i] << endl;
    }
}

int EpanetReader::nextPipeAtNode(int idx)
{
    if (junctions[idx].idxPipe[0] == -1)
        return 0;
    if (junctions[idx].idxPipe[1] == -1)
        return 1;
    if (junctions[idx].idxPipe[2] == -1)
        return 2;
    return -1;
}

string EpanetReader::getOtherNodeOfPipe(int idxPipe, string Node)
{
    if (pipes[idxPipe].Node1 == Node)
        return pipes[idxPipe].Node2;
    else
        return pipes[idxPipe].Node1;
}

void EpanetReader::printEdgesAndCons()
{
    for (int i = 0; i < edges.size(); i++)
    {
        cout << edges[i]->Info();
    }
}

void EpanetReader::unifyJunctions()
{
    for (int i = 0; i < junctions.size(); i++)
    {
        junctions[i].type = 0;
    }

    for (int i = 0; i < reservoirs.size(); i++)
    {
        cout << "Reservoir " << reservoirs[i].ID << "\n";
        JunctionReader J;
        J.type = 1;
        J.Head = reservoirs[i].Head;
        J.ID = reservoirs[i].ID;
        junctions.push_back(J);
    }

    cout << "Njunctions = " << junctions.size() << "\n";
}

int EpanetReader::findNodeByID(const std::string ID)
{
    for (int i = 0; i < junctions.size(); i++)
    {
        if (junctions[i].ID == ID)
        {
            return i;
        }
    }
    return -1;
}

int EpanetReader::findReservoirByID(const std::string ID)
{
    for (int i = 0; i < reservoirs.size(); i++)
    {
        if (reservoirs[i].ID == ID)
        {
            return i;
        }
    }
    return -1;
}

bool EpanetReader::checkIfIncludedInVector(vector<int> list, int elem)
{
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i] == elem)
            return true;
    }
    return false;
}

vector<int> EpanetReader::findConnectingPipes(const std::string node)
{
    vector<int> connectingPipes;
    for (int i = 0; i < pipes.size(); i++)
    {
        if (node == pipes[i].Node1 || node == pipes[i].Node2)
        {
            // cout << "Found at " << i << "\n";
            connectingPipes.push_back(i);
        }
    }
    return connectingPipes;
}

JunctionReader EpanetReader::parseJunction(const std::string &line)
{
    istringstream iss(line);
    JunctionReader junction;
    iss >> junction.ID >> junction.Elev >> junction.Demand >> junction.Pattern;
    return junction;
}

ReservoirReader EpanetReader::parseReservoir(const std::string &line)
{
    istringstream iss(line);
    ReservoirReader reservoir;
    iss >> reservoir.ID >> reservoir.Head >> reservoir.Pattern;
    return reservoir;
}

TankReader EpanetReader::parseTank(const std::string &line)
{
    istringstream iss(line);
    TankReader tank;
    iss >> tank.ID >> tank.Elevation >> tank.InitLevel >> tank.MinLevel >> tank.MaxLevel >> tank.Diameter >> tank.MinVol >> tank.VolCurve >> tank.Overflow;
    return tank;
}

PipeReader EpanetReader::parsePipe(const std::string &line)
{
    istringstream iss(line);
    PipeReader pipe;
    string Status;
    iss >> pipe.ID >> pipe.Node1 >> pipe.Node2 >> pipe.Length >> pipe.Diameter >> pipe.Roughness >> pipe.MinorLoss >> Status;

    if (Status == "Open")
    {
        return pipe;
    }
    else
    {
        pipe.ID = "ClosedPipe - exclude";
        return pipe;
    }
}

PumpReader EpanetReader::parsePump(const std::string &line)
{
    istringstream iss(line);
    PumpReader pump;
    string tmp;
    iss >> pump.ID >> pump.Node1 >> pump.Node2 >> tmp >> pump.Parameters;
    return pump;
}

ValveReader EpanetReader::parseValve(const std::string &line)
{
    istringstream iss(line);
    ValveReader valve;
    string tmp;
    iss >> valve.ID >> valve.Node1 >> valve.Node2 >> valve.Diameter >> valve.Type >> valve.Setting >> valve.MinorLoss;
    return valve;
}

