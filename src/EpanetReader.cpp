#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "PSToolboxBaseEdge.h"
#include "EpanetReader.h"
#include "Connector.h"
#include "SCP.h"
#include "BioConnector.h"
#include "BioPipe.h"

using namespace std;


void EpanetReader::readFromFile(const std::string &filename, const string &unit) {
    // Convert to m3/s as demand (m3/s) = demand_from_datafile * demand_conv
    // Must be called before reading the file.
    if (unit == "lps") {
        demand_conv = 1 / 1000.;
        diameter_conv = 1 / 1000.;
        length_conv = 1.;
    } else {
        if (unit == "m3/h") {
            demand_conv = 1. / 3600.;
            diameter_conv = 1 / 1000.;
            length_conv = 1.;
        } else {
            if (unit == "gpm") {
                // 1 gallon=3.79 litre
                cout << endl << "EpanetReader::readFromFile unit detected: gpm, converting to m3/h" << endl;
                demand_conv = 3.79 / 1000. / 60.; //3.79 / 1000. / 60. is the conversion;
                diameter_conv = 0.0254;
                length_conv = 0.3048;
            } else {
                cout << endl <<
                        "ERROR: UNKNOWN conversion factor for EpanetReader::readFromFile. Possible choices: lps, m3/h, gpm";
                cin.get();
            }
        }
    }
    readFromFile(filename);
}

void EpanetReader::readFromFile(const std::string &filename) {
    ifstream file(filename);
    string line;
    string section;

    while (std::getline(file, line)) {
        // cout << "len=" << line.length() << "\t";
        if (line.find("[JUNCTIONS]") != string::npos) {
            section = "JUNCTIONS";
        } else if (line.find("[RESERVOIRS]") != string::npos) {
            section = "RESERVOIRS";
        } else if (line.find("[TANKS]") != string::npos) {
            section = "TANKS";
        } else if (line.find("[PIPES]") != string::npos) {
            section = "PIPES";
        } else if (line.find("[PUMPS]") != string::npos) {
            section = "PUMPS";
        } else if (line.find("[VALVES]") != string::npos) {
            section = "VALVES";
            cout << "Valve found!\n";
        } else if (line.find("[TAGS]") != string::npos) {
            section = "TAGS";
            cout << "Tags found!\n";
        } else if (line.find("[DEMANDS]") != string::npos) {
            section = "DEMANDS";
            cout << "DEMANDS found!\n";
        } else if (line.find("[STATUS]") != string::npos) {
            section = "STATUS";
            cout << "STATUS found!\n";
        } else if (line.find("[PATTERNS]") != string::npos) {
            section = "PATTERNS";
            cout << "PATTERNS found!\n";
        } else if (line.find("[CURVES]") != string::npos) {
            section = "CURVES";
            cout << "CURVES found!\n";
        } else if (line.find("[CONTROLS]") != string::npos) {
            section = "[CONTROLS]";
            cout << "[CONTROLS] found!\n";
        } else if (line.find("[RULES]") != string::npos) {
            section = "[RULES]";
            cout << "[RULES] found!\n";
        } else if (line.find("[ENERGY]") != string::npos) {
            section = "[ENERGY]";
            cout << "[ENERGY] found!\n";
        } else if (line.find("[EMITTERS]") != string::npos) {
            section = "[EMITTERS]";
            cout << "[EMITTERS] found!\n";
        } else if (line.find("[QUALITY]") != string::npos) {
            section = "C[QUALITY]";
            cout << "[QUALITY] found!\n";
        } else if (line.find("[SOURCES]") != string::npos) {
            section = "[SOURCES]";
            cout << "[SOURCES] found!\n";
        } else if (line.find("[REACTIONS]") != string::npos) {
            section = "[REACTIONS]";
            cout << "[REACTIONS] found!\n";
        } else if (line.find("[MIXING]") != string::npos) {
            section = "[MIXING]";
            cout << "[MIXING] found!\n";
        } else if (line.find("[TIMES]") != string::npos) {
            section = "[TIMES]";
            cout << "[TIMES] found!\n";
        } else if (line.find("[REPORT]") != string::npos) {
            section = "[REPORT]";
            cout << "[REPORT] found!\n";
        } else if (line.find("[OPTIONS]") != string::npos) {
            section = "[OPTIONS]";
            cout << "[OPTIONS] found!\n";
        } else if (line.find("[COORDINATES]") != string::npos) {
            section = "[COORDINATES]";
            cout << "[COORDINATES] found!\n";
        } else if (line.find("[VERTICES]") != string::npos) {
            section = "[]";
            cout << "[VERTICES] found!\n";
        } else if (line.find("[LABELS]") != string::npos) {
            section = "[LABELS]";
            cout << "[LABELS] found!\n";
        } else if (line.find("[BACKDROP]") != string::npos) {
            section = "[BACKDROP]";
            cout << "[BACKDROP] found!\n";
        } else if (line.length() <= 1 || line[0] == ';') {
            cout << "Line skipped"
                    << "\n";
            continue; // Skip empty lines or comments
        } else {
            // cout << "Section = " << section<< "\n";
            if (section == "JUNCTIONS") {
                JunctionReader junction = parseJunction(line);
                cout << "Junction ID = " << junction.ID << "\telev = " << junction.Elev << "\tDemand = " << junction.
                        Demand << "\t pattern = " << junction.Pattern << "\n";

                junctions.push_back(junction);
            } else if (section == "RESERVOIRS") {
                ReservoirReader reservoir = parseReservoir(line);
                cout << "Reservoir ID = " << reservoir.ID << "\thead = " << reservoir.Head << "\tPattern = " <<
                        reservoir.Pattern << "\n";

                reservoirs.push_back(reservoir);
            } else if (section == "TANKS") {
                TankReader tank = parseTank(line);
                cout << "Tank ID = " << tank.ID << "\tCurve = " << tank.VolCurve << "\t Overflow: " << tank.Overflow <<
                        "\n";

                tanks.push_back(tank);
            } else if (section == "PIPES") {
                PipeReader pipe = parsePipe(line);
                if (pipe.Status == "Closed") {
                    cout << endl << "Pipe " << pipe.ID << "is closed, skipping...." << endl;
                    //cin.get();
                } else {
                    cout << "Pipe ID = " << pipe.ID << "\tNode1 = " << pipe.Node1 << "\t Node2: " << pipe.Node2 << "\n";

                    pipes.push_back(pipe);
                }
            } else if (section == "PUMPS") {
                PumpReader pump = parsePump(line);
                cout << "Pump ID = " << pump.ID << "\tNode1 = " << pump.Node1 << "\t Node2: " << pump.Node2 <<
                        "\tPars: " << pump.Parameters << "\n";

                pumps.push_back(pump);
            } else if (section == "VALVES") {
                ValveReader valve = parseValve(line);
                cout << "Valve ID = " << valve.ID << "\tNode1 = " << valve.Node1 << "\t Node2: " << valve.Node2 << "\n";

                valves.push_back(valve);
            } else {
                cout << "Line skipped in section " << section << "\n";
                continue; // Skip empty lines or comments
            }
        }
    }
}

void EpanetReader::convertToRunner2() {
    // allow the uniform treatment of junctionsans reservoirs
    unifyJunctions();
    cout << "Junctions unified!\n";

    for (int i = 0; i < pipes.size(); i++) // go through all the pipes
    {
        double a = pipes[i].SpeedOfSound; //m/s
        double L = pipes[i].Length * length_conv; //m
        double D = pipes[i].Diameter * diameter_conv; //m
        //double lambda = pipes[i].Lambda; //hs model
        double lambda = pipes[i].Roughness; //hw model
        int idx_e = findNodeByID(pipes[i].Node1);
        int idx_v = findNodeByID(pipes[i].Node2);
        double he = junctions[idx_e].Elev; //m
        double hv = junctions[idx_v].Elev; //m

        if (junctions[idx_e].type == 1) //reservoir/tank at the start of pipe
        {
            he = hv;
        }
        if (junctions[idx_v].type == 1) //reservoir/tank at the end of pipe
        {
            hv = he;
        }

        if (true) {
            edges.push_back(new SCP(pipes[i].ID, pipes[i].Node1, pipes[i].Node2, 1000, a, L, D, lambda, he, hv, true));
            edges[i]->Set_string_prop("lambda_model", "hw");
            cout << "Pipe Call SCP (" << pipes[i].ID << ",\t" << pipes[i].Node1 << ",\t" << pipes[i].Node2 << ",\t1000,"
                    << a << ",\tL=" << L << ",\tD=" << D << ",\tlambda=" << lambda << ",\the=" << he << ",\thv=" << hv
                    << ",false)" << endl;
        }
        /*else if(pipes[i].type == "CV")
        {
            double zetaMin = 0.0;
            double zetaMax = 100000.0;
            edges.push_back(new CheckValve(pipes[i].ID, pipes[i].Node1, pipes[i].Node2, 1000, zetaMin, zetaMax, true));

        } */
    }
    cout << "SCP pipes ready!\n";

    // go through all the pipes and find their respective nodes
    for (int i = 0; i < pipes.size(); i++) {
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

    // go through all the valves and find their respective nodes
    for (int i = 0; i < valves.size(); i++) {
        // find the node index of pipe ends
        int idx1 = findNodeByID(valves[i].Node1); // idx is start (1)
        int idx2 = findNodeByID(valves[i].Node2); // idx is end (0)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(false);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(true);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);
    }

    // go through all the valves and find their respective nodes
    for (int i = 0; i < pumps.size(); i++) {
        // find the node index of pipe ends
        int idx1 = findNodeByID(pumps[i].Node1); // idx is start (1)
        int idx2 = findNodeByID(pumps[i].Node2); // idx is end (0)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(false);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(true);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);
    }


    cout << "Pipe end nodes found!\n";

    // create the connectivities
    for (int i = 0; i < junctions.size(); i++) {
        int size = junctions[i].idxPipe.size();
        cout << "i = " << i << "\tsize =" << size << "\t";

        if (size == 0) {
            cout << endl << endl << " Warning! junction[" << i << "].ID = " << junctions[i].ID;
            cout << endl << "Orphan node: no pipe connected." << endl;
            cin.get();
            //while(true) 1;
        }
        if (size == 1) {
            // Connect pipe-BC
            int idx1 = junctions[i].idxPipe[0];
            bool end1 = junctions[i].end[0];
            string name = junctions[i].ID;

            // check type
            if (junctions[i].type == 0) //"free" end
            {
                double demand = junctions[i].Demand * demand_conv;
                double D = pipes[idx1].Diameter * diameter_conv;
                double A = 0.25 * D * D * 3.14159265;
                double v = demand / A; //????
                if (!end1) v = -v;
                cons.push_back(new Connector(name, edges[idx1], !end1, "Velocity", v, demand, true));
                cout << "Node " << name << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ", Velocity: " << v
                        << ") - demand: " << demand << "\n";
                //cout<<"demand = "<<demand<<endl;
                //cin.get();
            }
            if (junctions[i].type == 1) // reservoir or tank
            {
                double head = junctions[i].Head * length_conv;
                int idx_other = findNodeByID(getOtherNodeOfPipe(idx1, junctions[i].ID));
                double elev = junctions[idx_other].Elev * length_conv;
                // Reservoirs, tanks
                double p = (head - elev) * 1000.0 * 9.81; //????
                // Head
                //double p = head*1000.0*9.81; //????
                double demand = 0;
                cons.push_back(new Connector(name, edges[idx1], !end1, "Pressure", p, demand, true));
                cout << "Node " << name << "\t call Connector(" << pipes[idx1].ID << "," << end1 << ",Pressure," << p <<
                        ") - demand: " << demand << "\n";
                cout << "\t head=" << head << ", elevation=" << elev << ", head-elevation=" << head - elev << "\n";
                cout << "\t p=" << p / 1000 / 9.81 << " vom" << endl;
                cout << "demand = " << demand << endl;

                //cin.get();
            }
        }

        if (size > 1) {
            string name = junctions[i].ID;
            vector<PSToolboxBaseEdge *> edges_act;
            vector<bool> front_act;
            vector<int> idx_act;
            double demand = junctions[i].Demand * demand_conv;

            for (int j = 0; j < size; j++) {
                int idx = junctions[i].idxPipe[j];
                edges_act.push_back(edges[idx]);
                front_act.push_back(!junctions[i].end[j]);
                idx_act.push_back(idx);
            }

            cons.push_back(new Connector(name, edges, front_act, idx_act, demand, false));
            cout << "Node " << name << "\t call Connector \t size = " << size << " \t" << ")\n";
            //cout<<"demand = "<<demand<<endl;
            //            cin.get();
        }

        /*if (size == 2)
        {
            // Connect pipe-pipe
            int idx1 = junctions[i].idxPipe[0];
            int idx2 = junctions[i].idxPipe[1];
            bool end1 = junctions[i].end[0];
            bool end2 = junctions[i].end[1];
            double demand = junctions[i].Demand / 3600.0;
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
            double demand = junctions[i].Demand / 3600.0;
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
            double demand = junctions[i].Demand / 3600.0;
            string name = junctions[i].ID;
            vector<int> idx_edge;
            idx_edge.push_back(idx1);
            idx_edge.push_back(idx2);
            idx_edge.push_back(idx3);
            idx_edge.push_back(idx4);
            cons.push_back(new Connector(name,edges[idx1],!end1,edges[idx2],!end2,edges[idx3],!end3,edges[idx4],!end4,demand,false,idx_edge));
            cout << "Node " << junctions[i].ID << "\t call Connector(" << pipes[idx1].ID << "," << end1 << "," << pipes[idx2].ID << "," << end2 << "," << pipes[idx3].ID << "," << end3  << "," << pipes[idx4].ID << "," << end4  << "," << demand << ")\n";
        }*/
    }

    cout << "Connections ready!\n";
}

void EpanetReader::convertToBioRunner(const string &bio_type) {
    // allow the uniform treatment of junctions and reservoirs
    unifyJunctions();
    cout << "Junctions unified!\n";

    for (auto &pipe: pipes) // go through all the pipes
    {
        double L = pipe.Length * length_conv; //m
        double D = pipe.Diameter * diameter_conv; //m
        //int idx_e = findNodeByID(pipes[i].Node1);
        //int idx_v = findNodeByID(pipes[i].Node2);
        //        double he = junctions[idx_e].Elev; //m
        //        double hv = junctions[idx_v].Elev; //m

        //        if (junctions[idx_e].type == 1) //reservoir/tank at the start of pipe
        //        {
        //            he = hv;
        //        }
        //        if (junctions[idx_v].type == 1) //reservoir/tank at the end of pipe
        //        {
        //            hv = he;
        //        }

        if (true) {
            edges.push_back(new BioPipe(pipe.ID, pipe.Node1, pipe.Node2, L, D, bio_type, true));
            cout << "Pipe Call BioPipe (" << pipe.ID << ",\t" << pipe.Node1 << ",\t" << pipe.Node2 <<
                    ",\tL=" << L << ",\tD=" << D << ",\tbio_type=" << bio_type << ",\t,false)" << endl;
        }
    }
    cout << "BioPipes ready!\n";

    // go through all the pipes and find their respective nodes
    for (int i = 0; i < pipes.size(); i++) {
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

    // go through all the valves and find their respective nodes
    for (int i = 0; i < valves.size(); i++) {
        // find the node index of pipe ends
        int idx1 = findNodeByID(valves[i].Node1); // idx is start (1)
        int idx2 = findNodeByID(valves[i].Node2); // idx is end (0)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(false);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(true);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);

        cout << endl << "ERROR! Found Valve " << valves[i].ID << ". For BioPipe computations, this must be removed." <<
                endl << endl;
        cin.get();
    }

    // go through all the pumps and find their respective nodes
    for (int i = 0; i < pumps.size(); i++) {
        // find the node index of pipe ends
        int idx1 = findNodeByID(pumps[i].Node1); // idx is start (1)
        int idx2 = findNodeByID(pumps[i].Node2); // idx is end (0)

        junctions[idx1].idxPipe.push_back(i);
        junctions[idx1].end.push_back(false);
        junctions[idx2].idxPipe.push_back(i);
        junctions[idx2].end.push_back(true);

        con_at_edge_start.push_back(idx1);
        con_at_edge_end.push_back(idx2);
        cout << endl << "ERROR! Found Pump " << pumps[i].ID << ". For BioPipe computations, this must be removed." <<
                endl << endl;
        cin.get();
    }

    cout << "Pipe end nodes found!\n";

    // create the connectivities
    for (int i = 0; i < junctions.size(); i++) {
        int size = junctions[i].idxPipe.size();
        cout << "junction " << junctions[i].ID << " (i = " << i << ")\tnumber of connected pipes =" << size << "\t";

        if (size == 0) {
            cout << endl << endl << " Warning! junction[" << i << "].ID = " << junctions[i].ID;
            cout << endl << "Orphan node: no pipe connected." << endl;
            cin.get();
            //while(true) 1;
        }
        if (size == 1) {
            // Connect pipe-BC
            int idx1 = junctions[i].idxPipe[0];
            bool end1 = junctions[i].end[0];
            string name = junctions[i].ID;

            // check type
            if (junctions[i].type == 0) //"free" end
            {
                double demand = junctions[i].Demand * demand_conv;
                double D = pipes[idx1].Diameter * diameter_conv;
                double A = 0.25 * D * D * 3.14159265;
                double v = demand / A; //????
                if (!end1) v = -v;
                bio_cons.push_back(new BioConnector(name, edges[idx1], !end1, "Velocity", v, demand, true));
                cout << "Node " << name << "\t call BioConnector(" << pipes[idx1].ID << "," << end1 << ",Velocity," << v
                        << ") - demand: " << demand << "\n";
                cin.get();
            }
            if (junctions[i].type == 1) // reservoir or tank
            {
                double head = junctions[i].Head * length_conv;
                int idx_other = findNodeByID(getOtherNodeOfPipe(idx1, junctions[i].ID));
                double elev = junctions[idx_other].Elev * length_conv;
                // Reservoirs, tanks
                double p = (head - elev) * 1000.0 * 9.81; //????
                // Head
                //double p = head*1000.0*9.81; //????
                double demand = 0;
                bio_cons.push_back(new BioConnector(name, edges[idx1], !end1, "Pressure", p, demand, true));
                cout << "Node " << name << "\t call BIOConnector(" << pipes[idx1].ID << "," << end1 << ",Pressure," << p
                        << ") - demand: " << demand << "\n";
            }
        }

        if (size > 1) {
            string name = junctions[i].ID;
            vector<PSToolboxBaseEdge *> edges_act;
            vector<bool> front_act;
            vector<int> idx_act;
            double demand = junctions[i].Demand * demand_conv;

            for (int j = 0; j < size; j++) {
                int idx = junctions[i].idxPipe[j];
                edges_act.push_back(edges[idx]);
                front_act.push_back(!junctions[i].end[j]);
                idx_act.push_back(idx);
            }

            bio_cons.push_back(new BioConnector(name, edges, front_act, idx_act, demand, false));
            cout << "Node " << name << "\t call Connector \t size = " << size << " \t" << ")\n";
        }
    }
    cout << "BioConnections ready!\n";

    //cin.get();
}


void EpanetReader::PrintData() {
    cout << endl << endl << "------------------CPP FILE--------------------" << endl << endl;
    cout << "Size of edges:\t" << edges.size() << endl;
    for (int i = 0; i < edges.size(); i++) {
        cout << "edge i = " << i << "\t" << edges[i]->name << endl;
    }

    cout << " Size of cons:\t" << cons.size() << endl;
    for (int i = 0; i < cons.size(); i++) {
        cout << "connector i = " << i << "\t BC: " << cons[i]->BC_type << "\t size of edgeidx: " << cons[i]->edges_idx.
                size() << endl;
    }

    cout << "Size of start:\t" << con_at_edge_start.size() << endl;
    for (int i = 0; i < con_at_edge_start.size(); i++) {
        cout << "i = " << i << "\t" << con_at_edge_start[i] << endl;
    }

    cout << "  Size of end:\t" << con_at_edge_end.size() << endl;
    for (int i = 0; i < con_at_edge_end.size(); i++) {
        cout << "i = " << i << "\t" << con_at_edge_end[i] << endl;
    }
}

int EpanetReader::nextPipeAtNode(int idx) {
    if (junctions[idx].idxPipe[0] == -1)
        return 0;
    if (junctions[idx].idxPipe[1] == -1)
        return 1;
    if (junctions[idx].idxPipe[2] == -1)
        return 2;
    return -1;
}

string EpanetReader::getOtherNodeOfPipe(int idxPipe, string Node) {
    if (pipes[idxPipe].Node1 == Node)
        return pipes[idxPipe].Node2;
    else
        return pipes[idxPipe].Node1;
}

void EpanetReader::printEdgesAndCons() {
    for (auto &edge: edges)
        cout << edge->Info();
}

void EpanetReader::unifyJunctions() {
    for (int i = 0; i < junctions.size(); i++) {
        junctions[i].type = 0;
    }


    for (int i = 0; i < reservoirs.size(); i++) {
        JunctionReader J;
        J.type = 1;
        J.Head = reservoirs[i].Head;
        J.ID = reservoirs[i].ID;
        junctions.push_back(J);
        cout << " Reservoir " << J.ID << " converted to node with constant (total) head " << J.Head << "m." << endl <<
                endl;
    }

    for (int i = 0; i < tanks.size(); i++) {
        JunctionReader J;
        J.type = 1;
        J.Head = tanks[i].Elevation + tanks[i].InitLevel;
        J.ID = tanks[i].ID;
        junctions.push_back(J);
        cout << " Tank " << J.ID << " converted to node with constant (total) head " << J.Head << "m." << endl << endl;
        //cin.get();
    }


    cout << "Njunctions = " << junctions.size() << "\n";
}

int EpanetReader::findNodeByID(const std::string ID) {
    for (int i = 0; i < junctions.size(); i++) {
        if (junctions[i].ID == ID) {
            return i;
        }
    }
    cout << "Node " << ID << " not found!" << endl;
    return -1;
}

int EpanetReader::findPipeByID(const std::string ID) {
    for (int i = 0; i < pipes.size(); i++) {
        if (pipes[i].ID == ID) {
            return i;
        }
    }
    cout << "Pipe " << ID << " not found!" << endl;
    return -1;
}


int EpanetReader::findReservoirByID(const std::string ID) {
    for (int i = 0; i < reservoirs.size(); i++) {
        if (reservoirs[i].ID == ID) {
            return i;
        }
    }
    return -1;
}

int EpanetReader::findPumpByID(const std::string ID) {
    for (int i = 0; i < pumps.size(); i++) {
        if (pumps[i].ID == ID) {
            return i;
        }
    }
    return -1;
}

bool EpanetReader::checkIfIncludedInVector(vector<int> list, int elem) {
    for (int i = 0; i < list.size(); i++) {
        if (list[i] == elem)
            return true;
    }
    return false;
}

vector<int> EpanetReader::findConnectingPipes(const std::string node) {
    vector<int> connectingPipes;
    for (int i = 0; i < pipes.size(); i++) {
        if (node == pipes[i].Node1 || node == pipes[i].Node2) {
            // cout << "Found at " << i << "\n";
            connectingPipes.push_back(i);
        }
    }
    return connectingPipes;
}

vector<int> EpanetReader::findConnectingFrontPipes(const std::string node) {
    vector<int> connectingPipes;
    for (int i = 0; i < pipes.size(); i++) {
        if (node == pipes[i].Node1) {
            // cout << "Found at " << i << "\n";
            connectingPipes.push_back(i);
        }
    }
    return connectingPipes;
}

vector<int> EpanetReader::findConnectingBackPipes(const std::string node) {
    vector<int> connectingPipes;
    for (int i = 0; i < pipes.size(); i++) {
        if (node == pipes[i].Node2) {
            // cout << "Found at " << i << "\n";
            connectingPipes.push_back(i);
        }
    }
    return connectingPipes;
}

JunctionReader EpanetReader::parseJunction(const std::string &line) {
    istringstream iss(line);
    JunctionReader junction;
    iss >> junction.ID >> junction.Elev >> junction.Demand >> junction.Pattern;
    return junction;
}

ReservoirReader EpanetReader::parseReservoir(const std::string &line) {
    istringstream iss(line);
    ReservoirReader reservoir;
    iss >> reservoir.ID >> reservoir.Head >> reservoir.Pattern;
    return reservoir;
}

TankReader EpanetReader::parseTank(const std::string &line) {
    istringstream iss(line);
    TankReader tank;
    iss >> tank.ID >> tank.Elevation >> tank.InitLevel >> tank.MinLevel >> tank.MaxLevel >> tank.Diameter >> tank.MinVol
            >> tank.VolCurve >> tank.Overflow;
    return tank;
}

PipeReader EpanetReader::parsePipe(const std::string &line) {
    istringstream iss(line);
    PipeReader pipe;
    string Status;
    iss >> pipe.ID >> pipe.Node1 >> pipe.Node2 >> pipe.Length >> pipe.Diameter >> pipe.Roughness >> pipe.MinorLoss >>
            pipe.Status;

    // Closed pipes are handled in the caller function, no need to do anything here.

    return pipe;
}

PumpReader EpanetReader::parsePump(const std::string &line) {
    istringstream iss(line);
    PumpReader pump;
    string tmp;
    iss >> pump.ID >> pump.Node1 >> pump.Node2 >> tmp >> pump.Parameters;
    return pump;
}

ValveReader EpanetReader::parseValve(const std::string &line) {
    istringstream iss(line);
    ValveReader valve;
    string tmp;
    iss >> valve.ID >> valve.Node1 >> valve.Node2 >> valve.Diameter >> valve.Type >> valve.Setting >> valve.MinorLoss;
    return valve;
}

