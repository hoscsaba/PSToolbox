#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include "BioConnector.h"

using namespace std;

struct JunctionReader {
    string ID;
    double Elev;
    double Head;
    double Demand;
    string Pattern;
    vector<int> idxPipe;
    vector<bool> end;
    int type; // 0:junction, 1:reservoir
};

struct ReservoirReader {
    string ID;
    double Head;
    int Pattern;
};

struct TankReader {
    string ID;
    double Elevation;
    int InitLevel;
    int MinLevel;
    double MaxLevel;
    double Diameter;
    double MinVol;
    string VolCurve;
    string Overflow;
};

struct PipeReader {
    string ID;
    string Node1;
    string Node2;
    double Length;
    double Diameter;
    double Roughness;
    double Lambda;
    double Delta;
    double SpeedOfSound;
    double MinorLoss;
    string Status;
};

struct PumpReader {
    string ID;
    string Node1;
    string Node2;
    string Parameters;
};

struct ValveReader {
    string ID;
    string Node1;
    string Node2;
    double Diameter;
    string Type;
    int Setting;
    double MinorLoss;
};

class EpanetReader {
public:
    vector<JunctionReader> junctions;
    vector<ReservoirReader> reservoirs;
    vector<TankReader> tanks;
    vector<PipeReader> pipes;
    vector<PumpReader> pumps;
    vector<ValveReader> valves;

    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> cons;
    vector<BioConnector *> bio_cons;
    vector<int> con_at_edge_start;
    vector<int> con_at_edge_end;

    double demand_conv = 1.;
    double diameter_conv = 1.;
    double length_conv = 1.;

    void readFromFile(const std::string &filename);

    void readFromFile(const std::string &filename, const string &demand_unit);

    void convertToRunner2();

    void convertToBioRunner(const string &bio_type);

    int nextPipeAtNode(int idx);

    string getOtherNodeOfPipe(int idxPipe, string Node);

    void printEdgesAndCons();

    void PrintData();

    void unifyJunctions();

    int findNodeByID(const std::string ID);

    int findPipeByID(const std::string ID);

    int findPumpByID(const std::string ID);

    int findReservoirByID(const std::string ID);

    bool checkIfIncludedInVector(vector<int> list, int elem);

    vector<int> findConnectingPipes(const std::string node);

    vector<int> findConnectingFrontPipes(const std::string node);

    vector<int> findConnectingBackPipes(const std::string node);

    JunctionReader parseJunction(const std::string &line);

    ReservoirReader parseReservoir(const std::string &line);

    TankReader parseTank(const std::string &line);

    PipeReader parsePipe(const std::string &line);

    PumpReader parsePump(const std::string &line);

    ValveReader parseValve(const std::string &line);
};

