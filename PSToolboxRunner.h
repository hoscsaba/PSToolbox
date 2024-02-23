#include <vector>
#include "PSToolboxBaseEdge.h"
#include "Connector.h"

class PSToolboxRunner
{
  public:
    PSToolboxRunner(vector<PSToolboxBaseEdge *>, vector<Connector *>,vector<int>,vector<int>);
    void Run(double);
    void Set_Save_data(bool newval){save_data=newval;};
  private:
    vector<PSToolboxBaseEdge *> edges;
    vector<Connector *> cons;
    vector<int> con_at_edge_start;
    vector<int> con_at_edge_end;
    bool save_data;
};
