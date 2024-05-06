#pragma once

#include <vector>
#include <iostream>
#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"
using namespace std;

class PSToolboxPlotter
{
  private:
    std::vector<double> column1, column2, column3, column4, column5, column6, column7, column8, column9, column10;

  public:
    PSToolboxPlotter(string fname);
    ~PSToolboxPlotter();

    void Plot();
};
  
