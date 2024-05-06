#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "PSToolboxPlotter.h"
#include "/Users/hoscsaba/program/matplotlib-cpp/matplotlibcpp.h"

using namespace std; 
namespace plt=matplotlibcpp; 

 PSToolboxPlotter::PSToolboxPlotter(string fname){

  std::ifstream inputFile(fname); // Assuming the file is named "data.txt"

  if (!inputFile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
   cin.get(); 
  }

  // Read the header line
  std::string headerLine;
  std::getline(inputFile, headerLine);

  // Read each line and extract values
  std::string line;
  while (std::getline(inputFile, line)) {
    std::stringstream ss(line);
    std::string value;

    // Assuming the values are separated by semicolons
    std::getline(ss, value, ';');
    column1.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column2.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column3.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column4.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column5.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column6.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column7.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column8.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column9.push_back(std::stod(value));

    std::getline(ss, value, ';');
    column10.push_back(std::stod(value));
  }
 }

PSToolboxPlotter::~PSToolboxPlotter(){}

void PSToolboxPlotter::Plot(){
  plt::clf();
  plt::subplot(2, 1, 1);
  plt::plot(column1,column2,"-");
  plt::plot(column1,column3,"--");
  plt::ylabel("p, bar");

  plt::subplot(2, 1, 2);
  plt::plot(column1, column4,"-");
  plt::plot(column1,column5,"--");
  plt::ylabel("v, m/s");

  plt::show();
};

