#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Domain.h>
#include <Element.h>
#include <ElementGroup.h>
#include <Node.h>
#include "Elements/T3.h"
#include "Elements/Bar.h"
//#include "Elements/Q4.h"
using namespace std;
// the main function to output a data file for Tecplot
void TecplotOutPutter();

// write nodal and element data to the Tecplot data file
void OutputTecplotDatafile(vector<vector<double> >& nodalData);

// initialise nodal data, add magnified nodal displacements to coordinates
vector<vector<double>> getNodalData();

// get element data for T3 elements
vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement);

// get element data for Q4 elements
//vector<vector<double>> getQ4ElementStressData(CElementGroup& eleGrp, double* displacement);

// calculate nodal stress for T3 elements
void T3NodalStress(const vector<vector<double>> elementdata, vector<vector<double>>& nodestress);

// calculate nodal stress for Q4 elements
//void Q4NodalStress(const vector<vector<double>> elementdata, vector<vector<double>>& nodestress);