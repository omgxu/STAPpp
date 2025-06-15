/*****************************************************************************/
/*  STAP++ : MindlinPlate element made by Yu Jing                            */
/*****************************************************************************/

#pragma once

#include "Element.h"
#include <vector>

using namespace std;

//! MindlinPlate element class
class CMindlinPlate : public CElement
{
public:

//!	Constructor
	CMindlinPlate();

//!	Desconstructor
	~CMindlinPlate();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};

void gauss(int ngp, std::vector<double>& w, std::vector<double>& gp);