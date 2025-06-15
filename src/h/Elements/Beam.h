

//by hly




#pragma once

#include "Element.h"

using namespace std;

//! Beam element class
class CBeam : public CElement
{
public:

	//!	Constructor
	CBeam();

	//!	Desconstructor
	~CBeam();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output);

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};
