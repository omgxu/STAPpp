#pragma once

#include "Element.h"
#include <cmath>

using namespace std;

//! Bar element class
class CQ9 : public CElement
{
public:
    //!	Constructor
    CQ9();

    //!	Desconstructor
    ~CQ9();

    //!	Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //!	Write element data to stream
    virtual void Write(COutputter& output);

    //! Generate location matrix: the global equation number that corresponding to each DOF of the
    //! element
    //	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

    //!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //!	Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement);

    //!	Return the size of the element stiffness matrix (stored as an array column by column)
    virtual unsigned int SizeOfStiffnessMatrix();

    //! Calculate Gauss point displacement
	virtual void CalculateGaussPointDisplacement(double* gaussDisp, double *Displacement);

#ifdef _VIB_
//!	Calculate element mass matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementMass(double* mass); 
#endif

};