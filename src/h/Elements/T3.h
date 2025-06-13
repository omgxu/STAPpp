/*****************************************************************************/
/*  STAP++ : T3 element made by Yu Jing                                      */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! T3 element class
class CT3 : public CElement
{
public:
    //!	Constructor
    CT3();

    //!	Desconstructor
    ~CT3();

    //!	Read element data from stream Input
    // virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //!	Write element data to stream
    // virtual void Write(COutputter& output, unsigned int Ele);
    virtual void Write(COutputter& output);

    //! Generate location matrix: the global equation number that corresponding to each DOF of the element
	virtual void GenerateLocationMatrix();

    //!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //!	Calculate element stress
    /* virtual void ElementStress(double* stress, double* Displacement,
                               double* GaussPosition = nullptr,
                               double* GaussDisplacements = nullptr, double* weights = nullptr); */
    virtual void ElementStress(double* stress, double* Displacement);

	// //!	Calculate the values required in the POSTPROCESS 
	// virtual void ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions);

    // //!	Return the size of the element stiffness matrix (stored as an array column by column)
    // virtual unsigned int SizeOfStiffnessMatrix();

#ifdef _VIB_
//!	Calculate element mass matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementMass(double* mass); 
#endif

};