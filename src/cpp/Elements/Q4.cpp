#include "Elements/Q4.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Constants for Gauss quadrature
const double GAUSS_POINT_COORD = 1.0 / sqrt(3.0);
const int NUM_GAUSS_POINTS = 4;
const int NUM_NODES = 4;
const int DOF_PER_NODE = 3;

// Constructor - Initialize member variables
CQ4::CQ4()
{
    NEN_ = NUM_NODES;  // Number of nodes per element
    nodes_ = new CNode*[NEN_];
    
    ND_ = NUM_NODES * DOF_PER_NODE;  // Total degrees of freedom
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

// Destructor
CQ4::~CQ4()
{
}

/**
 * Read element data from input stream
 * @param Input Input file stream
 * @param MaterialSets Material properties
 * @param NodeList List of nodes
 * @return True if successful
 */
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int materialSetIndex;
    unsigned int nodeIndices[NUM_NODES];

    Input >> nodeIndices[0] >> nodeIndices[1] >> nodeIndices[2] >> nodeIndices[3] >> materialSetIndex;

    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + materialSetIndex - 1;
    
    for(int i = 0; i < NUM_NODES; i++) {
        nodes_[i] = NodeList + nodeIndices[i] - 1;
    }

    return true;
}

/**
 * Write element data to output stream
 * @param output Output stream
 */ 
void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset
           << endl;
}

/**
 * Generate location matrix for element DOFs
 * Maps local DOF indices to global equation numbers
 * Note: Equation numbers start from 1
 */
void CQ4::GenerateLocationMatrix()
{
    unsigned int index = 0;
    for (unsigned int node = 0; node < NEN_; node++) {
        for (unsigned int dof = 0; dof < DOF_PER_NODE; dof++) {
            LocationMatrix_[index++] = nodes_[node]->bcode[dof];
        }
    }
}

// Return size of element stiffness matrix stored as array
unsigned int CQ4::SizeOfStiffnessMatrix() 
{ 
    return ND_ * (ND_ + 1) / 2; 
}

/**
 * Calculate B matrix and Jacobian determinant
 * @param B Output B matrix 
 * @param eta Natural coordinate eta
 * @param psi Natural coordinate psi
 * @param nodeCoordX X coordinates of nodes
 * @param nodeCoordY Y coordinates of nodes
 * @return Determinant of Jacobian matrix
 */
double CalculateStrainDisplacementMatrix(double B[8], const double eta, const double psi, 
                                      const double nodeCoordX[4], const double nodeCoordY[4])
{
    // Shape function derivatives
    const double shapeFuncDerivatives[8] = {
        (eta - 1) / 4, (1 - eta) / 4,  (1 + eta) / 4, (-eta - 1) / 4,  // dN/dEta
        (psi - 1) / 4, (-psi - 1) / 4, (1 + psi) / 4, (1 - psi) / 4    // dN/dPsi
    };

    // Calculate Jacobian matrix
    double jacobian[4] = {
        0.0, 0.0,  // J11, J12
        0.0, 0.0   // J21, J22
    };

    // Calculate Jacobian components
    for(int i = 0; i < NUM_NODES; i++) {
        jacobian[0] += shapeFuncDerivatives[i] * nodeCoordX[i];        // J11
        jacobian[1] += shapeFuncDerivatives[i] * nodeCoordY[i];        // J12
        jacobian[2] += shapeFuncDerivatives[i+4] * nodeCoordX[i];      // J21
        jacobian[3] += shapeFuncDerivatives[i+4] * nodeCoordY[i];      // J22
    }

    // Calculate determinant and inverse
    const double detJ = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
    const double invJ[4] = {
         jacobian[3] / detJ, -jacobian[1] / detJ,
        -jacobian[2] / detJ,  jacobian[0] / detJ
    };

    // Calculate B matrix components
    for(int i = 0; i < NUM_NODES; i++) {
        B[i] = invJ[0] * shapeFuncDerivatives[i] + invJ[1] * shapeFuncDerivatives[i+4];
        B[i+4] = invJ[2] * shapeFuncDerivatives[i] + invJ[3] * shapeFuncDerivatives[i+4];
    }

    return detJ;
}

/**
 * Accumulate stiffness matrix contributions for 2D 4Q element
 * @param eta Natural coordinate eta
 * @param psi Natural coordinate psi
 * @param weight Weight for Gauss quadrature
 * @param xe X coordinates of nodes
 * @param ye Y coordinates of nodes
 * @param ke Output stiffness matrix contributions
 * @param E Young's modulus
 * @param v Poisson's ratio
 */
void AccumulateEtaPsi(const double& eta, const double& psi, const double& weight, const double* xe,
                      const double* ye, double* ke, const double E, const double v)
{
    double B[8];
    double DetJe = CalculateStrainDisplacementMatrix(B, eta, psi, xe, ye);
    const double d33 = (1 - v) / 2.0;
    const double cof = E / (1 - v * v) * std::abs(DetJe) * weight;

    // see 4Q.nb and 4Q-form-key.py
    ke[0] += cof * (B[0] * B[0] + d33 * B[4] * B[4]);
    ke[1] += cof * (d33 * B[0] * B[4] + v * B[0] * B[4]);
    ke[2] += cof * (d33 * B[0] * B[0] + B[4] * B[4]);
    ke[3] += cof * (B[0] * B[1] + d33 * B[4] * B[5]);
    ke[4] += cof * (v * B[1] * B[4] + d33 * B[0] * B[5]);
    ke[5] += cof * (B[1] * B[1] + d33 * B[5] * B[5]);
    ke[6] += cof * (d33 * B[1] * B[4] + v * B[0] * B[5]);
    ke[7] += cof * (d33 * B[0] * B[1] + B[4] * B[5]);
    ke[8] += cof * (d33 * B[1] * B[5] + v * B[1] * B[5]);
    ke[9] += cof * (d33 * B[1] * B[1] + B[5] * B[5]);
    ke[10] += cof * (B[0] * B[2] + d33 * B[4] * B[6]);
    ke[11] += cof * (v * B[2] * B[4] + d33 * B[0] * B[6]);
    ke[12] += cof * (B[1] * B[2] + d33 * B[5] * B[6]);
    ke[13] += cof * (v * B[2] * B[5] + d33 * B[1] * B[6]);
    ke[14] += cof * (B[2] * B[2] + d33 * B[6] * B[6]);
    ke[15] += cof * (d33 * B[2] * B[4] + v * B[0] * B[6]);
    ke[16] += cof * (d33 * B[0] * B[2] + B[4] * B[6]);
    ke[17] += cof * (d33 * B[2] * B[5] + v * B[1] * B[6]);
    ke[18] += cof * (d33 * B[1] * B[2] + B[5] * B[6]);
    ke[19] += cof * (d33 * B[2] * B[6] + v * B[2] * B[6]);
    ke[20] += cof * (d33 * B[2] * B[2] + B[6] * B[6]);
    ke[21] += cof * (B[0] * B[3] + d33 * B[4] * B[7]);
    ke[22] += cof * (v * B[3] * B[4] + d33 * B[0] * B[7]);
    ke[23] += cof * (B[1] * B[3] + d33 * B[5] * B[7]);
    ke[24] += cof * (v * B[3] * B[5] + d33 * B[1] * B[7]);
    ke[25] += cof * (B[2] * B[3] + d33 * B[6] * B[7]);
    ke[26] += cof * (v * B[3] * B[6] + d33 * B[2] * B[7]);
    ke[27] += cof * (B[3] * B[3] + d33 * B[7] * B[7]);
    ke[28] += cof * (d33 * B[3] * B[4] + v * B[0] * B[7]);
    ke[29] += cof * (d33 * B[0] * B[3] + B[4] * B[7]);
    ke[30] += cof * (d33 * B[3] * B[5] + v * B[1] * B[7]);
    ke[31] += cof * (d33 * B[1] * B[3] + B[5] * B[7]);
    ke[32] += cof * (d33 * B[3] * B[6] + v * B[2] * B[7]);
    ke[33] += cof * (d33 * B[2] * B[3] + B[6] * B[7]);
    ke[34] += cof * (d33 * B[3] * B[7] + v * B[3] * B[7]);
    ke[35] += cof * (d33 * B[3] * B[3] + B[7] * B[7]);
}


// Convert 2D stiffness matrix to 3D stiffness matrix
void Convert2d23d(const double* k, double* Matrix, const double i[3], const double j[3])
{
    // to see how these are generated, see ../../memo/4Q.nb and 4Q2d23d.py
    Matrix[0] = i[0] * (i[0] * k[0] + j[0] * k[1]) + j[0] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[1] = i[1] * (i[1] * k[0] + j[1] * k[1]) + j[1] * (i[1] * k[1] + j[1] * k[2]);
    Matrix[2] = i[1] * (i[0] * k[0] + j[0] * k[1]) + j[1] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[3] = i[2] * (i[2] * k[0] + j[2] * k[1]) + j[2] * (i[2] * k[1] + j[2] * k[2]);
    Matrix[4] = i[2] * (i[1] * k[0] + j[1] * k[1]) + j[2] * (i[1] * k[1] + j[1] * k[2]);
    Matrix[5] = i[2] * (i[0] * k[0] + j[0] * k[1]) + j[2] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[6] = i[0] * (i[0] * k[5] + j[0] * k[8]) + j[0] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[7] = i[0] * (i[2] * k[3] + j[2] * k[4]) + j[0] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[8] = i[0] * (i[1] * k[3] + j[1] * k[4]) + j[0] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[9] = i[0] * (i[0] * k[3] + j[0] * k[4]) + j[0] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[10] = i[1] * (i[1] * k[5] + j[1] * k[8]) + j[1] * (i[1] * k[8] + j[1] * k[9]);
    Matrix[11] = i[1] * (i[0] * k[5] + j[0] * k[8]) + j[1] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[12] = i[1] * (i[2] * k[3] + j[2] * k[4]) + j[1] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[13] = i[1] * (i[1] * k[3] + j[1] * k[4]) + j[1] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[14] = i[1] * (i[0] * k[3] + j[0] * k[4]) + j[1] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[15] = i[2] * (i[2] * k[5] + j[2] * k[8]) + j[2] * (i[2] * k[8] + j[2] * k[9]);
    Matrix[16] = i[2] * (i[1] * k[5] + j[1] * k[8]) + j[2] * (i[1] * k[8] + j[1] * k[9]);
    Matrix[17] = i[2] * (i[0] * k[5] + j[0] * k[8]) + j[2] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[18] = i[2] * (i[2] * k[3] + j[2] * k[4]) + j[2] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[19] = i[2] * (i[1] * k[3] + j[1] * k[4]) + j[2] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[20] = i[2] * (i[0] * k[3] + j[0] * k[4]) + j[2] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[21] = i[0] * (i[0] * k[14] + j[0] * k[19]) + j[0] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[22] = i[0] * (i[2] * k[12] + j[2] * k[13]) + j[0] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[23] = i[0] * (i[1] * k[12] + j[1] * k[13]) + j[0] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[24] = i[0] * (i[0] * k[12] + j[0] * k[13]) + j[0] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[25] = i[0] * (i[2] * k[10] + j[2] * k[11]) + j[0] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[26] = i[0] * (i[1] * k[10] + j[1] * k[11]) + j[0] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[27] = i[0] * (i[0] * k[10] + j[0] * k[11]) + j[0] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[28] = i[1] * (i[1] * k[14] + j[1] * k[19]) + j[1] * (i[1] * k[19] + j[1] * k[20]);
    Matrix[29] = i[1] * (i[0] * k[14] + j[0] * k[19]) + j[1] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[30] = i[1] * (i[2] * k[12] + j[2] * k[13]) + j[1] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[31] = i[1] * (i[1] * k[12] + j[1] * k[13]) + j[1] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[32] = i[1] * (i[0] * k[12] + j[0] * k[13]) + j[1] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[33] = i[1] * (i[2] * k[10] + j[2] * k[11]) + j[1] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[34] = i[1] * (i[1] * k[10] + j[1] * k[11]) + j[1] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[35] = i[1] * (i[0] * k[10] + j[0] * k[11]) + j[1] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[36] = i[2] * (i[2] * k[14] + j[2] * k[19]) + j[2] * (i[2] * k[19] + j[2] * k[20]);
    Matrix[37] = i[2] * (i[1] * k[14] + j[1] * k[19]) + j[2] * (i[1] * k[19] + j[1] * k[20]);
    Matrix[38] = i[2] * (i[0] * k[14] + j[0] * k[19]) + j[2] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[39] = i[2] * (i[2] * k[12] + j[2] * k[13]) + j[2] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[40] = i[2] * (i[1] * k[12] + j[1] * k[13]) + j[2] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[41] = i[2] * (i[0] * k[12] + j[0] * k[13]) + j[2] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[42] = i[2] * (i[2] * k[10] + j[2] * k[11]) + j[2] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[43] = i[2] * (i[1] * k[10] + j[1] * k[11]) + j[2] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[44] = i[2] * (i[0] * k[10] + j[0] * k[11]) + j[2] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[45] = i[0] * (i[0] * k[27] + j[0] * k[34]) + j[0] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[46] = i[0] * (i[2] * k[25] + j[2] * k[26]) + j[0] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[47] = i[0] * (i[1] * k[25] + j[1] * k[26]) + j[0] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[48] = i[0] * (i[0] * k[25] + j[0] * k[26]) + j[0] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[49] = i[0] * (i[2] * k[23] + j[2] * k[24]) + j[0] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[50] = i[0] * (i[1] * k[23] + j[1] * k[24]) + j[0] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[51] = i[0] * (i[0] * k[23] + j[0] * k[24]) + j[0] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[52] = i[0] * (i[2] * k[21] + j[2] * k[22]) + j[0] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[53] = i[0] * (i[1] * k[21] + j[1] * k[22]) + j[0] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[54] = i[0] * (i[0] * k[21] + j[0] * k[22]) + j[0] * (i[0] * k[28] + j[0] * k[29]);
    Matrix[55] = i[1] * (i[1] * k[27] + j[1] * k[34]) + j[1] * (i[1] * k[34] + j[1] * k[35]);
    Matrix[56] = i[1] * (i[0] * k[27] + j[0] * k[34]) + j[1] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[57] = i[1] * (i[2] * k[25] + j[2] * k[26]) + j[1] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[58] = i[1] * (i[1] * k[25] + j[1] * k[26]) + j[1] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[59] = i[1] * (i[0] * k[25] + j[0] * k[26]) + j[1] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[60] = i[1] * (i[2] * k[23] + j[2] * k[24]) + j[1] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[61] = i[1] * (i[1] * k[23] + j[1] * k[24]) + j[1] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[62] = i[1] * (i[0] * k[23] + j[0] * k[24]) + j[1] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[63] = i[1] * (i[2] * k[21] + j[2] * k[22]) + j[1] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[64] = i[1] * (i[1] * k[21] + j[1] * k[22]) + j[1] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[65] = i[1] * (i[0] * k[21] + j[0] * k[22]) + j[1] * (i[0] * k[28] + j[0] * k[29]);
    Matrix[66] = i[2] * (i[2] * k[27] + j[2] * k[34]) + j[2] * (i[2] * k[34] + j[2] * k[35]);
    Matrix[67] = i[2] * (i[1] * k[27] + j[1] * k[34]) + j[2] * (i[1] * k[34] + j[1] * k[35]);
    Matrix[68] = i[2] * (i[0] * k[27] + j[0] * k[34]) + j[2] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[69] = i[2] * (i[2] * k[25] + j[2] * k[26]) + j[2] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[70] = i[2] * (i[1] * k[25] + j[1] * k[26]) + j[2] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[71] = i[2] * (i[0] * k[25] + j[0] * k[26]) + j[2] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[72] = i[2] * (i[2] * k[23] + j[2] * k[24]) + j[2] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[73] = i[2] * (i[1] * k[23] + j[1] * k[24]) + j[2] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[74] = i[2] * (i[0] * k[23] + j[0] * k[24]) + j[2] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[75] = i[2] * (i[2] * k[21] + j[2] * k[22]) + j[2] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[76] = i[2] * (i[1] * k[21] + j[1] * k[22]) + j[2] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[77] = i[2] * (i[0] * k[21] + j[0] * k[22]) + j[2] * (i[0] * k[28] + j[0] * k[29]);
}

/**
 * Convert 3D coordinates to 2D local coordinate system
 * @param nodes Array of 4 nodes of quadrilateral element
 * @param n Element normal vector (output)  
 * @param i Local x-axis direction (output)
 * @param j Local y-axis direction (output)
 * @param xe Node x-coordinates in local system (output)
 * @param ye Node y-coordinates in local system (output)
 */
void Convert3d22d4Q(CNode* const nodes[4], double n[3], double i[3], double j[3], 
                    double xe[4], double ye[4])
{
    // Get references to nodes
    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];
    const CNode& n4 = *nodes[3];

    // Calculate edge vectors
    const double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    const double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};

    // Calculate normal vector as normalized cross product p31 x p21
    n[0] = p31[1] * p21[2] - p31[2] * p21[1];
    n[1] = p31[2] * p21[0] - p31[0] * p21[2];
    n[2] = p31[0] * p21[1] - p31[1] * p21[0];
    normalize(n);

    // Set i-axis parallel to p21 (local x-axis)
    i[0] = p21[0];
    i[1] = p21[1];
    i[2] = p21[2];
    normalize(i);

    // Calculate j-axis as n x i (local y-axis)
    j[0] = n[1] * i[2] - n[2] * i[1];
    j[1] = n[2] * i[0] - n[0] * i[2];
    j[2] = n[0] * i[1] - n[1] * i[0];

    // Transform node coordinates to local system
    for(int k = 0; k < 4; k++) {
        const CNode& node = *nodes[k];
        xe[k] = i[0] * node.XYZ[0] + i[1] * node.XYZ[1] + i[2] * node.XYZ[2];
        ye[k] = j[0] * node.XYZ[0] + j[1] * node.XYZ[1] + j[2] * node.XYZ[2];
    }
}

/**
 * Calculate element stiffness matrix
 * Upper triangular matrix stored column by column from diagonal
 */
void CQ4::ElementStiffness(double* Matrix) 
{
    clear(Matrix, SizeOfStiffnessMatrix());

    // Convert coordinates to 2D system
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // Initialize stiffness matrix components
    double ke[36] = {0.0};  // 8x8 matrix stored as upper triangular

    // Setup Gauss quadrature points
    const double pos = 1.0 / std::sqrt(3.0);
    const double gaussPoints[2] = {-pos, pos};
    const double gaussWeights[2][2] = {{1.0, 1.0}, {1.0, 1.0}};

    // Get material properties
    const CQ4Material* material = static_cast<CQ4Material*>(ElementMaterial_);
    const double& E = material->E;
    const double& v = material->nu;

    // Accumulate contributions from all Gauss points
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            AccumulateEtaPsi(gaussPoints[i], gaussPoints[j], 
                           gaussWeights[i][j], xe, ye, ke, E, v);
        }
    }

    // Transform 2D stiffness to 3D system
    Convert2d23d(ke, Matrix, i, j);
}

/**
 * Calculate stress at a point
 * @param eta Natural coordinate eta
 * @param psi Natural coordinate psi 
 * @param xe Node x-coordinates in local system
 * @param ye Node y-coordinates in local system
 * @param E Young's modulus
 * @param v Poisson's ratio
 * @param de Local nodal displacements
 * @param stress Output stress components [σx, σy, τxy]
 */
void CalculateStressAt(double eta, double psi, double xe[4], double ye[4],
                      double E, double v, const double de[8], double* stress)
{
    // Calculate B matrix
    double B[8];
    CalculateStrainDisplacementMatrix(B, eta, psi, xe, ye);

    // Calculate stress components
    const double d33 = (1.0 - v) / 2.0;
    const double coef = E / (1.0 - v * v);

    // Normal stress σx
    stress[0] = coef * (B[0] * de[0] + B[1] * de[2] + B[2] * de[4] + B[3] * de[6] +
                v * (B[4] * de[1] + B[5] * de[3] + B[6] * de[5] + B[7] * de[7]));

    // Normal stress σy  
    stress[1] = coef * (B[4] * de[1] + B[5] * de[3] + B[6] * de[5] + B[7] * de[7] +
                v * (B[0] * de[0] + B[1] * de[2] + B[2] * de[4] + B[3] * de[6]));

    // Shear stress τxy
    stress[2] = coef * d33 * (B[4] * de[0] + B[0] * de[1] + B[5] * de[2] + B[1] * de[3] +
                              B[6] * de[4] + B[2] * de[5] + B[7] * de[6] + B[3] * de[7]);
}


/**
 * Calculate stress at element's Gauss points
 * @param stress[12] Array to store stress components at 4 Gauss points
 * @param Displacement Global displacement vector
 */
void CQ4::ElementStress(double stress[12], double* Displacement)
{
    // Convert coordinates from 3D to 2D local system
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // Get nodal displacements from global vector 
    double d[12];
    for (unsigned index = 0; index < 12; ++index)
    {
        d[index] = LocationMatrix_[index] ? Displacement[LocationMatrix_[index] - 1] : 0.0;
    }

    // Convert nodal displacements to 2D local system
    double de[8] = {
        d[0] * i[0] + d[1] * i[1] + d[2] * i[2],   d[0] * j[0] + d[1] * j[1] + d[2] * j[2],
        d[3] * i[0] + d[4] * i[1] + d[5] * i[2],   d[3] * j[0] + d[4] * j[1] + d[5] * j[2],
        d[6] * i[0] + d[7] * i[1] + d[8] * i[2],   d[6] * j[0] + d[7] * j[1] + d[8] * j[2],
        d[9] * i[0] + d[10] * i[1] + d[11] * i[2], d[9] * j[0] + d[10] * j[1] + d[11] * j[2]
    };

    // Get material properties
    const CQ4Material* material = static_cast<CQ4Material*>(ElementMaterial_);
    const double& E = material->E;  // Young's modulus
    const double& v = material->nu; // Poisson's ratio

    // Setup Gauss points
    const double pos = 1.0 / std::sqrt(3.0);
    const double etas[2] = {-pos, pos};
    const double psis[2] = {-pos, pos};

    // Point order (counter-clockwise):
    const int order[4][2] = {
        {0, 0},  // Lower left  (-pos,-pos)
        {1, 0},  // Lower right (+pos,-pos)
        {1, 1},  // Upper right (+pos,+pos)
        {0, 1}   // Upper left  (-pos,+pos) 
    };

    for(int gp = 0; gp < 4; gp++) {
        CalculateStressAt(etas[order[gp][0]], psis[order[gp][1]], xe, ye, E, v, de, stress + gp * 3);
    }

    /*
    // Calculate stresses at each Gauss point
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            const int index = (i * 2 + j) * 3;  // Each point has 3 stress components
            CalculateStressAt(etas[i], psis[j], xe, ye, E, v, de, stress + index);
        }
    }
    */
}

/**
 * Calculate coordinates of Gauss points in global system
 * @param gaussCoords Output array for Gauss point coordinates (size = 12)
 */
void CQ4::CalculateGaussPointCoordinates(double* gaussCoords)
{
    // Setup local coordinate system
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // Gauss point natural coordinates
    const double pos = 1.0 / std::sqrt(3.0);
    const double eta[2] = {-pos, pos};
    const double psi[2] = {-pos, pos};

    // Point order (counter-clockwise):
    const int order[4][2] = {
        {0, 0},  // Lower left  (-pos,-pos)
        {1, 0},  // Lower right (+pos,-pos)
        {1, 1},  // Upper right (+pos,+pos)
        {0, 1}   // Upper left  (-pos,+pos) 
    };

    int index = 0;
    for (int gp = 0; gp < 4; ++gp)
    {
        const int e = order[gp][0];
        const int p = order[gp][1];

        // Calculate shape functions
        const double N1 = 0.25 * (1 - eta[e]) * (1 - psi[p]);
        const double N2 = 0.25 * (1 + eta[e]) * (1 - psi[p]);
        const double N3 = 0.25 * (1 + eta[e]) * (1 + psi[p]);
        const double N4 = 0.25 * (1 - eta[e]) * (1 + psi[p]);

        // Calculate local coordinates
        const double x_local = N1 * xe[0] + N2 * xe[1] + N3 * xe[2] + N4 * xe[3];
        const double y_local = N1 * ye[0] + N2 * ye[1] + N3 * ye[2] + N4 * ye[3];

        // Transform to global coordinates
        gaussCoords[index++] = x_local * i[0] + y_local * j[0];  // X
        gaussCoords[index++] = x_local * i[1] + y_local * j[1];  // Y
        gaussCoords[index++] = x_local * i[2] + y_local * j[2];  // Z
    }
}

/**
 * Calculate displacements at Gauss points in global system
 * @param gaussDisp Output array for Gauss point displacements (size = 12)
 * @param Displacement Global displacement vector
 */
void CQ4::CalculateGaussPointDisplacement(double* gaussDisp, double* Displacement)
{
    // Setup coordinate systems and get local displacements
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // Get nodal displacements
    double d[12];
    for (unsigned index = 0; index < 12; ++index)
    {
        d[index] = LocationMatrix_[index] ? Displacement[LocationMatrix_[index] - 1] : 0.0;
    }

    // Convert displacements to local system
    double de[8] = {
        d[0] * i[0] + d[1] * i[1] + d[2] * i[2],   d[0] * j[0] + d[1] * j[1] + d[2] * j[2],
        d[3] * i[0] + d[4] * i[1] + d[5] * i[2],   d[3] * j[0] + d[4] * j[1] + d[5] * j[2],
        d[6] * i[0] + d[7] * i[1] + d[8] * i[2],   d[6] * j[0] + d[7] * j[1] + d[8] * j[2],
        d[9] * i[0] + d[10] * i[1] + d[11] * i[2], d[9] * j[0] + d[10] * j[1] + d[11] * j[2]
    };

    // Setup Gauss point coordinates
    const double pos = 1.0 / std::sqrt(3.0);
    const double eta[2] = {-pos, pos};
    const double psi[2] = {-pos, pos};
    
    // Point order (counter-clockwise):
    const int order[4][2] = {
        {0, 0},  // Lower left  (-pos,-pos)
        {1, 0},  // Lower right (+pos,-pos)
        {1, 1},  // Upper right (+pos,+pos)
        {0, 1}   // Upper left  (-pos,+pos) 
    };

    int index = 0;
    for (int gp = 0; gp < 4; ++gp)
    {
        const int e = order[gp][0];
        const int p = order[gp][1];

        // Calculate shape functions
        const double N1 = 0.25 * (1 - eta[e]) * (1 - psi[p]);
        const double N2 = 0.25 * (1 + eta[e]) * (1 - psi[p]);
        const double N3 = 0.25 * (1 + eta[e]) * (1 + psi[p]);
        const double N4 = 0.25 * (1 - eta[e]) * (1 + psi[p]);

        // Calculate local displacements
        const double ux_local = N1 * de[0] + N2 * de[2] + N3 * de[4] + N4 * de[6];
        const double uy_local = N1 * de[1] + N2 * de[3] + N3 * de[5] + N4 * de[7];

        // Transform to global system
        gaussDisp[index++] = ux_local * i[0] + uy_local * j[0];  // X
        gaussDisp[index++] = ux_local * i[1] + uy_local * j[1];  // Y
        gaussDisp[index++] = ux_local * i[2] + uy_local * j[2];  // Z
    }
}