#include "Elements/Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ4::~CQ4()
{
}

//	Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// node indexs

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = NodeList + N1 - 1;
	nodes_[1] = NodeList + N2 - 1;
    nodes_[2] = NodeList + N3 - 1;
    nodes_[3] = NodeList + N4 - 1;

	return true;
}

//	Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset
           << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//  which is 77 BTW
unsigned int CQ4::SizeOfStiffnessMatrix() { return 12 * (12 + 1) / 2; }

// returns |Je|
// generate B
double GenerateB(double B[8], const double eta, const double psi, const double xe[4],
                 const double ye[4])
{
    double GN4Q[8] = {
        (eta - 1) / 4, (1 - eta) / 4,  (1 + eta) / 4, (-eta - 1) / 4, // first row
        (psi - 1) / 4, (-psi - 1) / 4, (1 + psi) / 4, (1 - psi) / 4   // second row
    };
    // Je = GN4Q * [xe ye]
    double Je[4] = {
        GN4Q[0] * xe[0] + GN4Q[1] * xe[1] + GN4Q[2] * xe[2] + GN4Q[3] * xe[3],
        GN4Q[0] * ye[0] + GN4Q[1] * ye[1] + GN4Q[2] * ye[2] + GN4Q[3] * ye[3], // first row
        GN4Q[4] * xe[0] + GN4Q[5] * xe[1] + GN4Q[6] * xe[2] + GN4Q[7] * xe[3],
        GN4Q[4] * ye[0] + GN4Q[5] * ye[1] + GN4Q[6] * ye[2] + GN4Q[7] * ye[3] // second row
    };
    double DetJe = Je[0] * Je[3] - Je[1] * Je[2];
    double JeI[4] = {
        Je[3] / DetJe, -Je[1] / DetJe, // first row
        -Je[2] / DetJe, Je[0] / DetJe  // second row
    };

    // generate [B] = Je^-1 * GN
    B[0] = JeI[0] * GN4Q[0] + JeI[1] * GN4Q[4]; // Be[1, 1]
    B[1] = JeI[0] * GN4Q[1] + JeI[1] * GN4Q[5]; // Be[1, 2]
    B[2] = JeI[0] * GN4Q[2] + JeI[1] * GN4Q[6]; // Be[1, 3]
    B[3] = JeI[0] * GN4Q[3] + JeI[1] * GN4Q[7]; // first row end
    B[4] = JeI[2] * GN4Q[0] + JeI[3] * GN4Q[4]; // Be[2, 1]
    B[5] = JeI[2] * GN4Q[1] + JeI[3] * GN4Q[5]; // Be[2, 2]
    B[6] = JeI[2] * GN4Q[2] + JeI[3] * GN4Q[6]; // Be[2, 3]
    B[7] = JeI[2] * GN4Q[3] + JeI[3] * GN4Q[7]; // second row

    return DetJe;
}

void AccumulateEtaPsi(const double& eta, const double& psi, const double& weight, const double* xe,
                      const double* ye, double* ke, const double E, const double v)
{
    double B[8];
    double DetJe = GenerateB(B, eta, psi, xe, ye);
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

// convert ke' to ke with R (input as i and j)
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

// calculate n, i, j and xe, ye
void Convert3d22d4Q(CNode* const nodes[4], double n[3], double i[3], double j[3], double xe[4],
                         double ye[4])
{
    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];
    const CNode& n4 = *nodes[3];

    // make p31 p21
    const double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    const double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    // double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    n[0] = p31[1] * p21[2] - p31[2] * p21[1];
    n[1] = p31[2] * p21[0] - p31[0] * p21[2];
    n[2] = p31[0] * p21[1] - p31[1] * p21[0];
    normalize(n);

    // i = normalized p21
    // i is manually set parallel to p21 so that y21 = 0
    i[0] = p21[0];
    i[1] = p21[1];
    i[2] = p21[2];
    normalize(i);
    // j = n cross i
    j[0] = n[1] * i[2] - n[2] * i[1];
    j[1] = n[2] * i[0] - n[0] * i[2];
    j[2] = n[0] * i[1] - n[1] * i[0];

    // by here, a conversion matrix is formed,
    // as (x', y') = ((i0, i1, i2), (j0, j1, j2)) . (x, y, z)

    // generate xe, ye
    xe[0] = i[0] * n1.XYZ[0] + i[1] * n1.XYZ[1] + i[2] * n1.XYZ[2];
    xe[1] = i[0] * n2.XYZ[0] + i[1] * n2.XYZ[1] + i[2] * n2.XYZ[2];
    xe[2] = i[0] * n3.XYZ[0] + i[1] * n3.XYZ[1] + i[2] * n3.XYZ[2];
    xe[3] = i[0] * n4.XYZ[0] + i[1] * n4.XYZ[1] + i[2] * n4.XYZ[2];

    ye[0] = j[0] * n1.XYZ[0] + j[1] * n1.XYZ[1] + j[2] * n1.XYZ[2];
    ye[1] = j[0] * n2.XYZ[0] + j[1] * n2.XYZ[1] + j[2] * n2.XYZ[2];
    ye[2] = j[0] * n3.XYZ[0] + j[1] * n3.XYZ[1] + j[2] * n3.XYZ[2];
    ye[3] = j[0] * n4.XYZ[0] + j[1] * n4.XYZ[1] + j[2] * n4.XYZ[2];
}

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    // =========================== 3d to 2d ============================
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // =========================== assembly Ke' =========================
    // generate GN4Q for eta, psi

    double ke[36];
    // for ke', size = (8, 8), number of elements = 8*9/2 = 36
    // notice ke is cached column by column, from up to down.
    clear(ke, 36);
    const double pos = 1 / std::sqrt(3.0f);
    const double etas[2] = {-pos, pos};
    const double psis[2] = {-pos, pos};
    const double weights[2][2] = {{1.0, 1.0}, {1.0, 1.0}};

    const CQ4Material* material =
        static_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
    const double& E = material->E;
    const double& v = material->nu;
    AccumulateEtaPsi(etas[0], psis[0], weights[0][0], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[0], psis[1], weights[0][1], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[1], psis[0], weights[1][0], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[1], psis[1], weights[1][1], xe, ye, ke, E, v);

    // ======================== assembly Ke (2d to 3d) ======================

    Convert2d23d(ke, Matrix, i, j);
    return;
}

void CalculateStressAt(double eta, double psi, double xe[4], double ye[4], double E, double v,
                       const double de[8], double* stress)
{
    // generate B first
    double B[8];
    GenerateB(B, eta, psi, xe, ye);

    // see ../../memo/4Q/4Q.nb and ../../memo/4Q/4Q-calc-stress.py
    double d33 = (1.f - v) / 2.0f;
    const double cof = E / (1 - v * v);
    stress[0] = cof * (B[0] * de[0] + B[1] * de[2] + B[2] * de[4] + B[3] * de[6] +
                       v * (B[4] * de[1] + B[5] * de[3] + B[6] * de[5] + B[7] * de[7]));
    stress[1] =
        cof * (B[4] * de[1] + B[5] * de[3] + B[6] * de[5] +
               v * (B[0] * de[0] + B[1] * de[2] + B[2] * de[4] + B[3] * de[6]) + B[7] * de[7]);
    stress[2] = cof * (d33 * (B[4] * de[0] + B[0] * de[1] + B[5] * de[2] + B[1] * de[3] +
                              B[6] * de[4] + B[2] * de[5] + B[7] * de[6] + B[3] * de[7]));
}


//	Calculate element stress
void CQ4::ElementStress(double stress[12], double* Displacement)
{
    // =========================== 3d to 2d ============================
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // form d first.
    // d represent 3d displacements at boundary nodes.
    double d[12];
    for (unsigned index = 0; index < 12; ++index)
    {
        if (LocationMatrix_[index])
        {
            d[index] = Displacement[LocationMatrix_[index] - 1];
        }
        else
        {
            d[index] = 0.0;
        }
    }

    // generate de, convert from 3d to 2d.
    double de[8] = {
        d[0] * i[0] + d[1] * i[1] + d[2] * i[2],   d[0] * j[0] + d[1] * j[1] + d[2] * j[2],
        d[3] * i[0] + d[4] * i[1] + d[5] * i[2],   d[3] * j[0] + d[4] * j[1] + d[5] * j[2],
        d[6] * i[0] + d[7] * i[1] + d[8] * i[2],   d[6] * j[0] + d[7] * j[1] + d[8] * j[2],
        d[9] * i[0] + d[10] * i[1] + d[11] * i[2], d[9] * j[0] + d[10] * j[1] + d[11] * j[2]};

    // ======================= calculate stress ========================
    CQ4Material* material =
        static_cast<CQ4Material*>(ElementMaterial_); // Pointer to material of the element
    const double& E = material->E;
    const double& v = material->nu;

    double pos = 1 / std::sqrt(3.0f);
    double etas[2] = {-pos, pos};
    double psis[2] = {-pos, pos};

    // calculate stresses
    CalculateStressAt(etas[0], psis[0], xe, ye, E, v, de, stress + 0);
    CalculateStressAt(etas[0], psis[1], xe, ye, E, v, de, stress + 3);
    CalculateStressAt(etas[1], psis[0], xe, ye, E, v, de, stress + 6);
    CalculateStressAt(etas[1], psis[1], xe, ye, E, v, de, stress + 9);
}

void CQ4::CalculateGaussPointCoordinates(double* gaussCoords)
{
    double pos = 1 / std::sqrt(3.0);
    double eta[2] = {-pos, pos};
    double psi[2] = {-pos, pos};

    // 转换节点坐标到局部坐标系
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // 遍历高斯点，标准顺序（逆时针）：
    // (-pos,-pos), (+pos,-pos), (-pos,+pos), (+pos,+pos)
    const int order[4][2] = {
        {0, 0},  // 第1个高斯点 (-pos,-pos)  左下
        {1, 0},  // 第2个高斯点 (+pos,-pos)  右下
        {0, 1},  // 第3个高斯点 (-pos,+pos)  左上
        {1, 1}   // 第4个高斯点 (+pos,+pos)  右上
    };

    int index = 0;
    for (int gp = 0; gp < 4; ++gp)
    {
        int e = order[gp][0];
        int p = order[gp][1];

        // 形函数
        double N1 = 0.25 * (1 - eta[e]) * (1 - psi[p]);
        double N2 = 0.25 * (1 + eta[e]) * (1 - psi[p]);
        double N3 = 0.25 * (1 + eta[e]) * (1 + psi[p]);
        double N4 = 0.25 * (1 - eta[e]) * (1 + psi[p]);

        // 高斯点在局部坐标系中的坐标
        double x_local = N1 * xe[0] + N2 * xe[1] + N3 * xe[2] + N4 * xe[3];
        double y_local = N1 * ye[0] + N2 * ye[1] + N3 * ye[2] + N4 * ye[3];

        // 转换到三维全局坐标系
        gaussCoords[index++] = x_local * i[0] + y_local * j[0];
        gaussCoords[index++] = x_local * i[1] + y_local * j[1];
        gaussCoords[index++] = x_local * i[2] + y_local * j[2];
    }
}

void CQ4::CalculateGaussPointDisplacement(double* gaussDisp, double* Displacement)
{
    double pos = 1 / std::sqrt(3.0);
    double eta[2] = {-pos, pos};
    double psi[2] = {-pos, pos};

    // 转换节点位移到局部坐标系
    double n[3], i[3], j[3], xe[4], ye[4];
    Convert3d22d4Q(nodes_, n, i, j, xe, ye);

    // 获取节点位移
    double d[12];
    for (unsigned index = 0; index < 12; ++index)
    {
        if (LocationMatrix_[index])
        {
            d[index] = Displacement[LocationMatrix_[index] - 1];
        }
        else
        {
            d[index] = 0.0;
        }
    }

    // 转换到局部坐标系的节点位移
    double de[8] = {
        d[0] * i[0] + d[1] * i[1] + d[2] * i[2],   d[0] * j[0] + d[1] * j[1] + d[2] * j[2],
        d[3] * i[0] + d[4] * i[1] + d[5] * i[2],   d[3] * j[0] + d[4] * j[1] + d[5] * j[2],
        d[6] * i[0] + d[7] * i[1] + d[8] * i[2],   d[6] * j[0] + d[7] * j[1] + d[8] * j[2],
        d[9] * i[0] + d[10] * i[1] + d[11] * i[2], d[9] * j[0] + d[10] * j[1] + d[11] * j[2]
    };

    // 遍历高斯点，标准顺序（逆时针）：
    // (-pos,-pos), (+pos,-pos), (-pos,+pos), (+pos,+pos)
    const int order[4][2] = {
        {0, 0},  // 第1个高斯点 (-pos,-pos)  左下
        {1, 0},  // 第2个高斯点 (+pos,-pos)  右下
        {0, 1},  // 第3个高斯点 (-pos,+pos)  左上
        {1, 1}   // 第4个高斯点 (+pos,+pos)  右上
    };

    int index = 0;
    for (int gp = 0; gp < 4; ++gp)
    {
        int e = order[gp][0];
        int p = order[gp][1];

        // 形函数
        double N1 = 0.25 * (1 - eta[e]) * (1 - psi[p]);
        double N2 = 0.25 * (1 + eta[e]) * (1 - psi[p]);
        double N3 = 0.25 * (1 + eta[e]) * (1 + psi[p]);
        double N4 = 0.25 * (1 - eta[e]) * (1 + psi[p]);

        // 高斯点在局部坐标系中的位移
        double ux_local = N1 * de[0] + N2 * de[2] + N3 * de[4] + N4 * de[6];
        double uy_local = N1 * de[1] + N2 * de[3] + N3 * de[5] + N4 * de[7];

        // 转换到三维全局坐标系
        gaussDisp[index++] = ux_local * i[0] + uy_local * j[0];
        gaussDisp[index++] = ux_local * i[1] + uy_local * j[1];
        gaussDisp[index++] = ux_local * i[2] + uy_local * j[2];
    }
}