/*****************************************************************************/
/* STAP++ : MindlinPlate element made by Yu Jing                             */
/*****************************************************************************/

#include "Elements/MindlinPlate.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CMindlinPlate::CMindlinPlate()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;   // Dimension of the location matrix
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CMindlinPlate::~CMindlinPlate()
{
	delete[] nodes_;
    delete[] LocationMatrix_;
}

//	Read element data from stream Input
bool CMindlinPlate::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Four node number

	Input >> N1 >> N2 >> N3 >> N4 >>MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CMindlinPlate::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9)  << nodes_[1]->NodeNumber 
           << setw(9)  << nodes_[2]->NodeNumber
           << setw(9)  << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CMindlinPlate::ElementStiffness(double* Matrix)
{
    // clear the stiffness matrix to initialize
	clear(Matrix, SizeOfStiffnessMatrix());
    // initialize the element stiffness matrix
    vector<vector<double>> ke(12, vector<double>(12, 0.0));

    // get coodinates of elements nodes
	const CNode& n1 = *nodes_[0];
    const CNode& n2 = *nodes_[1];
    const CNode& n3 = *nodes_[2];
    const CNode& n4 = *nodes_[3];

    // set reduced integration
    int ngpb = 2; //  Gauss integration of Bb
    int ngps = 1; //  Gauss integration of Bs

    // get gauss points and weights
    std::vector<double> wb, gpb;
    std::vector<double> ws, gps;
    gauss(ngpb, wb, gpb); // Gauss points and weights for Bb
    gauss(ngps, ws, gps); // Gauss points and weights for Bs

    /* TODO*/
    // compute element bending stiffness matrix Bb
    for (int i = 0; i < ngpb; ++i) {
        for (int j = 0; j < ngpb; ++j) {
            double eta = gpb[i];
            double psi = gpb[j];

            // 获取形函数导数矩阵和雅可比行列式
            vector<vector<double>> Bb, Bs;
            double detJ;
            BmatMindlinPlate(eta, psi, Bb, Bs, detJ);  // 假设这个函数已经实现

            // 计算Bb^T * Db * Bb
            vector<vector<double>> temp = matmul(transpose(Bb), material_->Db);
            vector<vector<double>> contrib = matmul(temp, Bb);

            // 累加到刚度矩阵
            for (int m = 0; m < 12; ++m) {
                for (int n = 0; n < 12; ++n) {
                    ke[m][n] += wb[i] * wb[j] * detJ * contrib[m][n];
                }
            }
        }
    }

    // compute element shear stiffness matrix Bs
    for (int i = 0; i < ngps; ++i) {
        for (int j = 0; j < ngps; ++j) {
            double eta = gps[i];
            double psi = gps[j];

            // 获取形函数导数矩阵和雅可比行列式
            vector<vector<double>> Bb, Bs;
            double detJ;
            BmatMindlinPlate(eta, psi, Bb, Bs, detJ);

            // 计算Bs^T * Ds * Bs
            vector<vector<double>> temp = matmul(transpose(Bs), material_->Ds);
            vector<vector<double>> contrib = matmul(temp, Bs);

            // 累加到刚度矩阵
            for (int m = 0; m < 12; ++m) {
                for (int n = 0; n < 12; ++n) {
                    ke[m][n] += ws[i] * ws[j] * detJ * contrib[m][n];
                }
            }
        }
    }

    // Store Ke's upper triangle part in Marix
    int index = 0;
    for (int j = 0; j < 12; ++j) {
        for (int i = 0; i <= j; ++i) {
            Matrix[index++] = ke[i][j];
        }
    }
}

//	Calculate element stress 
void CMindlinPlate::ElementStress(double* stress, double* Displacement)
{
	CMindlinPlate* material_ = dynamic_cast<CMindlinPlate*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	// sigma=EBd -- S=EB
	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}

// Function to calculate Gauss points and weights
// ngp: number of Gauss points
// w: output vector for weights
// gp: output vector for Gauss points in the parent element domain
void gauss(int ngp, std::vector<double>& w, std::vector<double>& gp) {
    // 清空输出向量
    w.clear();
    gp.clear();
    
    switch(ngp) {
        case 1:
            gp = {0.0};
            w = {2.0};
            break;
            
        case 2:
            gp = {-0.57735027, 0.57735027};
            w = {1.0, 1.0};
            break;
            
        case 3:
            gp = {-0.7745966692, 0.7745966692, 0.0};
            w = {0.5555555556, 0.5555555556, 0.8888888889};
            break;
            
        default:
            throw std::invalid_argument(
                "The given number (ngp = " + std::to_string(ngp) + 
                ") of Gauss points is too large and not implemented"
            );
    }
}