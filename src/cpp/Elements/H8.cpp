/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/H8.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>

using namespace std;
const int ngp = 2; // 高斯积分点个数

//	Constructor
CH8::CH8()
{
	NEN_ = 8;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CH8::~CH8()
{
	delete[] nodes_;
    delete[] LocationMatrix_;
}

//	Read element data from stream Input
bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// Eight node numbers
    Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;

	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    nodes_[4] = &NodeList[N5 - 1];
    nodes_[5] = &NodeList[N6 - 1];
    nodes_[6] = &NodeList[N7 - 1];
    nodes_[7] = &NodeList[N8 - 1];

	return true;
}

//	Write element data to stream
void CH8::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(9) << nodes_[4]->NodeNumber
           << setw(9) << nodes_[5]->NodeNumber
           << setw(9) << nodes_[6]->NodeNumber
           << setw(9) << nodes_[7]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

// dot product function of two 3D vectors
inline double dot(const double* p1, const double* p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

// matrix multiplication
vector<vector<double>> matmul(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    size_t m = A.size(), n = A[0].size(), p = B[0].size();
    vector<vector<double>> C(m, vector<double>(p, 0.0));
    for (size_t i = 0; i < m; ++i)
        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < p; ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// matrix & vector multiplication
// function overloading
std::vector<double> matmul(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    size_t m = A.size(), n = A[0].size();
    std::vector<double> y(m, 0.0);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

// transpose matrix
vector<vector<double>> transpose(const vector<vector<double>>& A) {
    if (A.empty() || A[0].empty()) return {};

    size_t rows = A.size();
    size_t cols = A[0].size();
    vector<vector<double>> AT(cols, vector<double>(rows)); 
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            AT[j][i] = A[i][j];

    return AT;
}

#include <vector>
#include <stdexcept>

// 获取高斯点和权重
void gauss(int ngp, std::vector<double>& w, std::vector<double>& gp)
{
    w.clear();
    gp.clear();
    if (ngp == 1) {
        gp = {0.0};
        w = {2.0};
    } else if (ngp == 2) {
        gp = {-0.57735027, 0.57735027};
        w = {1.0, 1.0};
    } else if (ngp == 3) {
        gp = {-0.7745966692, 0.7745966692, 0.0};
        w = {0.5555555556, 0.5555555556, 0.8888888889};
    } else {
        throw std::runtime_error("The given number (ngp = " + std::to_string(ngp) + ") of Gauss points is too large and not implemented");
    }
}

// 计算H8单元的形函数矩阵N
std::vector<std::vector<double>> NmatElastT3(double zeta, double eta, double psi)
{
    double N1 = 0.125 * (1 - psi) * (1 - eta) * (1 - zeta);
    double N2 = 0.125 * (1 + psi) * (1 - eta) * (1 - zeta);
    double N3 = 0.125 * (1 + psi) * (1 + eta) * (1 - zeta);
    double N4 = 0.125 * (1 - psi) * (1 + eta) * (1 - zeta);
    double N5 = 0.125 * (1 - psi) * (1 - eta) * (1 + zeta);
    double N6 = 0.125 * (1 + psi) * (1 - eta) * (1 + zeta);
    double N7 = 0.125 * (1 + psi) * (1 + eta) * (1 + zeta);
    double N8 = 0.125 * (1 - psi) * (1 + eta) * (1 + zeta);

    std::vector<std::vector<double>> N(3, std::vector<double>(24, 0.0));
    N[0][0] = N1; N[1][1] = N1; N[2][2] = N1;
    N[0][3]  = N2; N[1][4]  = N2; N[2][5]  = N2;
    N[0][6]  = N3; N[1][7]  = N3; N[2][8]  = N3;
    N[0][9]  = N4; N[1][10] = N4; N[2][11] = N4;
    N[0][12] = N5; N[1][13] = N5; N[2][14] = N5;
    N[0][15] = N6; N[1][16] = N6; N[2][17] = N6;
    N[0][18] = N7; N[1][19] = N7; N[2][20] = N7;
    N[0][21] = N8; N[1][22] = N8; N[2][23] = N8;
    return N;
}

// 计算T3单元的B矩阵和Jacobi行列式
std::pair<std::vector<std::vector<double>>, double>
BmatElastH8(double zeta, double eta, double psi, const std::vector<std::vector<double>>& C)
{
    // 1. 计算8个节点的形函数对zeta, eta, psi的导数
    double dNdzeta[8], dNdeta[8], dNdpsi[8];
    dNdzeta[0] = -0.125 * (1 - psi) * (1 - eta);
    dNdzeta[1] = -0.125 * (1 + psi) * (1 - eta);
    dNdzeta[2] = -0.125 * (1 + psi) * (1 + eta);
    dNdzeta[3] = -0.125 * (1 - psi) * (1 + eta);
    dNdzeta[4] =  0.125 * (1 - psi) * (1 - eta);
    dNdzeta[5] =  0.125 * (1 + psi) * (1 - eta);
    dNdzeta[6] =  0.125 * (1 + psi) * (1 + eta);
    dNdzeta[7] =  0.125 * (1 - psi) * (1 + eta);

    dNdeta[0] = -0.125 * (1 - psi) * (1 - zeta);
    dNdeta[1] = -0.125 * (1 + psi) * (1 - zeta);
    dNdeta[2] =  0.125 * (1 + psi) * (1 - zeta);
    dNdeta[3] =  0.125 * (1 - psi) * (1 - zeta);
    dNdeta[4] = -0.125 * (1 - psi) * (1 + zeta);
    dNdeta[5] = -0.125 * (1 + psi) * (1 + zeta);
    dNdeta[6] =  0.125 * (1 + psi) * (1 + zeta);
    dNdeta[7] =  0.125 * (1 - psi) * (1 + zeta);

    dNdpsi[0] = -0.125 * (1 - eta) * (1 - zeta);
    dNdpsi[1] =  0.125 * (1 - eta) * (1 - zeta);
    dNdpsi[2] =  0.125 * (1 + eta) * (1 - zeta);
    dNdpsi[3] = -0.125 * (1 + eta) * (1 - zeta);
    dNdpsi[4] = -0.125 * (1 - eta) * (1 + zeta);
    dNdpsi[5] =  0.125 * (1 - eta) * (1 + zeta);
    dNdpsi[6] =  0.125 * (1 + eta) * (1 + zeta);
    dNdpsi[7] = -0.125 * (1 + eta) * (1 + zeta);

    // 2. 计算Jacobi矩阵 J (3x3)
    std::vector<std::vector<double>> J(3, std::vector<double>(3, 0.0));
    for (int i = 0; i < 8; ++i) {
        J[0][0] += dNdzeta[i] * C[i][0]; // dx/dzeta
        J[0][1] += dNdeta[i]  * C[i][0]; // dx/deta
        J[0][2] += dNdpsi[i]  * C[i][0]; // dx/dpsi

        J[1][0] += dNdzeta[i] * C[i][1]; // dy/dzeta
        J[1][1] += dNdeta[i]  * C[i][1]; // dy/deta
        J[1][2] += dNdpsi[i]  * C[i][1]; // dy/dpsi

        J[2][0] += dNdzeta[i] * C[i][2]; // dz/dzeta
        J[2][1] += dNdeta[i]  * C[i][2]; // dz/deta
        J[2][2] += dNdpsi[i]  * C[i][2]; // dz/dpsi
    }

    // 3. Jacobi行列式
    double detJ =
        J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    if (std::abs(detJ) < 1e-12)
        throw std::runtime_error("Jacobian determinant is zero!");

    // 4. Jacobi逆矩阵
    std::vector<std::vector<double>> Jinv(3, std::vector<double>(3, 0.0));
    double invdet = 1.0 / detJ;
    Jinv[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * invdet;
    Jinv[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]) * invdet;
    Jinv[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * invdet;
    Jinv[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]) * invdet;
    Jinv[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * invdet;
    Jinv[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]) * invdet;
    Jinv[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * invdet;
    Jinv[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]) * invdet;
    Jinv[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * invdet;

    // 5. 计算全局坐标下的形函数导数 dN/dx, dN/dy, dN/dz
    std::vector<std::vector<double>> dNdX(8, std::vector<double>(3, 0.0));
    for (int i = 0; i < 8; ++i) {
        // dN/dx = Jinv[0][0]*dNdzeta + Jinv[0][1]*dNdeta + Jinv[0][2]*dNdpsi
        dNdX[i][0] = Jinv[0][0]*dNdzeta[i] + Jinv[0][1]*dNdeta[i] + Jinv[0][2]*dNdpsi[i];
        // dN/dy = Jinv[1][0]*dNdzeta + Jinv[1][1]*dNdeta + Jinv[1][2]*dNdpsi
        dNdX[i][1] = Jinv[1][0]*dNdzeta[i] + Jinv[1][1]*dNdeta[i] + Jinv[1][2]*dNdpsi[i];
        // dN/dz = Jinv[2][0]*dNdzeta + Jinv[2][1]*dNdeta + Jinv[2][2]*dNdpsi[i];
        dNdX[i][2] = Jinv[2][0]*dNdzeta[i] + Jinv[2][1]*dNdeta[i] + Jinv[2][2]*dNdpsi[i];
    }

    // 6. 组装B矩阵（6×24）
    std::vector<std::vector<double>> B(6, std::vector<double>(24, 0.0));
    for (int i = 0; i < 8; ++i) {
        int col = i * 3;
        B[0][col    ] = dNdX[i][0]; // dN_i/dx
        B[1][col + 1] = dNdX[i][1]; // dN_i/dy
        B[2][col + 2] = dNdX[i][2]; // dN_i/dz
        B[3][col    ] = dNdX[i][1]; // dN_i/dy
        B[3][col + 1] = dNdX[i][0]; // dN_i/dx
        B[4][col + 1] = dNdX[i][2]; // dN_i/dz
        B[4][col + 2] = dNdX[i][1]; // dN_i/dy
        B[5][col    ] = dNdX[i][2]; // dN_i/dz
        B[5][col + 2] = dNdX[i][0]; // dN_i/dx
    }

    return std::make_pair(B, detJ);
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CH8::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

    const CNode& n1 = *nodes_[0];
    const CNode& n2 = *nodes_[1];
    const CNode& n3 = *nodes_[2];
    const CNode& n4 = *nodes_[3];
    const CNode& n5 = *nodes_[4];
    const CNode& n6 = *nodes_[5];
    const CNode& n7 = *nodes_[6];
    const CNode& n8 = *nodes_[7];

    std::vector<std::vector<double>> C = {
        {n1.XYZ[0], n1.XYZ[1], n1.XYZ[2]},
        {n2.XYZ[0], n2.XYZ[1], n2.XYZ[2]},
        {n3.XYZ[0], n3.XYZ[1], n3.XYZ[2]},
        {n4.XYZ[0], n4.XYZ[1], n4.XYZ[2]},
        {n5.XYZ[0], n5.XYZ[1], n5.XYZ[2]},
        {n6.XYZ[0], n6.XYZ[1], n6.XYZ[2]},
        {n7.XYZ[0], n7.XYZ[1], n7.XYZ[2]},
        {n8.XYZ[0], n8.XYZ[1], n8.XYZ[2]}
    };

    CH8Material& material =
        *static_cast<CH8Material*>(ElementMaterial_); 

    auto E = material.E;
    auto v = material.nu;
    double cofD = E / (1 - v * v);

    double cofD = E / ((1 + v) * (1 - 2 * v));
    vector<vector<double>> D(6, vector<double>(6, 0.0));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            D[i][j] = v;
    for (int i = 0; i < 3; ++i)
        D[i][i] = 1 - v;
    for (int i = 3; i < 6; ++i)
        D[i][i] = (1 - 2 * v) / 2.0;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            D[i][j] *= cofD;

    vector<vector<double>> K(6, vector<double>(6));

    std::vector<double> w, gp;
    gauss(ngp, w, gp); // 高斯积分点和权重

    for (int i = 0; i < ngp; ++i) {
        for (int j = 0; j < ngp; ++j) {
            for (int k = 0; k < ngp; ++k){
                double eta = gp[i];
                double psi = gp[j];
                double zeta = gp[k];
                auto [B, detJ]=BmatElastH8(zeta,eta,psi,C);
                auto KB = matmul(transpose(B), D);
                auto KBB = matmul(KB, B);
                for (size_t i = 0; i < K.size(); ++i)
                    for (size_t j = 0; j < K[0].size(); ++j)
                        K[i][j] += w[i]*w[j]*w[k]*detJ*KBB[i][j];
            }
        }
    
    for (int i1 = 0; i1 < 24; ++i1)
        for (int j1 = 0; j1 <= i1; ++j1)
            Matrix[i1 * (i1 + 1) / 2 + j1] = K[i1-j1][i1];
    }
}

//	Calculate element stress 
void CH8::ElementStress(double* stress, double* Displacement)
{
    // Nodal coordinates
    const CNode& n1 = *nodes_[0];
    const CNode& n2 = *nodes_[1];
    const CNode& n3 = *nodes_[2];
    const CNode& n4 = *nodes_[3];
    const CNode& n5 = *nodes_[4];
    const CNode& n6 = *nodes_[5];
    const CNode& n7 = *nodes_[6];
    const CNode& n8 = *nodes_[7];

    // Coordinate matrix C
    std::vector<std::vector<double>> C = {
        {n1.XYZ[0], n1.XYZ[1], n1.XYZ[2]},
        {n2.XYZ[0], n2.XYZ[1], n2.XYZ[2]},
        {n3.XYZ[0], n3.XYZ[1], n3.XYZ[2]},
        {n4.XYZ[0], n4.XYZ[1], n4.XYZ[2]},
        {n5.XYZ[0], n5.XYZ[1], n5.XYZ[2]},
        {n6.XYZ[0], n6.XYZ[1], n6.XYZ[2]},
        {n7.XYZ[0], n7.XYZ[1], n7.XYZ[2]},
        {n8.XYZ[0], n8.XYZ[1], n8.XYZ[2]}
    };

    CH8Material& material =*static_cast<CH8Material*>(ElementMaterial_); 

    // Elastic tensor matrix D
    auto E = material.E;
    auto v = material.nu;
    double cofD = E / (1 - v * v);

    double cofD = E / ((1 + v) * (1 - 2 * v));
    vector<vector<double>> D(6, vector<double>(6, 0.0));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            D[i][j] = v;
    for (int i = 0; i < 3; ++i)
        D[i][i] = 1 - v;
    for (int i = 3; i < 6; ++i)
        D[i][i] = (1 - 2 * v) / 2.0;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            D[i][j] *= cofD;

    // Calculate the stress of 8 Gaussian points in sequence
    std::vector<double> w, gp;
    gauss(ngp, w, gp); // 高斯积分点和权重
    std::vector<double> d(Displacement, Displacement + 24);

    for (int i = 0; i < ngp; ++i) {
        for (int j = 0; j < ngp; ++j) {
            for (int k = 0; k < ngp; ++k){
            double eta = gp[i];
            double psi = gp[j];
            double zeta = gp[k];
            auto [B, detJ]=BmatElastH8(zeta,eta,psi,C);
            auto epsilon = matmul(B, d);
            auto sstress = matmul(D, epsilon);
            stress[(4*i+2*j+k)*6] = sstress[0];
            stress[(4*i+2*j+k)*6+1] = sstress[1];
            stress[(4*i+2*j+k)*6+2] = sstress[2];
            stress[(4*i+2*j+k)*6+3] = sstress[3];
            stress[(4*i+2*j+k)*6+4] = sstress[4];
            stress[(4*i+2*j+k)*6+5] = sstress[5];
            }
        }
    }
}
