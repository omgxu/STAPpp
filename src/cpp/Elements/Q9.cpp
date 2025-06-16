#include "Elements/Q9.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <matrix.h>
using namespace std;

//	Constructor
CQ9::CQ9()
{
    NEN_ = 9; // Each element has 4 nodes
    nodes_ = new CNode*[NEN_];

    ND_ = 27; // 12 DOF in total
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = NULL;
}

//	Desconstructor解构器
CQ9::~CQ9() {}

//	Read element data from stream Input
bool CQ9::Read(ifstream& Input,  CMaterial* MaterialSets, CNode* NodeList)
{
    
    for (unsigned n = 0; n < 9; n++)
    {
        unsigned int N; // node indexs
        Input >> N;
        nodes_[n] = NodeList + N - 1;
    }

    unsigned int MSet; // Material property set number
    Input >> MSet;
    ElementMaterial_ = dynamic_cast<CQ9Material*>(MaterialSets) + MSet - 1;

    return true;
}

//	Write element data to stream
void CQ9::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber; // node indexs
    for (unsigned i = 1; i < 9; i++)
    {
        output << setw(9) << nodes_[i]->NodeNumber;
    }
    output << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CQ9::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
unsigned int CQ9::SizeOfStiffnessMatrix() { return ND_ * (ND_ + 1) / 2; }

// returns |Je|
// generate B 算Je（P170内容）
double GenerateB9Q(Matrix<double>& BB, const double xi, const double eta, const double xe[9],
                   const double ye[9])
{
    double GNData[18] = {((-1 + eta) * eta * (-1 + 2 * xi)) / 4.,
                         ((-1 + 2 * eta) * (-1 + xi) * xi) / 4.,
                         ((-1 + eta) * eta * (1 + 2 * xi)) / 4.,
                         ((-1 + 2 * eta) * xi * (1 + xi)) / 4.,
                         (eta * (1 + eta) * (1 + 2 * xi)) / 4.,
                         ((1 + 2 * eta) * xi * (1 + xi)) / 4.,
                         (eta * (1 + eta) * (-1 + 2 * xi)) / 4.,
                         ((1 + 2 * eta) * (-1 + xi) * xi) / 4.,
                         -((-1 + eta) * eta * xi),
                         -((-1 + 2 * eta) * (-1 + std::pow(xi, 2))) / 2.,
                         -((-1 + std::pow(eta, 2)) * (1 + 2 * xi)) / 2.,
                         -(eta * xi * (1 + xi)),
                         -(eta * (1 + eta) * xi),
                         -((1 + 2 * eta) * (-1 + std::pow(xi, 2))) / 2.,
                         -((-1 + std::pow(eta, 2)) * (-1 + 2 * xi)) / 2.,
                         -(eta * (-1 + xi) * xi),
                         2 * (-1 + std::pow(eta, 2)) * xi,
                         2 * eta * (-1 + std::pow(xi, 2))};

    Matrix<double> GN(2, 9, GNData);

    double dData[18];
    for (unsigned i = 0; i < 9; i++)
    {
        dData[i] = xe[i];
        dData[i + 9] = ye[i];
    }
    Matrix<double> d(9, 2, dData);

    // Je = GN * [xe ye]
    Matrix<double> Je = GN * d;

    BB = Je.inverse() * GN;

    double JeDet = Je.c_at(1, 1) * Je.c_at(2, 2) - Je.c_at(1, 2) * Je.c_at(2, 1);

    return JeDet;
}
// 求B（单元应变矩阵）（P180内容）
void AccumulateXiEta9Q(const double& xi, const double& eta, double weight, const double* xe,
                       const double* ye, Matrix<double>& ke, const double E, const double v)
{
    Matrix<double> BB(2, 9);
    double DetJe = GenerateB9Q(BB, xi, eta, xe, ye);

    Matrix<double> B(3, 18);
    for (unsigned i = 0; i < 9; i++)
    {
        double partialX = BB.c_at(1, 1 + i);
        double partialY = BB.c_at(2, 1 + i);
        B.at(1, 1 + 2 * i) = partialX;
        B.at(3, 1 + 2 * i + 1) = partialX;
        B.at(2, 1 + 2 * i + 1) = partialY;
        B.at(3, 1 + 2 * i) = partialY;
    }

    double DData[9] = {1, v, 0, v, 1, 0, 0, 0, (1. - v) / 2.};
    Matrix<double> D(3, 3, DData);
    double cof = E / (1 - v * v) * std::abs(DetJe) * weight;
    D = D * cof;

    Matrix<double>&& ke_ = B.transpose() * D * B;
    ke += ke_;
}

// convert ke' to ke with R (input as i and j)
void Convert2d23d9Q(Matrix<double>& k, double* matrix, const double i[3], const double j[3])
{
    auto f = [i, j](Matrix<double>::Pos_t row, Matrix<double>::Pos_t column) -> double {
        if ((row - 1) / 2 == (column - 1) / 3)
        {
            return ((row - 1) % 2 ? j : i)[(column - 1) % 3];
        }
        else
        {
            return double(0);
        }
    };
    Matrix<double> R(2 * 9, 3 * 9, f);
    Matrix<double> K = R.transpose() * k * R;
    for (unsigned column = 1; column <= 27; column++)
    {
        for (unsigned row = column; row >= 1; row--)
        {
            unsigned index = column * (column - 1) / 2 + (column - row);
            matrix[index] = K.c_at(row, column);
        }
    }
}

// calculate n, i, j and xe, ye
void Convert3d22d9Q(CNode* const nodes[9], double n[3], double i[3], double j[3], double xe[9],
                    double ye[9])
{
    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];

    // make p31 p21
    const double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    const double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};

    // n = p21 cross p31 (normalized)
    n[0] = -p31[1] * p21[2] + p31[2] * p21[1];
    n[1] = -p31[2] * p21[0] + p31[0] * p21[2];
    n[2] = -p31[0] * p21[1] + p31[1] * p21[0];
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

    // generate xe, ye
    for (unsigned n = 0; n < 9; n++)
    {
        xe[n] = i[0] * nodes[n]->XYZ[0] + i[1] * nodes[n]->XYZ[1] + i[2] * nodes[n]->XYZ[2];
        ye[n] = j[0] * nodes[n]->XYZ[0] + j[1] * nodes[n]->XYZ[1] + j[2] * nodes[n]->XYZ[2];
    }
}

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ9::ElementStiffness(double* matrix)
{
    clear(matrix, SizeOfStiffnessMatrix());

    // =========================== 3d to 2d ============================
    double n[3], i[3], j[3], xe[9], ye[9];
    Convert3d22d9Q(nodes_, n, i, j, xe, ye);

    // =========================== assembly Ke' =========================
    // generate GN4Q for eta, psi

    const double pos = std::sqrt(0.6f);
    const double xis[3] = {-pos, 0, pos};
    const double etas[3] = {-pos, 0, pos};
    const double weights[3] = {5. / 9., 8. / 9., 5. / 9.};

    const CQ9Material* material =
        dynamic_cast<CQ9Material*>(ElementMaterial_); // Pointer to material of the element
    const double& E = material->E;
    const double& v = material->nu;

    Matrix<double> ke(18, 18);
    for (unsigned px = 0; px < 3; px++)
        for (unsigned py = 0; py < 3; py++)
            AccumulateXiEta9Q(xis[px], etas[py], weights[px] * weights[py], xe, ye, ke, E, v);

    // ======================== assembly Ke (2d to 3d) ======================

    Convert2d23d9Q(ke, matrix, i, j);
    return;
}

void CalculateStressAt9Q(double xi, double eta, double xe[9], double ye[9], double E, double v,
                         const double de[18], double* stress)
{
    Matrix<double> BB(2, 9);
    GenerateB9Q(BB, xi, eta, xe, ye);
    Matrix<double> B(3, 18);
    for (unsigned i = 0; i < 9; i++)
    {
        double partialX = BB.c_at(1, 1 + i);
        double partialY = BB.c_at(2, 1 + i);
        B.at(1, 1 + 2 * i) = partialX;
        B.at(3, 1 + 2 * i + 1) = partialX;
        B.at(2, 1 + 2 * i + 1) = partialY;
        B.at(3, 1 + 2 * i) = partialY;
    }

    // sigma = D B d
    double DData[9] = {1, v, 0, v, 1, 0, 0, 0, (1. - v) / 2.};
    Matrix<double> D(3, 3, DData);

    Matrix<double> d(18, 1, de);

    Matrix<double>&& stress_ = D * B * d;

    stress[0] = stress_.c_at(1, 1);
    stress[1] = stress_.c_at(2, 1);
    stress[2] = stress_.c_at(3, 1);
}

void CalculateN9Q(double xi, double eta, double N[9])
{
    N[0] = ((-1 + eta) * eta * (-1 + xi) * xi) / 4.;
    N[1] = ((-1 + eta) * eta * xi * (1 + xi)) / 4.;
    N[2] = (eta * (1 + eta) * xi * (1 + xi)) / 4.;
    N[3] = (eta * (1 + eta) * (-1 + xi) * xi) / 4.;
    N[4] = ((-1 + eta) * eta * (1 - std::pow(xi, 2))) / 2.;
    N[5] = ((1 - std::pow(eta, 2)) * xi * (1 + xi)) / 2.;
    N[6] = (eta * (1 + eta) * (1 - std::pow(xi, 2))) / 2.;
    N[7] = ((1 - std::pow(eta, 2)) * (-1 + xi) * xi) / 2.;
    N[8] = (1 - std::pow(eta, 2)) * (1 - std::pow(xi, 2));
}

// generate 3d position and return weight
void CalculatePositionAt9Q(double eta, double psi, double xe[9], double ye[9], double i[3],
                           double j[3], double Positions[3])
{
    double N[9];
    CalculateN9Q(eta, psi, N);
    // generate local x
    double x2d = 0, y2d = 0;
    for (unsigned _ = 0; _ < 9; _++)
    {
        x2d += N[_] * xe[_];
        y2d += N[_] * ye[_];
    }

    // convert to 3d
    Positions[0] = i[0] * x2d + j[0] * y2d;
    Positions[1] = i[1] * x2d + j[1] * y2d;
    Positions[2] = i[2] * x2d + j[2] * y2d;
}

void CalculateDisplacementAt9Q(double eta, double psi, double de[8], double i[3], double j[3],
                               double Displacements[3])
{
    double N[4];
    CalculateN9Q(eta, psi, N);
    double ux = N[0] * de[0] + N[1] * de[2] + N[2] * de[4] + N[3] * de[6];
    double uy = N[0] * de[1] + N[1] * de[3] + N[2] * de[5] + N[3] * de[7];
    Displacements[0] = i[0] * ux + j[0] * uy;
    Displacements[1] = i[1] * ux + j[1] * uy;
    Displacements[2] = i[2] * ux + j[2] * uy;
}

double CalculateWeightAt9Q(double eta, double psi, double xe[4], double ye[4])
{
    
    return 0.0f;
}

void CQ9::ElementStress(double stress[27], double* Displacement)
{
    double Positions[27]; 
    // =========================== 3d to 2d ============================
    double n[3], i[3], j[3], xe[9], ye[9];
    Convert3d22d9Q(nodes_, n, i, j, xe, ye);
    double d[27];
    for (unsigned index = 0; index < 27; ++index)
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
    double de[2 * 9];
    clear(de, 2 * 9);
    for (unsigned nodeIndex = 0; nodeIndex < 9; nodeIndex++)
    {
        for (unsigned _ = 0; _ < 3; _++)
        {
            de[nodeIndex * 2 + 0] += d[3 * nodeIndex + _] * i[_];
            de[nodeIndex * 2 + 1] += d[3 * nodeIndex + _] * j[_];
        }
    }

    // ======================= calculate stress ========================
    CQ9Material* material =
        dynamic_cast<CQ9Material*>(ElementMaterial_); // Pointer to material of the element
    const double& E = material->E;
    const double& v = material->nu;

    const double pos = std::sqrt(0.6f);
    const double xis[3] = {-pos, 0, pos};
    const double etas[3] = {-pos, 0, pos};

    // calculate Positions
    for (unsigned px = 0; px < 3; px++)
    {

        for (unsigned py = 0; py < 3; py++)
        {
            CalculatePositionAt9Q(xis[px], etas[py], xe, ye, i, j, Positions + 9 * px + 3 * py);
            CalculateStressAt9Q(xis[px], etas[py], xe, ye, E, v, de, stress + 9 * px + 3 * py);
        }
    }
}

void CQ9::CalculateGaussPointDisplacement(double* gaussDisp, double* Displacement)
{
    double pos = std::sqrt(3.0) / std::sqrt(5.0);
    double eta[3] = {-pos, 0, pos};
    double psi[3] = {-pos, 0, pos};

    // 转换节点位移到局部坐标系
    double n[3], i[3], j[3], xe[9], ye[9];
    Convert3d22d9Q(nodes_, n, i, j, xe, ye);

    // 获取节点位移
    double d[27];
    for (unsigned index = 0; index < 27; ++index)
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
    double de[18] = {
        d[0] * i[0] + d[1] * i[1] + d[2] * i[2],   d[0] * j[0] + d[1] * j[1] + d[2] * j[2],
        d[3] * i[0] + d[4] * i[1] + d[5] * i[2],   d[3] * j[0] + d[4] * j[1] + d[5] * j[2],
        d[6] * i[0] + d[7] * i[1] + d[8] * i[2],   d[6] * j[0] + d[7] * j[1] + d[8] * j[2],
        d[9] * i[0] + d[10] * i[1] + d[11] * i[2], d[9] * j[0] + d[10] * j[1] + d[11] * j[2],
        d[12] * i[0] + d[13] * i[1] + d[14] * i[2], d[12] * j[0] + d[13] * j[1] + d[14] * j[2],
        d[15] * i[0] + d[16] * i[1] + d[17] * i[2], d[15] * j[0] + d[16] * j[1] + d[17] * j[2],
        d[18] * i[0] + d[19] * i[1] + d[20] * i[2], d[18] * j[0] + d[19] * j[1] + d[20] * j[2],
        d[21] * i[0] + d[22] * i[1] + d[23] * i[2], d[21] * j[0] + d[22] * j[1] + d[23] * j[2],
        d[24] * i[0] + d[25] * i[1] + d[26] * i[2], d[24] * j[0] + d[25] * j[1] + d[26] * j[2],
    };

    // 计算9个高斯点的位移
    int index = 0;
    for (int px = 0; px < 3; px++)
    {
        for (int py = 0; py < 3; py++)
        {
            double N[9];
            // 计算形函数值
            CalculateN9Q(eta[px], psi[py], N);

            // 计算局部坐标系下的位移
            double ux_local = 0.0, uy_local = 0.0;
            for (int i = 0; i < 9; i++)
            {
                ux_local += N[i] * de[2*i];
                uy_local += N[i] * de[2*i + 1];
            }

            // 转换到全局坐标系
            gaussDisp[index++] = ux_local * i[0] + uy_local * j[0];
            gaussDisp[index++] = ux_local * i[1] + uy_local * j[1];
            gaussDisp[index++] = ux_local * i[2] + uy_local * j[2];
        }
    }
}
