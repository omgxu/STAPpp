/*****************************************************************************/
/* STAP++ : T3 element made by Yu Jing                                       */
/* The input file must be 2-dimensional T3 element data                      */
/*****************************************************************************/

#include "Elements/T3.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

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


//  Constructor
CT3::CT3()
{
    NEN_ = 3; // Each element has 3 nodes
    nodes_ = new CNode*[NEN_];

    // ND_ = 9;
    ND_ = 9;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//  Desconstructor
CT3::~CT3()
{
    delete[] nodes_;
    delete[] LocationMatrix_;
}

//  Read element data from stream Input
bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;       // Material property set number
    unsigned int N1, N2, N3; // node number

    Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = static_cast<CT3Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];

    return true;
}

//  Write element data to stream
void CT3::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber 
           << setw(9)  << nodes_[1]->NodeNumber 
           << setw(9)  << nodes_[2]->NodeNumber 
           << setw(12) << ElementMaterial_->nset << endl;
}

void CT3::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    
    const CNode& n1 = *nodes_[0];
    const CNode& n2 = *nodes_[1];
    const CNode& n3 = *nodes_[2];

    double p13[2] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1]};
    double p12[2] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1]};
    double Area = 0.5 * (p13[0] * p12[1] - p13[1] * p12[0]);
//check if the area is zero
    if (Area == 0) {
        cerr << "Error: Area of the triangle is zero. Check node coordinates." << endl;
        return;
    }

    double x32 = n2.XYZ[0] - n3.XYZ[0];
    double y23 = n3.XYZ[1] - n2.XYZ[1];
    double x13 = n3.XYZ[0] - n1.XYZ[0];
    double y31 = n1.XYZ[1] - n3.XYZ[1];
    double x21 = n1.XYZ[0] - n2.XYZ[0];
    double y21 = n1.XYZ[1] - n2.XYZ[1];

// Pointer to material of the element
    CT3Material& material =
        *static_cast<CT3Material*>(ElementMaterial_); 

    auto E = material.E;
    auto v = material.nu;
    double cofD = E / (1 - v * v);

    vector<vector<double>> D(3, vector<double>(3));
    D[0][0] = 1; D[0][1] = v; D[0][2] = 0;
    D[1][0] = v; D[1][1] = 1; D[1][2] = 0;
    D[2][0] = 0; D[2][1] = 0; D[2][2] = (1 - v) / 2;
    for (size_t i1 = 0; i1 < 3; ++i1)
        for (size_t j1 = 0; j1 < 3; ++j1)
            D[i1][j1]*= cofD;
    
    // element stiffness matrix
    double b1 = -y23;
    double c1 = -x32;
    double b2 = -y31;
    double c2 = -x13;
    double b3 =  y21;
    double c3 = -x21;
    vector<vector<double>> B(3, vector<double>(6));
    B[0][0] = b1; B[0][1] = 0; B[0][2] = b2; B[0][3] = 0; B[0][4] = b3; B[0][5] = 0;
    B[1][0] = 0;  B[1][1] = c1; B[1][2] = 0;  B[1][3] = c2; B[1][4] = 0;  B[1][5] = c3;
    B[2][0] = c1; B[2][1] = b1; B[2][2] = c2; B[2][3] = b2; B[2][4] = c3; B[2][5] = b3;
    for (size_t i1 = 0; i1 < 3; ++i1)
        for (size_t j1 = 0; j1 < 6; ++j1)
            B[i1][j1]/= 2* Area;

    vector<vector<double>> K(6, vector<double>(6));
    K = matmul(matmul(transpose(B), D), B);
    for (size_t i1 = 0; i1 < 6; ++i1)
        for (size_t j1 = 0; j1 < 6; ++j1)
            K[i1][j1] *= Area;
    
    vector<vector<double>> LM(9, vector<double>(6));
    LM = {
        {1, 0, 0, 0, 0, 0}, 
        {0, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 0, 0}, 
        {0, 0, 1, 0, 0, 0}, 
        {0, 0, 0, 1, 0, 0}, 
        {0, 0, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 0}, 
        {0, 0, 0, 0, 0, 1}, 
        {0, 0, 0, 0, 0, 0}, 
    };
    // Convert the 2D stiffness matrix to 3D
    vector<vector<double>> MAT(9, vector<double>(9));
    MAT = matmul(LM, matmul(K, transpose(LM)));
    for (size_t i1 = 0; i1 < 9; ++i1)
        for (size_t j1 = 0; j1 <= i1; ++j1)
            Matrix[i1 * (i1 + 1) / 2 + j1] = MAT[i1-j1][i1];
}

//  Calculate element stress
void CT3::ElementStress(double* stress, double* Displacement)
{
    /*TODO: 修改成二维情况*/
    // Nodal coordinates
    const CNode& n1 = *nodes_[0];
    const CNode& n2 = *nodes_[1];
    const CNode& n3 = *nodes_[2];

    // Edge vector
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                   p31[0] * p21[1] - p31[1] * p21[0]};
    // generate area and normalize n at the same time
    double Area = std::sqrt(dot(n, n))/2.0;
    n[0] /= Area;
    n[1] /= Area;
    n[2] /= Area;
    // i = normalized p21, so that y21 = 0
    double i[3] = {p21[0], p21[1], p21[2]};
    normalize(i);
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};

    // generate M here
    double x32 = dot(p32, i);
    double y23 = -dot(p32, j);
    double x13 = -dot(p31, i);
    double y31 = dot(p31, j);
    double x21 = dot(p21, i);

    // form d first.
    // d represent 3d displacements at boundary nodes.
    double d[9];
    for (unsigned int index = 0; index < 9; ++index)
    {
        if (LocationMatrix_[index])
            d[index] = Displacement[LocationMatrix_[index] - 1];
        else
            d[index] = 0;
    }

    double de[6] = {
        i[0] * d[0] + i[1] * d[1] + i[2] * d[2],
        j[0] * d[0] + j[1] * d[1] + j[2] * d[2], // node 1 (2d)
        i[0] * d[3] + i[1] * d[4] + i[2] * d[5],
        j[0] * d[3] + j[1] * d[4] + j[2] * d[5], // node 2 (2d)
        i[0] * d[6] + i[1] * d[7] + i[2] * d[8],
        j[0] * d[6] + j[1] * d[7] + j[2] * d[8] // node 3 (2d)
    };

    CT3Material* material_ = dynamic_cast<CT3Material*>(ElementMaterial_);

    double v = material_->nu;
    double cof = material_->E / (4 * Area * (1 - v * v));

    stress[0] =
        2 * cof * (y23 * de[0] + y31 * de[2] + v * (x32 * de[1] + x13 * de[3] + x21 * de[5]));
    stress[1] = 2 * (v * y23 * de[0] + x32 * de[1] + v * y31 * de[2] + x13 * de[3] + x21 * de[5]);
    stress[2] = (1 - v) * (x32 * de[0] + y23 * de[1] + x13 * de[2] + y31 * de[3] + x21 * de[4]);

// #ifdef _TEST_
//     for (unsigned index = 0; index < 3; index++)
//     {
//         weights[index] = Area / 3;
//     }
//     for (unsigned GPindex = 0; GPindex < 3; ++GPindex)
//     {
//         for (unsigned dof = 0; dof < 3; dof++)
//         {
//             unsigned index = GPindex * 3 + dof;
//             GaussPosition[index] = nodes[GPindex]->XYZ[dof] * 2. / 3. +
//                                    nodes[(GPindex + 1) % 3]->XYZ[dof] / 6. +
//                                    nodes[(GPindex + 2) % 3]->XYZ[dof] / 6.;
//             GaussDisplacements[index] =
//                 d[index] * 2. / 3. + d[(index + 3) % 9] / 6. + d[(index + 6) % 9] / 6.;
//         }
//     }
// #endif
}

#ifdef _VIB_
void CT3::ElementMass(double* mass) {
}
#endif

// void CT3::ElementPostInfo(double* stress, double* Displacement, double* PrePositions,
//                                 double* PostPositions)
// {
//     ElementStress(stress, Displacement, nullptr, nullptr, nullptr);
//     for (unsigned index = 0; index < 9; ++index)
//     {
//         if (LocationMatrix_[index])
//         {
//             PrePositions[index] = nodes_[index / 3]->XYZ[index % 3];
//             PostPositions[index] = PrePositions[index] + Displacement[LocationMatrix_[index] - 1];
//         }
//         else
//         {
//             PrePositions[index] = PostPositions[index] = nodes_[index / 3]->XYZ[index % 3];
//         }
//     }
// }