
//by hly




#include "Elements/Beam.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// constructor
CBeam::CBeam()
{
	NEN_ = 2;	
	nodes_ = new CNode * [NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}
	//	Desconstructor
	CBeam::~CBeam()
{	
}

bool CBeam::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

void CBeam::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

void CBeam::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate beam length
    double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

    //计算局部坐标系的矩阵需要的元素
    const CBeamMaterial& material =
        static_cast<CBeamMaterial&>(*ElementMaterial_); 
    double Iz = (material.a * material.a * material.a) * material.b / 12 -
        pow(material.a - material.t1 - material.t3, 3.0) *
        (material.b - material.t2 - material.t4) / 12;
    double Iy = (material.b * material.b * material.b) * material.a / 12 -
        pow(material.b - material.t2 - material.t4, 3.0) *
        (material.a - material.t1 - material.t3) / 12;
    double Ip = Iz + Iy;
    double nx = DX[0] / L;
    double ny = DX[1] / L;
    double nz = DX[2] / L;
    double vx = material.n1;
    double vy = material.n2;
    double vz = material.n3;
    double wx = ny * vz - nz * vy;
    double wy = nz * vx - nx * vz;
    double wz = nx * vy - ny * vx;

    
    Matrix[0] = material.E * ((material.a * material.b * nx * nx) / L + (12 * Iz * vx * vx) / (L * L * L) + (12 * Iy * wx * wx) / (L * L * L));
    Matrix[1] = material.E * ((material.a * material.b * ny * ny) / L + (12 * Iz * vy * vy) / (L * L * L) + (12 * Iy * wy * wy) / (L * L * L));
    Matrix[2] = material.E * ((material.a * material.b * nx * ny) / L + (12 * Iz * vx * vy) / (L * L * L) + (12 * Iy * wx * wy) / (L * L * L));
    Matrix[3] = material.E * ((material.a * material.b * nz * nz) / L + (12 * Iz * vz * vz) / (L * L * L) + (12 * Iy * wz * wz) / (L * L * L));
    Matrix[4] = material.E * ((material.a * material.b * ny * nz) / L + (12 * Iz * vy * vz) / (L * L * L) + (12 * Iy * wy * wz) / (L * L * L));
    Matrix[5] = material.E * ((material.a * material.b * nx * nz) / L + (12 * Iz * vx * vz) / (L * L * L) + (12 * Iy * wx * wz) / (L * L * L));

    Matrix[6] = material.E * ((Ip * nx * nx) / ((2 + 2 * material.nu) * L) + (4 * Iy * vx * vx) / L + (4 * Iz * wx * wx) / L);
    Matrix[7] = 6 * material.E * (Iz * (vy * wz - vz * wy) - Iy * (vz * wy - vy * wz)) / (L * L);
    Matrix[8] = 6 * material.E * (Iz * (vx * wz - vz * wx) - Iy * (vz * wx - vx * wz)) / (L * L);
    Matrix[9] = 6 * material.E * (Iz - Iy) * vx * wx / (L * L);

    Matrix[10] = material.E * ((Ip * ny * ny) / ((2 + 2 * material.nu) * L) + (4 * Iy * vy * vy) / L + (4 * Iz * wy * wy) / L);
    Matrix[11] = material.E * ((Ip * nx * ny) / ((2 + 2 * material.nu) * L) + (4 * Iy * vx * vy) / L + (4 * Iz * wx * wy) / L);
    Matrix[12] = 6 * material.E * (Iz * (vy * wz - vz * wy) - Iy * (vz * wy - vy * wz)) / (L * L);
    Matrix[13] = 6 * material.E * (Iz - Iy) * vy * wy / (L * L);
    Matrix[14] = 6 * material.E * (Iz * (vy * wx - vx * wy) - Iy * (vx * wy - vy * wx)) / (L * L);

    Matrix[15] = material.E * ((Ip * nz * nz) / ((2 + 2 * material.nu) * L) + (4 * Iy * vz * vz) / L + (4 * Iz * wz * wz) / L);
    Matrix[16] = material.E * ((Ip * ny * nz) / ((2 + 2 * material.nu) * L) + (4 * Iy * vy * vz) / L + (4 * Iz * wy * wz) / L);
    Matrix[17] = material.E * ((Ip * nx * nz) / ((2 + 2 * material.nu) * L) + (4 * Iy * vx * vz) / L + (4 * Iz * wx * wz) / L);
    Matrix[18] = 6 * material.E * (Iz - Iy) * vz * wz / (L * L);
    Matrix[19] = 6 * material.E * (Iz * (vz * wy - vy * wz) - Iy * (vy * wz - vz * wy)) / (L * L);
    Matrix[20] = 6 * material.E * (Iz * (vx * wz - vz * wx) - Iy * (vz * wx - vx * wz)) / (L * L);

  
    Matrix[21] = Matrix[0];
    Matrix[22] = -Matrix[7];
    Matrix[23] = -Matrix[8];
    Matrix[24] = -Matrix[9];
    Matrix[25] = -Matrix[5];
    Matrix[26] = -Matrix[2];
    Matrix[27] = -Matrix[0];
    Matrix[28] = Matrix[1];
    Matrix[29] = Matrix[2];
    Matrix[30] = -Matrix[12];
    Matrix[31] = -Matrix[13];
    Matrix[32] = -Matrix[14];
    Matrix[33] = -Matrix[4];
    Matrix[34] = -Matrix[1];
    Matrix[35] = -Matrix[2];
    Matrix[36] = Matrix[3];
    Matrix[37] = Matrix[4];
    Matrix[38] = Matrix[5];
    Matrix[39] = -Matrix[18];
    Matrix[40] = -Matrix[19];
    Matrix[41] = -Matrix[20];
    Matrix[42] = -Matrix[3];
    Matrix[43] = -Matrix[4];
    Matrix[44] = -Matrix[5];

    
    Matrix[45] = material.E * (-(material.a * material.b * nx * nx) / L + (2 * Iy * vx * vx) / L + (2 * Iz * wx * wx) / L + Ip * nx * nx / ((2 + 2 * material.nu) * L));
    Matrix[46] = 6 * material.E * (Iy * (vz * wy - vy * wz) - Iz * (vy * wz - vz * wy)) / (L * L);
    Matrix[47] = 6 * material.E * (Iy * (vz * wx - vx * wz) - Iz * (vx * wz - vz * wx)) / (L * L);
    Matrix[48] = 6 * material.E * (Iy - Iz) * vx * wx / (L * L);
    Matrix[49] = material.E * (-(material.a * material.b * nx * nz) / L + (2 * Iy * vx * vz) / L + (2 * Iz * wx * wz) / L + Ip * nx * nz / ((2 + 2 * material.nu) * L));
    Matrix[50] = material.E * (-(material.a * material.b * nx * ny) / L + (2 * Iy * vx * vy) / L + (2 * Iz * wx * wy) / L + Ip * nx * ny / ((2 + 2 * material.nu) * L));
    Matrix[51] = material.E * (-(material.a * material.b * nx * nx) / L + (2 * Iy * vx * vx) / L + (2 * Iz * wx * wx) / L + Ip * nx * nx / ((2 + 2 * material.nu) * L));

  
    Matrix[52] = Matrix[7];
    Matrix[53] = Matrix[8];
    Matrix[54] = Matrix[9];
    Matrix[55] = material.E * (-(material.a * material.b * ny * ny) / L + (2 * Iy * vy * vy) / L + (2 * Iz * wy * wy) / L + Ip * ny * ny / ((2 + 2 * material.nu) * L));
    Matrix[56] = material.E * (-(material.a * material.b * nx * ny) / L + (2 * Iy * vx * vy) / L + (2 * Iz * wx * wy) / L + Ip * nx * ny / ((2 + 2 * material.nu) * L));
    Matrix[57] = 6 * material.E * (Iy * (vx * wz - vz * wx) - Iz * (vz * wx - vx * wz)) / (L * L);
    Matrix[58] = 6 * material.E * (Iy - Iz) * vy * wy / (L * L);
    Matrix[59] = 6 * material.E * (Iy * (vx * wy - vy * wx) - Iz * (vy * wx - vx * wy)) / (L * L);
    Matrix[60] = material.E * (-(material.a * material.b * ny * nz) / L + (2 * Iy * vy * vz) / L + (2 * Iz * wy * wz) / L + Ip * ny * nz / ((2 + 2 * material.nu) * L));
    Matrix[61] = material.E * (-(material.a * material.b * ny * ny) / L + (2 * Iy * vy * vy) / L + (2 * Iz * wy * wy) / L + Ip * ny * ny / ((2 + 2 * material.nu) * L));
    Matrix[62] = Matrix[50];
    Matrix[63] = Matrix[12];
    Matrix[64] = Matrix[13];
    Matrix[65] = Matrix[14];

    Matrix[66] = material.E * (-(material.a * material.b * nz * nz) / L + (2 * Iy * vz * vz) / L + (2 * Iz * wz * wz) / L + Ip * nz * nz / ((2 + 2 * material.nu) * L));
    Matrix[67] = material.E * (-(material.a * material.b * ny * nz) / L + (2 * Iy * vy * vz) / L + (2 * Iz * wy * wz) / L + Ip * ny * nz / ((2 + 2 * material.nu) * L));
    Matrix[68] = material.E * (-(material.a * material.b * nx * nz) / L + (2 * Iy * vx * vz) / L + (2 * Iz * wx * wz) / L + Ip * nx * nz / ((2 + 2 * material.nu) * L));
    Matrix[69] = 6 * material.E * (Iy - Iz) * vz * wz / (L * L);
    Matrix[70] = 6 * material.E * (Iy * (vy * wz - vz * wy) - Iz * (vz * wy - vy * wz)) / (L * L);
    Matrix[71] = 6 * material.E * (Iy * (vx * wz - vz * wx) - Iz * (vz * wx - vx * wz)) / (L * L);
    Matrix[72] = material.E * (-(material.a * material.b * nz * nz) / L + (2 * Iy * vz * vz) / L + (2 * Iz * wz * wz) / L + Ip * nz * nz / ((2 + 2 * material.nu) * L));
    Matrix[73] = Matrix[60];
    Matrix[74] = Matrix[49];
    Matrix[75] = Matrix[18];
    Matrix[76] = Matrix[19];
    Matrix[77] = Matrix[20];
}

void CBeam::ElementStress(double* stress, double* Displacement)
{
    const CBeamMaterial& material =
        static_cast<CBeamMaterial&>(*ElementMaterial_); // Pointer to material of the element
    clear(stress, 3);
    double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
    double L = 0; //	 beam length

    for (unsigned int i = 0; i < 3; i++)
    {
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
    }

    L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

    double S[6];
    for (unsigned int i = 0; i < 3; i++)
    {
        S[i] = -DX[i] * DX[i] * material.E / (L * L * L);
        S[i + 3] = -S[i];
    }

    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (LocationMatrix_[i * 6 + j])
            {
                double a = S[i * 3 + j] * Displacement[LocationMatrix_[i * 6 + j] - 1];
                stress[j] += S[i * 3 + j] * Displacement[LocationMatrix_[i * 6 + j] - 1];
            }
        }
    }
}
