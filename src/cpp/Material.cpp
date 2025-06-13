/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

/*********************************/
// Q4 material class by Jinhao Xu
/*********************************/

bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}
/*************************************************/
/* T3 Material class by Yu Jing                  */
/*************************************************/
//	Read material data from stream Input
bool CT3Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

void CT3Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}

/*************************************************/
/* H8 Material class by Yu Jing                  */
/*************************************************/
//	Read material data from stream Input
bool CH8Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

void CH8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}