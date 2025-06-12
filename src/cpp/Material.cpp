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

/***********************************/
// Q9 material class by Chenxuan Xu//
/***********************************/
//	读取材料属性
bool CQ9Material::Read(ifstream& Input)
{
	Input >> nset;	// 材料属性编号

	Input >> E >> nu;	// 杨氏模量、泊松比

	return true;
}

//	将材料属性输出至数据流
void CQ9Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}
