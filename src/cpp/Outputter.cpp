/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::Q4: // Q4 element
				OutputQ4Elements(EleGrp);
			case ElementTypes::T3: // 3T element
				OutputT3Elements(EleGrp);
			case ElementTypes::Beam:
				OutputBeamElements(EleGrp);
			case ElementTypes::Q9: // Q9 element
				OutputQ9Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}

//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}
void COutputter::OutputBeamElements(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::GetInstance();

    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    unsigned int NUMMAT = ElementGroup.GetNUMMAT();

    *this << " M A T E R I A L   D E F I N I T I O N" << endl
          << endl;
    *this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
    *this << " AND SECTION PROPERTIES  . . . . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
          << endl
          << endl;

    *this << "  SET       YOUNG'S     POISSON     SECTION      SECTION      " << endl
          << " NUMBER     MODULUS      RATIO         a            b         " << endl
          << "               E           nu" << endl;

    *this << setiosflags(ios::scientific) << setprecision(5);

    // ������в�������
    for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset + 1;
        ElementGroup.GetMaterial(mset).Write(*this);
    }

    *this << endl << endl
          << " E L E M E N T   I N F O R M A T I O N" << endl;

    *this << " ELEMENT     NODE     NODE       MATERIAL" << endl
          << " NUMBER-N      I        J       SET NUMBER" << endl;

    unsigned int NUME = ElementGroup.GetNUME();

    // ������е�Ԫ��Ϣ
    for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele + 1;
        ElementGroup[Ele].Write(*this);
    }

    *this << endl;
}


// 	Print Q4 element data
void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();
	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl
		  << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;
	*this << "  SET       YOUNG'S     POISSON'S" << endl
		  << " NUMBER     MODULUS        RATIO" << endl
		  << "               E              NU" << endl;
	*this << setiosflags(ios::scientific) << setprecision(5);
	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}
	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE       NODE       NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J         K         L       SET NUMBER" << endl;
	unsigned int NUME = ElementGroup.GetNUME();
	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}
	*this << endl;
}

//  Output T3 element data
void COutputter::OutputT3Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl
		  << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;
	*this << "  SET       YOUNG'S     POISSON'S     THICKNESS" << endl
		  << " NUMBER     MODULUS        RATIO          T" << endl
		  << "               E              NU            H" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}

void COutputter::OutputQ9Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL CONSTANTS . . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POISSON'S" << endl
		  << " NUMBER     MODULUS          RATIO" << endl
		  << "               E              nu" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE   MATERIAL" << endl
		  << " NUMBER-N      1        2        3        4        5        6       7         8       9    SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
	{
		*this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
	}
	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::GetInstance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;

			case ElementTypes::Q4: // Q4 element
				*this << "  ELEMENT                 GAUSS POINT 1                                GAUSS POINT 2                                GAUSS POINT 3                                GAUSS POINT 4" << endl
					  << "  NUMBER      SIGMA_XX      SIGMA_YY      SIGMA_XY         SIGMA_XX      SIGMA_YY      SIGMA_XY         SIGMA_XX      SIGMA_YY      SIGMA_XY         SIGMA_XX      SIGMA_YY      SIGMA_XY" << endl;

				double stressQ4[12];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(stressQ4, Displacement);

					*this << setw(5) << Ele + 1;

					for (int i = 0; i < 4; i++)
					{
						*this << setw(17) << stressQ4[i * 3]
							  << setw(14) << stressQ4[i * 3 + 1]
							  << setw(14) << stressQ4[i * 3 + 2];
					}

					*this << endl;
				}

				*this << "                X             Y             Z                X             Y             Z                X             Y             Z                X             Y             Z" << endl;

				double gaussCoordsQ4[12];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.CalculateGaussPointCoordinates(gaussCoordsQ4);

					*this << setw(5) << Ele + 1;

					for (int i = 0; i < 4; i++)
					{
						*this << setw(17) << gaussCoordsQ4[i * 3]
							  << setw(14) << gaussCoordsQ4[i * 3 + 1]
							  << setw(14) << gaussCoordsQ4[i * 3 + 2];
					}

					*this << endl;
				}

				*this << "               D_X           D_Y           D_Z              D_X           D_Y           D_Z              D_X           D_Y           D_Z              D_X           D_Y           D_Z" << endl;
				
				double gaussDispQ4[12];
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.CalculateGaussPointDisplacement(gaussDispQ4, Displacement);

					*this << setw(5) << Ele + 1;

					for (int i = 0; i < 4; i++)
					{
						*this << setw(17) << gaussDispQ4[i * 3]
							  << setw(14) << gaussDispQ4[i * 3 + 1]
							  << setw(14) << gaussDispQ4[i * 3 + 2];
					}

					*this << endl;
				}


				*this << endl;
				break;

			case ElementTypes::T3: // T3 element
			/* TODO */
				double stress3T[3];
				#ifndef _TEST_
				*this << "  ELEMENT        LOCAL          ELEMENT        STRESS" << endl
					  << "  NUMBER         SXX            SYY            SXY" << endl;
				#else
				double GPPosition[9];
				double GPDisplacement[9];
				double weights3T[3];
				*this << "  ELEMENT    GP               GAUSS POINTS POSITION    "
					  << "                GAUSS POINTS DISPLACEMENTS       " 
					  << "               GAUSS POINTS STRESSES              INTEGRATE"
					  << std::endl
					  << "   INDEX   INDEX          X            Y              Z"
					  << "                DX           DY           DZ     "
					  << "          SXX           SYY           SXY          WEIGHTS"
					  << std::endl;
				#endif

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					#ifndef _TEST_
					dynamic_cast<CT3&>(Element).ElementStress(stress3T, Displacement);
					#else
					dynamic_cast<CTriangle&>(Element).ElementStress(stress3T, Displacement, GPPosition, GPDisplacement, weights3T);
					#endif
					CT3Material material = *dynamic_cast<CT3Material*>(Element.GetElementMaterial());
					
					#ifndef _TEST_
					*this << setw(5) << Ele + 1 << setw(20) << stress3T[0] 
					      << setw(15) << stress3T[1] << setw(15) << stress3T[2] << endl;
					#else
					for (unsigned GPIndex=0; GPIndex<3; GPIndex++)
					{
						*this << setw(6) << Ele+1 << setw(8) << GPIndex+1 
							  << setw(18) << GPPosition[3*GPIndex] 
							  << setw(14) << GPPosition[3*GPIndex + 1] 
							  << setw(14) << GPPosition[3*GPIndex + 2]
							  << setw(17) << GPDisplacement[3*GPIndex]
							  << setw(14) << GPDisplacement[3*GPIndex + 1]
							  << setw(14) << GPDisplacement[3*GPIndex + 2]
							  << setw(17) << stress3T[0]
							  << setw(14) << stress3T[1]
							  << setw(14) << stress3T[2]
							  << setw(14) << weights3T[GPIndex]
							  << std::endl;
					}
					#endif
				}

				*this << endl;
				break;


            case ElementTypes::Beam:
				*this << "  ELEMENT             FX                FY                FZ                MX                MY                MZ" << endl
				<< "  NUMBER         (Axial)         (Shear-Y)         (Shear-Z)         (Torsion)         (Moment-Y)         (Moment-Z)" << endl;

				// ����ÿ��Beam��Ԫ��6������������Fx, Fy, Fz, Mx, My, Mz
				double internalForces[6];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
				CElement& Element = EleGrp[Ele];
				Element.ElementStress(internalForces, Displacement); // ����ٶ�ElementStress���6������

				*this << setw(5) << Ele + 1;
				for (int i = 0; i < 6; ++i)
				{
				*this << setw(18) << internalForces[i];
				}
				*this << endl;
				}

				*this << endl;

				break;
			case ElementTypes::Q9:
				*this << "    ELEMENT   GAUSS P           GUASS POINTS POSITIONS"
					<< "                       GUASS POINTS STRESSES"
					<< endl;
				*this << "     NUMBER    INDEX " 
					<< "               SX'X'         SY'Y'        SX'Y'"
					<< endl;
				double stresses9Q[3*9];
				double Positions9Q[3*9];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					static_cast<CQ9&>(
						EleGrp.GetElement(Ele)).ElementStress(stresses9Q, Displacement);

					for (unsigned i=0; i<9; ++i) { // 9 gauss points
						*this << setw(8) << Ele + 1;
						*this << setw(10) << i+1;
						// *this << setw(17) << Positions9Q[i*3] << setw(14) << Positions9Q[i*3+1] << setw(14) << Positions9Q[i*3+2];
						*this << setw(17) << stresses9Q[i*3] << setw(14) << stresses9Q[i*3+1] << setw(14) << stresses9Q[i*3+2];
						*this << std::endl;
					}
				}

                *this << "          D_X           D_Y           D_Z" << endl;
				
				double gaussDispQ9[27];
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.CalculateGaussPointDisplacement(gaussDispQ9, Displacement);

					*this << setw(5) << Ele + 1 << endl;

					for (int i = 0; i < 9; i++)
					{
						*this << setw(17) << gaussDispQ9[i * 3]
							  << setw(14) << gaussDispQ9[i * 3 + 1]
							  << setw(14) << gaussDispQ9[i * 3 + 2] << endl;
					}

					*this << endl;
				}

				*this << endl;
				break;
			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif