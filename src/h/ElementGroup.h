/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <fstream>

#include "Element.h"
#include "Elements/Bar.h"
#include "Elements/Q4.h"
#include "Elements/T3.h"
#include "Elements/Beam.h"
#include "Elements/Q9.h"
#include "Material.h"
#include "Node.h"
using namespace std;

//! Define set of element types
enum ElementTypes
{
    UNDEFINED = 0,
    Bar = 1,    // Bar element
    Q4 = 2,     // 4Q element
    T3 = 3,     // 3T element
    H8 = 4,     // 8H element
    Beam = 5,   // Beam element
    Plate = 6,  // Plate element
    Shell = 7,  // Shell elment
    Q9 = 8      // Q9 element
};

//! Element group class
class CElementGroup
{
private:
    //! List of all nodes in the domain, obtained from CDomain object
    static CNode* NodeList_;

    //! Element type of this group
    ElementTypes ElementType_;

    //! Size of an Element object in this group
    std::size_t ElementSize_;

    //! Number of elements in this group
    unsigned int NUME_;

    //! Element List in this group
    CElement* ElementList_;

    //! Number of material/section property sets in this group
    unsigned int NUMMAT_;

    //! Material list in this group
    CMaterial* MaterialList_;

    //! Size of an Material object in this group
    std::size_t MaterialSize_;

public:
    //! Constructor
    CElementGroup();

    //! Deconstructor
    ~CElementGroup();

    //! Read element group data from stream Input
    bool Read(ifstream& Input);

    //! Calculate the size of the derived element class and material class
    void CalculateMemberSize();

    //! Allocate array of derived elements
    void AllocateElements(std::size_t size);

    //! Allocate array of derived materials
    void AllocateMaterials(std::size_t size);

    //! Read element data from the input data file
    bool ReadElementData(ifstream& Input);

    //! Return element type of this group
    ElementTypes GetElementType() { return ElementType_; }

    //! Return the number of elements in the group
    unsigned int GetNUME() { return NUME_; }

    //! operator []
    //! For the sake of efficiency, the index bounds are not checked
    CElement& operator[](unsigned int i);

    CElement& GetElement(unsigned int i);

    //! Return the index-th material in this group
    CMaterial& GetMaterial(unsigned int i);

    //! Return the number of material/section property setss in this element group
    unsigned int GetNUMMAT() { return NUMMAT_; }
};
