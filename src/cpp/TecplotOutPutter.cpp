
#include "TecplotOutPutter.h"   

using namespace std;



void TecplotOutPutter() {

    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    //build nodal data matrix with serial number, coordinates, displacements, stress counter for average, and stress components
    vector<vector<double>> NodeStressData = getNodalData();

    // get number of element groups and element group list
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // traverse all element groups
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();

        // build element data matrix with element number, node numbers, and stress components
        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            T3NodalStress(ElementsData, NodeStressData); // 处理T3单元节点应力，整合进总体的节点信息数组
            break;

   //     case ElementTypes::Q4:  // Q4
   //         ElementsData = getQ4ElementStressData(eleGrp, displacement);
   //         Q4NodalStress(ElementsData, NodeStressData); // 处理Q4单元节点应力，整合进总体的节点信息数组
			//break;

        case ElementTypes::Bar:  // Bar(to be supplemented)
            // ElementsData = processBarElementGroup(eleGrp, displacement);
            break;

        default:  // other unhandled types
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    for (unsigned int i = 0; i < NodeStressData.size(); i++) {
        // calculate average stress with the stress counter
        if (NodeStressData[i][7] > 1.0) {
            NodeStressData[i][8] /= NodeStressData[i][7]; // σ_xx
            NodeStressData[i][9] /= NodeStressData[i][7]; // σ_yy
            NodeStressData[i][10] /= NodeStressData[i][7]; // σ_xy
        }
    }
    OutputTecplotDatafile(NodeStressData);
}



void OutputTecplotDatafile(vector<vector<double> >& nodalData) {

	ofstream outFile("output.dat");
	if (!outFile) {
		cerr << "Error opening output file." << endl;
		return;
	}

	// write title and variable names
	outFile << "TITLE = \"STAPpp\"\n";
	outFile << "VARIABLES = \"Node\", \"X\", \"Y\", \"Z\", \"Displacement_X\", \"Displacement_Y\", \"Displacement_Z\", \"Sigma_X\", \"Sigma_Y\", \"Tau_XY\"\n";
	
	

    // write nodal data
    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    

    // get number of element groups and element group list
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // traverse all element groups
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();


        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3
            
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            // write T3 element basic information
            outFile << "ZONE T=\"" << grpIdx << "\", ZONETYPE=FETRIANGLE, NODES=" << nodalData.size() << ", ELEMENTS=" << ElementsData.size() << ", DATAPACKING=POINT\n";
            // write node data
            for (const auto& node : nodalData) {
                for (unsigned int i = 1; i < node.size(); ++i) {
                    if (i != 7) {
                        outFile << node[i];
                        if (i < node.size() - 1) outFile << " ";
                    } // skip counter
                }
                outFile << "\n";
            }
            // write IEN
            for (const auto& element : ElementsData) {
				for (int i = 1; i < 4; ++i) {
					outFile << element[i];
					if (i < element.size() - 1) outFile << " ";
				}
				outFile << "\n";
			}
            break;

        //case ElementTypes::Q4:  // Q4
        //    
        //    ElementsData = getQ4ElementStressData(eleGrp, displacement);
        //    // write Q4 element basic information
        //    outFile << "ZONE T=\"" << grpIdx + 1 << "\", ZONETYPE=FEQUADRILATERAL, NODES=" << nodalData.size() << ", ELEMENTS=" << ElementsData.size() << ", DATAPACKING=POINT\n";
        //    // write node data
        //    for (const auto& node : nodalData) {
        //        for (unsigned int i = 1; i < node.size(); ++i) {
        //            if (i != 7) {
        //                outFile << node[i];
        //                if (i < node.size() - 1) outFile << " ";
        //            } // skip counter
        //        }
        //        outFile << "\n";
        //    }
        //    // write IEN
        //    for (const auto& element : ElementsData) {
        //        for (int i = 1; i < 5; ++i) {
        //            outFile << element[i];
        //            if (i < element.size() - 1) outFile << " ";
        //        }
        //        outFile << "\n";
        //    }
        //    break;

        case ElementTypes::Bar:  // Bar(to be supplemented)
            
            break;

        default:  // other unhandled types
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    
	outFile.close();
}



vector<vector<double>> getNodalData() {

    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMNP = FEMData->GetNUMNP();

    vector<vector<double>> nodalData; //node data matrix with serial number, coordinate x, coordinate y, coordinate z, displacement x, displacement y, displacement z, stress counter for average, sigma_x, sigma_y, tau_xy
    for (unsigned int np = 0; np < NUMNP; np++) {
        vector<double> nodeData;
        nodeData.push_back(np + 1); // serial number (from 1)
        nodeData.push_back(NodeList[np].XYZ[0]); // x
        nodeData.push_back(NodeList[np].XYZ[1]); // y
        nodeData.push_back(NodeList[np].XYZ[2]); // z
        // node displacement
        nodeData.push_back(Displacement[np * 6]); // u_x
        nodeData.push_back(Displacement[np * 6 + 1]); // u_y
        nodeData.push_back(Displacement[np * 6 + 2]); // u_z
        nodeData.push_back(0.0);// initialize stress counter for average
        nodeData.push_back(0.0); // initialize sigma_x
        nodeData.push_back(0.0); // initialize sigma_y
        nodeData.push_back(0.0); // initialize tau_xy
        nodalData.push_back(nodeData);
    }
    return nodalData;
}




vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement) {

    vector<vector<double>> t3ElementData;
    unsigned int numElements = eleGrp.GetNUME();

    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
        // get current T3 element (insure it is of type T3)
        CT3* t3Element = dynamic_cast<CT3*>(&eleGrp[eleIdx]);
        if (!t3Element) continue;

        vector<double> elementInfo;

        // add element number (group index)
        elementInfo.push_back(eleIdx + 1);

        // add node numbers
        CNode** nodes_= t3Element->GetNodes();
        elementInfo.push_back(nodes_[0]->NodeNumber);  // node1 of the T3 element
        elementInfo.push_back(nodes_[1]->NodeNumber);  // node2 of the T3 element
        elementInfo.push_back(nodes_[2]->NodeNumber);  // node3 of the T3 element

        // calculate and add stress components
        double stress[3] = { 0 };  // store stress components
        t3Element->ElementStress(stress, displacement);

        elementInfo.push_back(stress[0]);  // sigma_x
        elementInfo.push_back(stress[1]);  // sigma_y
        elementInfo.push_back(stress[2]);  // tau_xy

        t3ElementData.push_back(elementInfo);
    }

    return t3ElementData;
}



//vector<vector<double>> getQ4ElementStressData(CElementGroup& eleGrp, double* displacement) {
//
//    vector<vector<double>> q4ElementData;
//    unsigned int numElements = eleGrp.GetNUME();
//
//    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
//        // get current Q4 element (ensure it is of type Q4)
//        CQ4* q4Element = dynamic_cast<CQ4*>(&eleGrp[eleIdx]);
//        if (!q4Element) continue;
//
//        vector<double> elementInfo;
//
//        // add element number (group index)
//        elementInfo.push_back(eleIdx + 1);
//
//        // add node numbers
//        CNode** nodes_ = q4Element->GetNodes();
//        elementInfo.push_back(nodes_[0]->NodeNumber);  // node1 of the Q4 element
//        elementInfo.push_back(nodes_[1]->NodeNumber);  // node2 of the Q4 element
//        elementInfo.push_back(nodes_[2]->NodeNumber);  // node3 of the Q4 element
//        elementInfo.push_back(nodes_[3]->NodeNumber);  // node4 of the Q4 element
//
//        // add stress components (4 Gauss points * 3 components)
//        double stress[12] = { 0 };  // store stress components
//        q4Element->ElementStress(stress, displacement);
//
//        // push stress components into elementInfo
//        for (int i = 0; i < 12; i++) {
//            elementInfo.push_back(stress[i]);
//        }
//
//        q4ElementData.push_back(elementInfo);
//    }
//
//    return q4ElementData;
//}







// T3 element strains are discontinuous between elements. 
// Node stresses are reconstructed by averaging stresses from surrounding elements.
void T3NodalStress(
    const vector<vector<double>> elementdata,
    vector<vector<double>>& nodestress
) {
    unsigned int nnp = nodestress.size();
    unsigned int nel = elementdata.size();

    // apply displacement amplification to each node
    for (unsigned int i = 0; i < nnp; i++) {
        // add node number, X, Y coordinates
        nodestress[i][1] += 50 * nodestress[i][4];
        nodestress[i][2] += 50 * nodestress[i][5];
    }
    // traverse each element
    for (unsigned int j = 0; j < nel; j++) {
        // get element stress components
        double sigma_xx = elementdata[j][4];
        double sigma_yy = elementdata[j][5];
        double sigma_xy = elementdata[j][6];

        // add up stresses to each node
        for (int k = 1; k < 4; k++) {
            int nodeIndex = elementdata[j][k] - 1; // node numbers start from 1, convert to index
            nodestress[nodeIndex][8] += sigma_xx; // sigma_xx
            nodestress[nodeIndex][9] += sigma_yy; // sigma_yy
            nodestress[nodeIndex][10] += sigma_xy; // tau_xy
            nodestress[nodeIndex][7] += 1.0; // count + 1
        }
    }
}



// Q4 elements are constant-strain elements with discontinuous strains between elements.
// Node stresses are reconstructed by averaging stresses from the nearest Gauss points of surrounding elements.
//void Q4NodalStress(
//    const vector<vector<double>> elementdata,
//    vector<vector<double>>& nodestress
//) {
//
//    int nnp = nodestress.size();
//
//    int nel = elementdata.size();
//
//    // apply displacement amplification to each node
//    for (int i = 0; i < nnp; i++) {
//        // add node number, X, Y coordinates
//        nodestress[i][1] += 50 * nodestress[i][4]; // X坐标
//        nodestress[i][2] += 50 * nodestress[i][5]; // Y坐标
//    }
//    // traverse each element
//    for (int j = 0; j < nel; j++) {
//        // add up stresses to each node
//        for (int k = 1; k < 5; k++) {
//            int nodeIndex = elementdata[j][k] - 1; // node numbers start from 1, convert to index
//            nodestress[nodeIndex][8] += elementdata[j][4 * k + 1]; // sigma_xx
//            nodestress[nodeIndex][9] += elementdata[j][4 * k + 2]; // sigma_yy
//            nodestress[nodeIndex][10] += elementdata[j][4 * k + 3]; // tau_xy
//            nodestress[nodeIndex][7] += 1.0; // count + 1
//        }
//    }
//}


