#include <vector>
#include <iostream>
#include "Domain.h"
#include "ElementGroup.h"
#include "Elements/T3.h"  // ���� T3 ��Ԫ����

using namespace std;

// ǰ������
vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement);

// ��ȡ���нڵ���Ϣ
vector<vector<double>> TecplotOutPutter() {
    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    vector<vector<double>> NodeStressData = getNodalData();

    // ��ȡ��Ԫ���б�
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // �������е�Ԫ��
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();

        // ���ݵ�Ԫ����ѡ����ʽ
        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3��Ԫ����
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            T3NodalStress(ElementsData, NodeStressData); // ����T3��Ԫ�ڵ�Ӧ�������Ͻ�����Ľڵ���Ϣ����
            break;

        case ElementTypes::Q4:  // Q4��Ԫ����
            ElementsData = getQ4ElementStressData(eleGrp, displacement);
            Q4NodalStress(ElementsData, NodeStressData); // ����Q4��Ԫ�ڵ�Ӧ�������Ͻ�����Ľڵ���Ϣ����
			break;

        case ElementTypes::Bar:  // �˵�Ԫʾ�����ɸ�����Ҫ��չ��
            // ElementsData = processBarElementGroup(eleGrp, displacement);
            break;

        default:  // ����δ��������
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    for (int i = 0; i < nnp; i++) {
        // ����ƽ��Ӧ��,�ü�������ƽ��
        if (NodeStressData[i][7] > 0) {
            NodeStressData[i][8] /= NodeStressData[i * 9 + 8]; // ��_xx
            NodeStressData[i][9] /= NodeStressData[i * 9 + 8]; // ��_yy
            NodeStressData[i][10] /= NodeStressData[i * 9 + 8]; // ��_xy
        }
    }
    return NodeStressData;
}



void OutputTecplotDatafile(vector<vector<double> >& nodalData) {
	// ���Tecplot�����ļ�
	ofstream outFile("output.dat");
	if (!outFile) {
		cerr << "Error opening output file." << endl;
		return;
	}

	// д�������
	outFile << "TITLE = \"Nodal Data\"\n";
	outFile << "VARIABLES = \"Node\", \"X\", \"Y\", \"Z\", \"Displacement_X\", \"Displacement_Y\", \"Displacement_Z\", \"Sigma_X\", \"Sigma_Y\", \"Tau_XY\"\n";
	
	// д��ڵ�����
	for (const auto& node : nodalData) {
		for (size_t i = 0; i < node.size(); ++i) {
			outFile << node[i];
			if (i < node.size() - 1) outFile << " "; // �ָ���
		}
		outFile << "\n"; // ����
	}

    //д�뵥Ԫ��Ϣ
    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    vector<vector<double>> NodeStressData = getNodalData();

    // ��ȡ��Ԫ���б�
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // �������е�Ԫ��
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();

        // ���ݵ�Ԫ����ѡ����ʽ
        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3��Ԫ����
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            T3NodalStress(ElementsData, NodeStressData); // ����T3��Ԫ�ڵ�Ӧ�������Ͻ�����Ľڵ���Ϣ����
            break;

        case ElementTypes::Q4:  // Q4��Ԫ����
            ElementsData = getQ4ElementStressData(eleGrp, displacement);
            Q4NodalStress(ElementsData, NodeStressData); // ����Q4��Ԫ�ڵ�Ӧ�������Ͻ�����Ľڵ���Ϣ����
            break;

        case ElementTypes::Bar:  // �˵�Ԫʾ�����ɸ�����Ҫ��չ��
            // ElementsData = processBarElementGroup(eleGrp, displacement);
            break;

        default:  // ����δ��������
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    outFile << "ZONE T=\"Element Data\"\n";

	outFile.close();
}


// ��ȡ�ڵ�����
vector<vector<double>> getNodalData() {
    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMNP = FEMData->GetNUMNP();

    vector<vector<double>> nodalData; //�ڵ��š�x���ꡢy���ꡢz���ꡢxλ�ơ�yλ�ơ�zλ�ơ���_x����_y����_xy����
    for (unsigned int np = 0; np < NUMNP; np++) {
        vector<double> nodeData;
        nodeData.push_back(np + 1); // �ڵ���
        nodeData.push_back(NodeList[np].XYZ[0]); // x����
        nodeData.push_back(NodeList[np].XYZ[1]); // y����
        nodeData.push_back(NodeList[np].XYZ[2]); // z����
        // ��ӽڵ�λ�ƣ�����Displacement����������Ϊnp*6, np*6+1��np*6+2��Ӧx��y��zλ�ƣ�
        nodeData.push_back(Displacement[np * 6]);
        nodeData.push_back(Displacement[np * 6 + 1]);
        nodeData.push_back(Displacement[np * 6 + 2]);
        nodeData.push_back(0.0);//��ʼ��������
        nodeData.push_back(0.0); // ��ʼ����_x
        nodeData.push_back(0.0); // ��ʼ����_y
        nodeData.push_back(0.0); // ��ʼ����_xy
        nodalData.push_back(nodeData);
    }
    return nodalData;
}




// ����T3��Ԫ��ĺ���
vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement) {
    vector<vector<double>> t3ElementData;
    unsigned int numElements = eleGrp.GetNUME();

    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
        // ��ȡ��ǰ��Ԫ��ȷ����T3���ͣ�
        CT3* t3Element = dynamic_cast<CT3*>(&eleGrp[eleIdx]);
        if (!t3Element) continue;

        vector<double> elementInfo;

        // 1. ��ӵ�Ԫ��ţ����ڱ�ţ�
        elementInfo.push_back(eleIdx + 1);

        // 2. ��ӽڵ���
        elementInfo.push_back(t3Element->nodes_[0]->NodeNumber);  // �ڵ�1
        elementInfo.push_back(t3Element->nodes_[1]->NodeNumber);  // �ڵ�2
        elementInfo.push_back(t3Element->nodes_[2]->NodeNumber);  // �ڵ�3

        // 3. ���㲢���Ӧ������
        double stress[3] = { 0 };  // �洢Ӧ������
        t3Element->ElementStress(stress, displacement);

        elementInfo.push_back(stress[0]);  // ��_x
        elementInfo.push_back(stress[1]);  // ��_y
        elementInfo.push_back(stress[2]);  // ��_xy

        t3ElementData.push_back(elementInfo);
    }

    return t3ElementData;
}


// ����Q4��Ԫ��ĺ���
vector<vector<double>> getQ4ElementStressData(CElementGroup& eleGrp, double* displacement) {
    vector<vector<double>> q4ElementData;
    unsigned int numElements = eleGrp.GetNUME();

    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
        // ��ȡ��ǰQ4��Ԫ
        CQ4* q4Element = dynamic_cast<CQ4*>(&eleGrp[eleIdx]);
        if (!q4Element) continue;  // ȷ��ת���ɹ�

        vector<double> elementInfo;

        // 1. ��ӵ�Ԫ��ţ����ڱ�ţ�
        elementInfo.push_back(eleIdx + 1);

        // 2. ��ӽڵ���
        elementInfo.push_back(q4Element->nodes_[0]->NodeNumber);  // �ڵ�1
        elementInfo.push_back(q4Element->nodes_[1]->NodeNumber);  // �ڵ�2
        elementInfo.push_back(q4Element->nodes_[2]->NodeNumber);  // �ڵ�3
        elementInfo.push_back(q4Element->nodes_[3]->NodeNumber);  // �ڵ�4

        // 3. ���㲢���Ӧ��������ÿ����Ԫ��4����˹�㣬ÿ����3��Ӧ��������
        double stress[12] = { 0 };  // �洢12��Ӧ��������4����˹�� �� 3������
        q4Element->ElementStress(stress, displacement);

        // �������Ӧ������
        for (int i = 0; i < 12; i++) {
            elementInfo.push_back(stress[i]);
        }

        q4ElementData.push_back(elementInfo);
    }

    return q4ElementData;
}







//T3��ԪӦ���ڵ�Ԫ�䲻��������Ҫͨ����Χ��Ԫ��Ӧ��ƽ��ֵ���ع��ڵ�Ӧ��
void T3NodalStress(
    const vector<vector<double>> elementdata,
    vector<vector<double>>& nodestress
) {
    // ��ȡ�ڵ�����
    int nnp = nodestress.size();
    //��ȡ��Ԫ����
    int nel = elementdata.size();

    // Ϊÿ���ڵ����λ�ƷŴ�
    for (int i = 0; i < nnp; i++) {
        // ��ӽڵ�š�X��Y����
        nodestress[i][1] += 50 * nodestress[i][4]); // X����
        nodestress[i][2] += 50 * nodestress[i][5]); // Y����
    }
    //������Ԫ
    for (int j = 0; j < nel; j++) {
        // ��ȡ��ԪӦ������
        double sigma_xx = elementdata[j][4];
        double sigma_yy = elementdata[j][5];
        double sigma_xy = elementdata[j][6];

        // ��ÿ���ڵ����Ӧ���ۼ�
        for (int k = 1; k < 4; k++) {
            int nodeIndex = nodeList[k] - 1; // �ڵ��Ŵ�1��ʼ��ת��Ϊ����
            nodestress[nodeIndex][8] += sigma_xx; // ��_xx
            nodestress[nodeIndex][9] += sigma_yy; // ��_yy
            nodestress[nodeIndex][10] += sigma_xy; // ��_xy
            nodestress[nodeIndex][7] += 1.0; // ����λ��1
        }
    }
}



//Q4��ԪΪ��Ӧ�䵥Ԫ��Ӧ���ڵ�Ԫ�䲻��������Ҫͨ����Χ��Ԫ������ĸ�˹���Ӧ��ƽ��ֵ���ع��ڵ�Ӧ��
void Q4NodalStress(
    const vector<vector<double>> elementdata,
    vector<vector<double>>& nodestress
) {
    // ��ȡ�ڵ�����
    int nnp = nodestress.size();
    //��ȡ��Ԫ����
    int nel = elementdata.size();

    // Ϊÿ���ڵ����λ�ƷŴ�
    for (int i = 0; i < nnp; i++) {
        // ��ӽڵ�š�X��Y����
        nodestress[i][1] += 50 * nodestress[i][4]); // X����
        nodestress[i][2] += 50 * nodestress[i][5]); // Y����
    }
    //������Ԫ
    for (int j = 0; j < nel; j++) {
        // ��ÿ���ڵ����Ӧ���ۼ�
        for (int k = 1; k < 5; k++) {
            int nodeIndex = nodeList[k] - 1; // �ڵ��Ŵ�1��ʼ��ת��Ϊ����
            nodestress[nodeIndex][8] += elementdata[j][4 * k + 1]; // ��_xx
            nodestress[nodeIndex][9] += elementdata[j][4 * k + 2]; // ��_yy
            nodestress[nodeIndex][10] += elementdata[j][4 * k + 3]; // ��_xy
            nodestress[nodeIndex][7] += 1.0; // ����λ��1
        }
    }
}


