

////ʱ��ԭ����δ��ɣ��д���һ������////


#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <TecplotOutputter.h>
#include <cmath>

using namespace std;

void NodalStress(

) {
    //�ӻ�ȡ�ڵ���Ϣ�͵�Ԫ��Ϣ
    vector<vector<double>> NodalData = getNodalData();

    vector<double> NodalStress;
    switch (ElementType)
    {
    case ElementTypes::Bar: // Bar element
        BarNodalStress();
        break;
    case ElementTypes::Q4: // 3T element
        Q4NodalStress(NodalData, ElementsData, NodalStress);
        break;
    case ElementTypes::T3: // 3T element
        vector<vector<double>> ElementsData = getT3ElementStressData();
        T3NodalStress(NodalData, ElementsData, NodalStress);
        break;
    case ElementTypes::H8: // 3T element
        H8NodalStress();
        break;
    case ElementTypes::Beam: // 3T element
        BeamNodalStress();
        break;
    case ElementTypes::Q9: // 3T element
        Q9NodalStress();
        break;
    default:
        *this << ElementType << " has not been implemented yet." << endl;
        break;
    }
}

//Q4��ԪΪ��Ӧ�䵥Ԫ��Ӧ���ڵ�Ԫ�䲻��������Ҫͨ����Χ��Ԫ��Ӧ��ƽ��ֵ���ع��ڵ�Ӧ��
void Q4NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {
    // ����������
    nodestress.clear();

    // ��ȡ�ڵ�����
    int nnp = nodedata.size();
    //��ȡ��Ԫ����
    int nel = elementdata.size();

    // Ϊÿ���ڵ��ʼ��Ӧ��ֵ
    for (int i = 0; i < nnp; i++) {
        // ��ӽڵ�š�X��Y����
        nodestress.push_back(nodedata[i][0]); // �ڵ��
        nodestress.push_back(nodedata[i][1] + 50 * nodedata[i][3]); // X����
        nodestress.push_back(nodedata[i][2] + 50 * nodedata[i][4]); // Y����

        // ��ʼ��Ӧ������Ϊ0
        nodestress.push_back(0.0); // ��_xx
        nodestress.push_back(0.0); // ��_yy
        nodestress.push_back(0.0); // ��_xy
        //��ӽڵ�X��Yλ��
        nodestress.push_back(nodedata[i][3]); // Xλ��
        nodestress.push_back(nodedata[i][4]); // Yλ��
        nodestress.push_back(0.0); //��ʼ������λ
    }
    //������Ԫ
    for (int j = 0; j < nel; j++) {
        // ��ȡ��ԪӦ������
        double sigma_xx = elementdata[j][5];
        double sigma_yy = elementdata[j][6];
        double sigma_xy = elementdata[j][7];

        // ��ÿ���ڵ����Ӧ���ۼ�
        for (int k = 1; k < 5; k++) {
            int nodeIndex = nodeList[k] - 1; // �ڵ��Ŵ�1��ʼ��ת��Ϊ����
            nodestress[nodeIndex * 9 + 3] += sigma_xx; // ��_xx
            nodestress[nodeIndex * 9 + 4] += sigma_yy; // ��_yy
            nodestress[nodeIndex * 9 + 5] += sigma_xy; // ��_xy
            nodestress[nodeIndex * 9 + 8] += 1.0; // ����λ��1
        }
    }
    for (int i = 0; i < nnp; i++) {
        // ����ƽ��Ӧ��
        if (nodestress[i * 9 + 8] > 0) {
            nodestress[i * 9 + 3] /= nodestress[i * 9 + 8]; // ��_xx
            nodestress[i * 9 + 4] /= nodestress[i * 9 + 8]; // ��_yy
            nodestress[i * 9 + 5] /= nodestress[i * 9 + 8]; // ��_xy
        }
    }
}


//T3��ԪӦ���ڵ�Ԫ�䲻��������Ҫͨ����Χ��Ԫ��Ӧ��ƽ��ֵ���ع��ڵ�Ӧ��
void T3NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {
    // ����������
    nodestress.clear();

    // ��ȡ�ڵ�����
    int nnp = nodedata.size();
    //��ȡ��Ԫ����
    int nel = elementdata.size();

    // Ϊÿ���ڵ��ʼ��Ӧ��ֵ
    for (int i = 0; i < nnp; i++) {
        // ��ӽڵ�š�X��Y����
        nodestress.push_back(nodedata[i][0]); // �ڵ��
        nodestress.push_back(nodedata[i][1] + 50 * nodedata[i][3]); // X����
        nodestress.push_back(nodedata[i][2] + 50 * nodedata[i][4]); // Y����

        // ��ʼ��Ӧ������Ϊ0
        nodestress.push_back(0.0); // ��_xx
        nodestress.push_back(0.0); // ��_yy
        nodestress.push_back(0.0); // ��_xy
        //��ӽڵ�X��Yλ��
        nodestress.push_back(nodedata[i][3]); // Xλ��
        nodestress.push_back(nodedata[i][4]); // Yλ��
        nodestress.push_back(0.0); //��ʼ������λ
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
            nodestress[nodeIndex * 9 + 3] += sigma_xx; // ��_xx
            nodestress[nodeIndex * 9 + 4] += sigma_yy; // ��_yy
            nodestress[nodeIndex * 9 + 5] += sigma_xy; // ��_xy
            nodestress[nodeIndex * 9 + 8] += 1.0; // ����λ��1
        }
    }
    for (int i = 0; i < nnp; i++) {
        // ����ƽ��Ӧ��
        if (nodestress[i * 9 + 8] > 0) {
            nodestress[i * 9 + 3] /= nodestress[i * 9 + 8]; // ��_xx
            nodestress[i * 9 + 4] /= nodestress[i * 9 + 8]; // ��_yy
            nodestress[i * 9 + 5] /= nodestress[i * 9 + 8]; // ��_xy
        }
    }
}


//H8��ԪӦ������������ֱ��ͨ���ڵ㴦��λ���ݶ������ڵ�Ӧ��
void H8NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {

}


//Beam��ԪӦ���ڵ�Ԫ�䲻��������Ҫͨ����Χ��Ԫ��Ӧ��ƽ��ֵ���ع��ڵ�Ӧ��
void BeamNodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {

}


//Q9��ԪӦ���ڵ�Ԫ������������ֱ��ͨ���ڵ㴦��λ���ݶ������ڵ�Ӧ��
void Q9NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {

}




vector<vector<double>> getNodalData() {
    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMNP = FEMData->GetNUMNP();

    vector<vector<double>> nodalData; //�ڵ��š�x���ꡢy���ꡢxλ�ơ�yλ��
    for (unsigned int np = 0; np < NUMNP; np++) {
        vector<double> nodeData;
        nodeData.push_back(np + 1); // �ڵ���
        nodeData.push_back(NodeList[np].GetX()); // x����
        nodeData.push_back(NodeList[np].GetY()); // y����
        // ��ӽڵ�λ�ƣ�����Displacement����������Ϊnp*3, np*3+1��Ӧx��yλ�ƣ�
        nodeData.push_back(Displacement[np * 3]);
        nodeData.push_back(Displacement[np * 3 + 1]);
        nodalData.push_back(nodeData);
    }
    return nodalData;
}

vector<vector<double>> getT3ElementStressData() {
    CDomain* FEMData = CDomain::GetInstance();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMEG = FEMData->GetNUMEG();
    vector<vector<double>> elementStressData; //��Ԫ��š��ڵ����б���ԪӦ��

    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes ElementType = EleGrp.GetElementType();

        for (unsigned int Ele = 0; Ele < NUME; Ele++) {
            CElement& Element = EleGrp[Ele];
            vector<double> elementData;
            elementData.push_back(Ele + 1); // ��Ԫ���

            // ��ȡ��Ԫ�ڵ���
            vector<int> nodeList = Element.GetNodeList();
            for (int node : nodeList) {
                elementData.push_back(node);
            }

            // ���㵥ԪӦ��
            double* stress;
            CT3::ElementStress(stress, Displacement);
            for (int i = 0; i < 3; i++)
            {
                elementData.push_back(stress[i]); // ���Ӧ������
            }

            elementStressData.push_back(elementData);
        }
    }
    return elementStressData;
}



vector<vector<double>> getQ4ElementStressData() {
    CDomain* FEMData = CDomain::GetInstance();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMEG = FEMData->GetNUMEG();
    vector<vector<double>> elementStressData; //��Ԫ��š��ڵ����б���ԪӦ��

    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes ElementType = EleGrp.GetElementType();

        for (unsigned int Ele = 0; Ele < NUME; Ele++) {
            CElement& Element = EleGrp[Ele];
            vector<double> elementData;
            elementData.push_back(Ele + 1); // ��Ԫ���

            // ��ȡ��Ԫ�ڵ���
            vector<int> nodeList = Element.GetNodeList();
            for (int node : nodeList) {
                elementData.push_back(node);
            }

            // ���㵥ԪӦ��
            double* stress;
            CT3::ElementStress(stress, Displacement);
            for (int i = 0; i < 3; i++)
            {
                elementData.push_back(stress[i]); // ���Ӧ������
            }

            elementStressData.push_back(elementData);
        }
    }
    return elementStressData;
}