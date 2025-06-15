

////时间原因尚未完成，有待进一步完善////


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
    //从获取节点信息和单元信息
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

//Q4单元为常应变单元，应变在单元间不连续，需要通过周围单元的应力平均值来重构节点应力
void Q4NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {
    // 清空输出数组
    nodestress.clear();

    // 获取节点数量
    int nnp = nodedata.size();
    //获取单元数量
    int nel = elementdata.size();

    // 为每个节点初始化应力值
    for (int i = 0; i < nnp; i++) {
        // 添加节点号、X、Y坐标
        nodestress.push_back(nodedata[i][0]); // 节点号
        nodestress.push_back(nodedata[i][1] + 50 * nodedata[i][3]); // X坐标
        nodestress.push_back(nodedata[i][2] + 50 * nodedata[i][4]); // Y坐标

        // 初始化应力分量为0
        nodestress.push_back(0.0); // σ_xx
        nodestress.push_back(0.0); // σ_yy
        nodestress.push_back(0.0); // σ_xy
        //添加节点X、Y位移
        nodestress.push_back(nodedata[i][3]); // X位移
        nodestress.push_back(nodedata[i][4]); // Y位移
        nodestress.push_back(0.0); //初始化计数位
    }
    //遍历单元
    for (int j = 0; j < nel; j++) {
        // 获取单元应力分量
        double sigma_xx = elementdata[j][5];
        double sigma_yy = elementdata[j][6];
        double sigma_xy = elementdata[j][7];

        // 对每个节点进行应力累加
        for (int k = 1; k < 5; k++) {
            int nodeIndex = nodeList[k] - 1; // 节点编号从1开始，转换为索引
            nodestress[nodeIndex * 9 + 3] += sigma_xx; // σ_xx
            nodestress[nodeIndex * 9 + 4] += sigma_yy; // σ_yy
            nodestress[nodeIndex * 9 + 5] += sigma_xy; // σ_xy
            nodestress[nodeIndex * 9 + 8] += 1.0; // 计数位加1
        }
    }
    for (int i = 0; i < nnp; i++) {
        // 计算平均应力
        if (nodestress[i * 9 + 8] > 0) {
            nodestress[i * 9 + 3] /= nodestress[i * 9 + 8]; // σ_xx
            nodestress[i * 9 + 4] /= nodestress[i * 9 + 8]; // σ_yy
            nodestress[i * 9 + 5] /= nodestress[i * 9 + 8]; // σ_xy
        }
    }
}


//T3单元应变在单元间不连续，需要通过周围单元的应力平均值来重构节点应力
void T3NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {
    // 清空输出数组
    nodestress.clear();

    // 获取节点数量
    int nnp = nodedata.size();
    //获取单元数量
    int nel = elementdata.size();

    // 为每个节点初始化应力值
    for (int i = 0; i < nnp; i++) {
        // 添加节点号、X、Y坐标
        nodestress.push_back(nodedata[i][0]); // 节点号
        nodestress.push_back(nodedata[i][1] + 50 * nodedata[i][3]); // X坐标
        nodestress.push_back(nodedata[i][2] + 50 * nodedata[i][4]); // Y坐标

        // 初始化应力分量为0
        nodestress.push_back(0.0); // σ_xx
        nodestress.push_back(0.0); // σ_yy
        nodestress.push_back(0.0); // σ_xy
        //添加节点X、Y位移
        nodestress.push_back(nodedata[i][3]); // X位移
        nodestress.push_back(nodedata[i][4]); // Y位移
        nodestress.push_back(0.0); //初始化计数位
    }
    //遍历单元
    for (int j = 0; j < nel; j++) {
        // 获取单元应力分量
        double sigma_xx = elementdata[j][4];
        double sigma_yy = elementdata[j][5];
        double sigma_xy = elementdata[j][6];

        // 对每个节点进行应力累加
        for (int k = 1; k < 4; k++) {
            int nodeIndex = nodeList[k] - 1; // 节点编号从1开始，转换为索引
            nodestress[nodeIndex * 9 + 3] += sigma_xx; // σ_xx
            nodestress[nodeIndex * 9 + 4] += sigma_yy; // σ_yy
            nodestress[nodeIndex * 9 + 5] += sigma_xy; // σ_xy
            nodestress[nodeIndex * 9 + 8] += 1.0; // 计数位加1
        }
    }
    for (int i = 0; i < nnp; i++) {
        // 计算平均应力
        if (nodestress[i * 9 + 8] > 0) {
            nodestress[i * 9 + 3] /= nodestress[i * 9 + 8]; // σ_xx
            nodestress[i * 9 + 4] /= nodestress[i * 9 + 8]; // σ_yy
            nodestress[i * 9 + 5] /= nodestress[i * 9 + 8]; // σ_xy
        }
    }
}


//H8单元应变连续，可以直接通过节点处的位移梯度来求解节点应力
void H8NodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {

}


//Beam单元应变在单元间不连续，需要通过周围单元的应力平均值来重构节点应力
void BeamNodalStress(
    const vector<vector<double>> nodedata,
    const vector<vector<double>> elementdata,
    vector<double>& nodestress
) {

}


//Q9单元应变在单元间连续，可以直接通过节点处的位移梯度来求解节点应力
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

    vector<vector<double>> nodalData; //节点编号、x坐标、y坐标、x位移、y位移
    for (unsigned int np = 0; np < NUMNP; np++) {
        vector<double> nodeData;
        nodeData.push_back(np + 1); // 节点编号
        nodeData.push_back(NodeList[np].GetX()); // x坐标
        nodeData.push_back(NodeList[np].GetY()); // y坐标
        // 添加节点位移（假设Displacement数组中索引为np*3, np*3+1对应x和y位移）
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
    vector<vector<double>> elementStressData; //单元编号、节点编号列表、单元应力

    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes ElementType = EleGrp.GetElementType();

        for (unsigned int Ele = 0; Ele < NUME; Ele++) {
            CElement& Element = EleGrp[Ele];
            vector<double> elementData;
            elementData.push_back(Ele + 1); // 单元编号

            // 获取单元节点编号
            vector<int> nodeList = Element.GetNodeList();
            for (int node : nodeList) {
                elementData.push_back(node);
            }

            // 计算单元应力
            double* stress;
            CT3::ElementStress(stress, Displacement);
            for (int i = 0; i < 3; i++)
            {
                elementData.push_back(stress[i]); // 添加应力分量
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
    vector<vector<double>> elementStressData; //单元编号、节点编号列表、单元应力

    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes ElementType = EleGrp.GetElementType();

        for (unsigned int Ele = 0; Ele < NUME; Ele++) {
            CElement& Element = EleGrp[Ele];
            vector<double> elementData;
            elementData.push_back(Ele + 1); // 单元编号

            // 获取单元节点编号
            vector<int> nodeList = Element.GetNodeList();
            for (int node : nodeList) {
                elementData.push_back(node);
            }

            // 计算单元应力
            double* stress;
            CT3::ElementStress(stress, Displacement);
            for (int i = 0; i < 3; i++)
            {
                elementData.push_back(stress[i]); // 添加应力分量
            }

            elementStressData.push_back(elementData);
        }
    }
    return elementStressData;
}