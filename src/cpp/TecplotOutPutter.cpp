#include <vector>
#include <iostream>
#include "Domain.h"
#include "ElementGroup.h"
#include "Elements/T3.h"  // 包含 T3 单元定义

using namespace std;

// 前置声明
vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement);

// 获取所有节点信息
vector<vector<double>> TecplotOutPutter() {
    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    vector<vector<double>> NodeStressData = getNodalData();

    // 获取单元组列表
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // 遍历所有单元组
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();

        // 根据单元类型选择处理方式
        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3单元处理
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            T3NodalStress(ElementsData, NodeStressData); // 处理T3单元节点应力，整合进总体的节点信息数组
            break;

        case ElementTypes::Q4:  // Q4单元处理
            ElementsData = getQ4ElementStressData(eleGrp, displacement);
            Q4NodalStress(ElementsData, NodeStressData); // 处理Q4单元节点应力，整合进总体的节点信息数组
			break;

        case ElementTypes::Bar:  // 杆单元示例（可根据需要扩展）
            // ElementsData = processBarElementGroup(eleGrp, displacement);
            break;

        default:  // 其他未处理类型
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    for (int i = 0; i < nnp; i++) {
        // 计算平均应力,用计数器做平均
        if (NodeStressData[i][7] > 0) {
            NodeStressData[i][8] /= NodeStressData[i * 9 + 8]; // σ_xx
            NodeStressData[i][9] /= NodeStressData[i * 9 + 8]; // σ_yy
            NodeStressData[i][10] /= NodeStressData[i * 9 + 8]; // σ_xy
        }
    }
    return NodeStressData;
}



void OutputTecplotDatafile(vector<vector<double> >& nodalData) {
	// 输出Tecplot数据文件
	ofstream outFile("output.dat");
	if (!outFile) {
		cerr << "Error opening output file." << endl;
		return;
	}

	// 写入标题行
	outFile << "TITLE = \"Nodal Data\"\n";
	outFile << "VARIABLES = \"Node\", \"X\", \"Y\", \"Z\", \"Displacement_X\", \"Displacement_Y\", \"Displacement_Z\", \"Sigma_X\", \"Sigma_Y\", \"Tau_XY\"\n";
	
	// 写入节点数据
	for (const auto& node : nodalData) {
		for (size_t i = 0; i < node.size(); ++i) {
			outFile << node[i];
			if (i < node.size() - 1) outFile << " "; // 分隔符
		}
		outFile << "\n"; // 换行
	}

    //写入单元信息
    CDomain* FEMData = CDomain::GetInstance();
    double* displacement = FEMData->GetDisplacement();
    vector<vector<double>> NodeStressData = getNodalData();

    // 获取单元组列表
    unsigned int numGroups = FEMData->GetNUMEG();
    CElementGroup* eleGrpList = FEMData->GetEleGrpList();

    // 遍历所有单元组
    for (unsigned int grpIdx = 0; grpIdx < numGroups; grpIdx++) {
        CElementGroup& eleGrp = eleGrpList[grpIdx];
        ElementTypes eleType = eleGrp.GetElementType();

        // 根据单元类型选择处理方式
        vector<vector<double>> ElementsData;

        switch (eleType) {
        case ElementTypes::T3:  // T3单元处理
            ElementsData = getT3ElementStressData(eleGrp, displacement);
            T3NodalStress(ElementsData, NodeStressData); // 处理T3单元节点应力，整合进总体的节点信息数组
            break;

        case ElementTypes::Q4:  // Q4单元处理
            ElementsData = getQ4ElementStressData(eleGrp, displacement);
            Q4NodalStress(ElementsData, NodeStressData); // 处理Q4单元节点应力，整合进总体的节点信息数组
            break;

        case ElementTypes::Bar:  // 杆单元示例（可根据需要扩展）
            // ElementsData = processBarElementGroup(eleGrp, displacement);
            break;

        default:  // 其他未处理类型
            cerr << "Unhandled element type: " << static_cast<int>(eleType)
                << " in group " << grpIdx + 1 << endl;
            break;
        }

    }
    outFile << "ZONE T=\"Element Data\"\n";

	outFile.close();
}


// 获取节点数据
vector<vector<double>> getNodalData() {
    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMNP = FEMData->GetNUMNP();

    vector<vector<double>> nodalData; //节点编号、x坐标、y坐标、z坐标、x位移、y位移、z位移、σ_x、σ_y、τ_xy……
    for (unsigned int np = 0; np < NUMNP; np++) {
        vector<double> nodeData;
        nodeData.push_back(np + 1); // 节点编号
        nodeData.push_back(NodeList[np].XYZ[0]); // x坐标
        nodeData.push_back(NodeList[np].XYZ[1]); // y坐标
        nodeData.push_back(NodeList[np].XYZ[2]); // z坐标
        // 添加节点位移（假设Displacement数组中索引为np*6, np*6+1，np*6+2对应x和y和z位移）
        nodeData.push_back(Displacement[np * 6]);
        nodeData.push_back(Displacement[np * 6 + 1]);
        nodeData.push_back(Displacement[np * 6 + 2]);
        nodeData.push_back(0.0);//初始化计数器
        nodeData.push_back(0.0); // 初始化σ_x
        nodeData.push_back(0.0); // 初始化σ_y
        nodeData.push_back(0.0); // 初始化τ_xy
        nodalData.push_back(nodeData);
    }
    return nodalData;
}




// 处理T3单元组的函数
vector<vector<double>> getT3ElementStressData(CElementGroup& eleGrp, double* displacement) {
    vector<vector<double>> t3ElementData;
    unsigned int numElements = eleGrp.GetNUME();

    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
        // 获取当前单元（确保是T3类型）
        CT3* t3Element = dynamic_cast<CT3*>(&eleGrp[eleIdx]);
        if (!t3Element) continue;

        vector<double> elementInfo;

        // 1. 添加单元编号（组内编号）
        elementInfo.push_back(eleIdx + 1);

        // 2. 添加节点编号
        elementInfo.push_back(t3Element->nodes_[0]->NodeNumber);  // 节点1
        elementInfo.push_back(t3Element->nodes_[1]->NodeNumber);  // 节点2
        elementInfo.push_back(t3Element->nodes_[2]->NodeNumber);  // 节点3

        // 3. 计算并添加应力分量
        double stress[3] = { 0 };  // 存储应力分量
        t3Element->ElementStress(stress, displacement);

        elementInfo.push_back(stress[0]);  // σ_x
        elementInfo.push_back(stress[1]);  // σ_y
        elementInfo.push_back(stress[2]);  // τ_xy

        t3ElementData.push_back(elementInfo);
    }

    return t3ElementData;
}


// 处理Q4单元组的函数
vector<vector<double>> getQ4ElementStressData(CElementGroup& eleGrp, double* displacement) {
    vector<vector<double>> q4ElementData;
    unsigned int numElements = eleGrp.GetNUME();

    for (unsigned int eleIdx = 0; eleIdx < numElements; eleIdx++) {
        // 获取当前Q4单元
        CQ4* q4Element = dynamic_cast<CQ4*>(&eleGrp[eleIdx]);
        if (!q4Element) continue;  // 确保转换成功

        vector<double> elementInfo;

        // 1. 添加单元编号（组内编号）
        elementInfo.push_back(eleIdx + 1);

        // 2. 添加节点编号
        elementInfo.push_back(q4Element->nodes_[0]->NodeNumber);  // 节点1
        elementInfo.push_back(q4Element->nodes_[1]->NodeNumber);  // 节点2
        elementInfo.push_back(q4Element->nodes_[2]->NodeNumber);  // 节点3
        elementInfo.push_back(q4Element->nodes_[3]->NodeNumber);  // 节点4

        // 3. 计算并添加应力分量（每个单元有4个高斯点，每个点3个应力分量）
        double stress[12] = { 0 };  // 存储12个应力分量（4个高斯点 × 3分量）
        q4Element->ElementStress(stress, displacement);

        // 添加所有应力分量
        for (int i = 0; i < 12; i++) {
            elementInfo.push_back(stress[i]);
        }

        q4ElementData.push_back(elementInfo);
    }

    return q4ElementData;
}







//T3单元应变在单元间不连续，需要通过周围单元的应力平均值来重构节点应力
void T3NodalStress(
    const vector<vector<double>> elementdata,
    vector<vector<double>>& nodestress
) {
    // 获取节点数量
    int nnp = nodestress.size();
    //获取单元数量
    int nel = elementdata.size();

    // 为每个节点添加位移放大
    for (int i = 0; i < nnp; i++) {
        // 添加节点号、X、Y坐标
        nodestress[i][1] += 50 * nodestress[i][4]); // X坐标
        nodestress[i][2] += 50 * nodestress[i][5]); // Y坐标
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
            nodestress[nodeIndex][8] += sigma_xx; // σ_xx
            nodestress[nodeIndex][9] += sigma_yy; // σ_yy
            nodestress[nodeIndex][10] += sigma_xy; // σ_xy
            nodestress[nodeIndex][7] += 1.0; // 计数位加1
        }
    }
}



//Q4单元为常应变单元，应变在单元间不连续，需要通过周围单元的最近的高斯点的应力平均值来重构节点应力
void Q4NodalStress(
    const vector<vector<double>> elementdata,
    vector<vector<double>>& nodestress
) {
    // 获取节点数量
    int nnp = nodestress.size();
    //获取单元数量
    int nel = elementdata.size();

    // 为每个节点添加位移放大
    for (int i = 0; i < nnp; i++) {
        // 添加节点号、X、Y坐标
        nodestress[i][1] += 50 * nodestress[i][4]); // X坐标
        nodestress[i][2] += 50 * nodestress[i][5]); // Y坐标
    }
    //遍历单元
    for (int j = 0; j < nel; j++) {
        // 对每个节点进行应力累加
        for (int k = 1; k < 5; k++) {
            int nodeIndex = nodeList[k] - 1; // 节点编号从1开始，转换为索引
            nodestress[nodeIndex][8] += elementdata[j][4 * k + 1]; // σ_xx
            nodestress[nodeIndex][9] += elementdata[j][4 * k + 2]; // σ_yy
            nodestress[nodeIndex][10] += elementdata[j][4 * k + 3]; // σ_xy
            nodestress[nodeIndex][7] += 1.0; // 计数位加1
        }
    }
}


