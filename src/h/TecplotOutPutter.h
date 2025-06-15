

////时间原因尚未完成，有待进一步完善////
#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

// 进行应力重构
void NodalStress(

);
void BarNodalStress(
    const vector<vector<double>> nodes,
    const vector<vector<double>> gaussPoints,
    vector<double>& nodeStresses
);
void Q4NodalStress(
    vector<vector<double>> nodedata,
    vector<vector<double>> elementdata,
    vector<double>& nodestress
);
void T3NodalStress(
    const vector<vector<double>> nodes,
    const vector<vector<double>> gaussPoints,
    vector<double>& nodeStresses
);
void H8NodalStress(
    const vector<vector<double>> nodes,
    const vector<vector<double>> gaussPoints,
    vector<double>& nodeStresses
);
void BeamRebuildEquation(
    const vector<vector<double>> nodes,
    const vector<vector<double>> gaussPoints,
    vector<double>& nodeStresses
);
void Q9NodalStress(
    const vector<vector<double>> nodes,
    const vector<vector<double>> gaussPoints,
    vector<double>& nodeStresses
);
// 其他成员函数可以在这里声明


vector<vector<double>> getT3ElementStressData();

vector<vector<double>> getNodalData();