import os
import re
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def run_stap(skip_existing=True):
    """运行STAPpp程序
    
    Args:
        skip_existing (bool): 如果为True，则跳过已存在对应.out文件的.dat文件
    """
    conv_dir = Path("convergence")
    
    # 获取所有.dat文件
    dat_files = sorted(conv_dir.glob("*.dat"))
    if not dat_files:
        print("未找到.dat文件！")
        return
        
    print(f"找到 {len(dat_files)} 个.dat文件")
    
    # 遍历运行每个文件
    for dat_file in dat_files:
        out_file = dat_file.with_suffix('.out')
        
        # 检查是否跳过已存在的文件
        if skip_existing and out_file.exists():
            print(f"跳过已存在的文件: {dat_file.name}")
            continue
            
        print(f"正在运行: {dat_file.name}")
        try:
            # 使用subprocess运行stap++命令，抑制输出
            result = subprocess.run(
                ["stap++", str(dat_file)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )
            print(f"成功完成: {dat_file.name}")
                
        except subprocess.CalledProcessError as e:
            print(f"运行失败: {dat_file.name}, 错误代码: {e.returncode}")
        except Exception as e:
            print(f"运行出错: {dat_file.name}")
            print(f"错误信息: {str(e)}")
    
def parse_output_file(filename):
    """解析.out文件获取节点位移和高斯点数据"""
    node_coords = []
    displacements = []
    gauss_data = []
    
    with open(filename) as f:
        lines = f.readlines()
        reading_coords = False
        reading_disp = False
        reading_gauss = False
        reading_gauss_coords = False
        reading_gauss_disp = False
        temp_gauss_coords = []
        temp_gauss_disp = []
        
        for line in lines:
            # 跳过空行
            if not line.strip():
                continue
            # 读取节点坐标
            if "N O D A L   P O I N T   D A T A" in line:
                reading_coords = True
                continue
            if reading_coords:
                if "EQUATION NUMBERS" in line:
                    reading_coords = False
                    continue
                if "NODE" in line:
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    x, y = float(parts[4]), float(parts[5])
                    node_coords.append([x, y])
                    
            # 读取位移
            if "D I S P L A C E M E N T S" in line:
                reading_disp = True
                continue
            if reading_disp:
                if "NODE" in line or not line.strip():  # 修改这里：忽略NODE行和空行
                    continue
                parts = line.split()
                if len(parts) >= 3:  # 修改这里：确保至少有3个字段
                    displacements.append([float(parts[1]), float(parts[2])])
                    if len(displacements) == len(node_coords):  # 如果已读取所有节点位移则退出
                        reading_disp = False

            # 读取高斯点数据
            if "S T R E S S  C A L C U L A T I O N S" in line:
                reading_gauss = True
                continue
                
            if reading_gauss:
                if "S O L U T I O N   T I M E   L O G" in line:
                    reading_gauss = False
                    reading_gauss_coords = False
                    reading_gauss_disp = False
                if "X" in line.strip() and "D_X" not in line.strip():
                    reading_gauss_coords = True
                    reading_gauss_disp = False
                    continue
                if "D_X" in line.strip():
                    reading_gauss_disp = True
                    reading_gauss_coords = False
                    continue
                if reading_gauss_coords:
                    coords = line.split()
                    for i in range(4):  # 4个高斯点
                        x = float(coords[i*3 + 1])
                        y = float(coords[i*3 + 2])
                        temp_gauss_coords.append((x, y))

                if reading_gauss_disp:
                    disp = line.split()
                    for i in range(4):  # 4个高斯点
                        ux = float(disp[i*3 + 1])
                        uy = float(disp[i*3 + 2])
                        temp_gauss_disp.append((ux, uy))

                    for coord, disp in zip(temp_gauss_coords, temp_gauss_disp):
                        gauss_data.append({
                            'coords': coord,
                            'disp': disp
                        })
                    temp_gauss_coords = []
                    temp_gauss_disp = []

    node_coords = np.array(node_coords)
    displacements = np.array(displacements)
    
    print(f"File: {filename}")
    print(f"  Number of nodes: {len(node_coords)}")
    print(f"  Number of displacements: {len(displacements)}")
    print(f"  Number of gauss points: {len(gauss_data)}")

    return node_coords, displacements, gauss_data

def compute_exact_solution(coords, config):
    """根据解析解表达式计算精确解"""
    exact_sol = config["exact_solution"]
    params = exact_sol["params"]
    P = params["P"]
    L = params["L"]
    E = config["materials"][0]["E"]
    nu = config["materials"][0]["nu"]
    
    x = coords[:, 0]
    y = coords[:, 1]
    
    # 计算精确解
    # 当y=0时位移为0（固定边界）
    dx = np.where(y > 0, -nu * P * (x - 0.5) / (E * L), 0.0)
    dy = np.where(y > 0, -P * y / (E * L), 0.0)
    
    return np.column_stack((dx, dy))

def compute_errors(gauss_data, config, h):
    """使用高斯点位移计算L2范数误差"""
    err_l2 = 0.0
    
    # 从配置文件读取参数
    exact_sol = config["exact_solution"]
    params = exact_sol["params"]
    P = params["P"]
    L = params["L"]
    E = config["materials"][0]["E"]
    nu = config["materials"][0]["nu"]

    for gp in gauss_data:
        x, y = gp['coords']
        ux_num, uy_num = gp['disp']
        
        # 计算精确解
        ux_exact = -nu * P * (x - 0.5) / (E * L) if y > 0 else 0.0
        uy_exact = -P * y / (E * L) if y > 0 else 0.0

        # 累加误差
        S = h ** 2 / 4
        err_l2 += S * ((ux_exact - ux_num) ** 2 + (uy_exact - uy_num) ** 2)
    
    return np.sqrt(err_l2)

def main():

    # 运行STAPpp程序
    run_stap(skip_existing=False)

    # 读取配置文件
    with open("config.json") as f:
        config = json.load(f)
    
    # 获取所有结果文件
    result_files = []
    h_values = []
    pattern = re.compile(r"model_(\d+)x\1.out")
    for f in Path("convergence").glob("*.out"):
        if m := pattern.match(f.name):
            n = int(m.group(1))
            h = 1.0/n  # 假设域的大小为1
            result_files.append(f)
            h_values.append(h)
    
    # 排序确保h值递减
    result_files = [x for _, x in sorted(zip(h_values, result_files), reverse=True)]
    h_values = sorted(h_values, reverse=True)
    
    # 计算误差
    err_l2 = []
    
    for f, h in zip(result_files, h_values):
        # 读取数值解和节点坐标
        coords, disp, gauss_data = parse_output_file(f)
        
        # 使用高斯积分计算L2误差
        e_l2 = compute_errors(gauss_data, config, h)
        err_l2.append(e_l2)
        print(f"h = {h}:")
        print(f"  Nodes: {len(coords)}")
        print(f"  L2 error: {e_l2}")
    
    # 绘制收敛曲线
    plt.figure(figsize=(10, 8))
    
    # 计算对数值
    log_h = np.log(h_values)
    log_err = np.log(err_l2)
    
    # 线性拟合得到收敛率
    slope, intercept = np.polyfit(log_h, log_err, 1)
    fit_err = np.exp(intercept + slope * log_h)
    
    # 绘制误差曲线和拟合线
    plt.loglog(h_values, err_l2, 'bo-', label='L2 norm error')
    plt.loglog(h_values, fit_err, 'r--', 
              label=f'Fitted line (rate = {slope:.2f})')
    
    # 参考线
    # h_ref = np.array([min(h_values), max(h_values)])
    # plt.loglog(h_ref, h_ref**2 * err_l2[-1]/h_values[-1]**2, 
    #           'k:', label='O(h²)')
    
    plt.grid(True)
    plt.xlabel('Element size (h)')
    plt.ylabel('Error')
    plt.legend()
    plt.title('Convergence Analysis of Q4 Element')
    plt.savefig('convergence.png')
    plt.show()
    plt.close()
    
    # 输出整体收敛率
    print("\nConvergence Analysis:")
    print(f"Overall convergence rate: {abs(slope):.2f}")
    print(f"R² value: {np.corrcoef(log_h, log_err)[0,1]**2:.4f}")

if __name__ == "__main__":
    main()
