import numpy as np
import argparse
import os

def generate_T3_mesh(nx, ny, width=6.0, height=1.0):
    """
    生成矩形区域的T3三角形单元网格
    nx, ny: x和y方向的节点数
    width, height: 矩形的宽度和高度
    """
    
    # 生成节点坐标
    x = np.linspace(0, width, nx+1)
    y = np.linspace(0, height, ny+1)
    nodes = []
    node_id = 1
    
    for j in range(ny+1):
        for i in range(nx+1):
            # 节点格式：编号 约束x 约束y 约束z x坐标 y坐标 z坐标
            constraint_x = 1 if i == 0 else 0  # 左边界固定
            constraint_y = 1 if i == 0 else 0
            nodes.append(f"{node_id}   {constraint_x}   {constraint_y}   1   {x[i]:.3f}   {y[j]:.3f}   0")
            node_id += 1

    # 生成单元
    elements = []
    elem_id = 1
    for j in range(ny):
        for i in range(nx):
            # 每个矩形分成两个三角形
            n1 = j * (nx+1) + i + 1
            n2 = n1 + 1
            n3 = n1 + nx + 1
            n4 = n2 + nx + 1
            
            # 第一个三角形
            elements.append(f"{elem_id}   {n1}   {n2}   {n3}   1")
            elem_id += 1
            
            # 第二个三角形
            elements.append(f"{elem_id}   {n2}   {n4}   {n3}   1")
            elem_id += 1

    filename = fr"C:\Users\ASUS\Desktop\ccode\STAPpp_git\STAPpp\src\data\T3\T3mesh-{ny}.dat"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    # 生成完整的输入文件
    with open(filename, 'w') as f:
        f.write(f"Generated T3 mesh - {ny} times {nx}\n\n")
        
        # 写入控制参数：节点数 材料数 约束数 载荷数
        f.write(f"{len(nodes)}   1   1   1\n\n")
        
        # 写入节点信息
        f.write("\n".join(nodes) + "\n\n")
        
        # 写入载荷信息
        n_leftup = (ny + 1) * (nx + 1)
        f.write("1   1\n")  # 载荷类型和载荷点数
        f.write(f"{n_leftup}   2   -1\n\n")  # 右边界施加载荷
        
        # 写入单元信息
        f.write(f"3   {len(elements)}   1\n")  # 单元类型 单元数 材料数
        f.write("1   1e3   0.3\n")  # 材料属性：E nu
        f.write("\n".join(elements))

if __name__ == "__main__":
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='生成T3三角形单元网格')
    parser.add_argument('--nx', type=int, default=4, help='x方向单元数')
    parser.add_argument('--ny', type=int, default=4, help='y方向单元数')
    parser.add_argument('--width', type=float, default=6.0, help='矩形宽度')
    parser.add_argument('--height', type=float, default=1.0, help='矩形高度')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 使用参数生成网格
    generate_T3_mesh(nx=args.nx, ny=args.ny, width=args.width, height=args.height)
    print(f"已生成网格文件 ''T3mesh-{args.ny}.dat'，网格大小：{args.nx}x{args.ny}")