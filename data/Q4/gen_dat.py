import os
import json
import numpy as np

def create_q4_element_data(nx, ny, material_group=1):
    """创建Q4单元数据"""
    elements = []

    for j in range(ny):
        for i in range(nx):
            node_bl = i + j * (nx + 1) + 1
            node_br = node_bl + 1
            node_tr = node_bl + nx + 2
            node_tl = node_bl + nx + 1

            element_number = i + j * nx + 1
            elements.append((element_number, node_bl, node_br, node_tr, node_tl, material_group))

    return elements

def generate_q4_model_file(filename, nx, ny, materials):
    """
    生成单个 Q4 模型文件
    :param filename: 输出文件名
    :param nx: x方向单元数
    :param ny: y方向单元数
    :param materials: 材料列表 [(id, E, nu), ...]
    """

    title = f"Q4 Convergence Model {nx}x{ny}"
    num_nodes_x = nx + 1
    num_nodes_y = ny + 1
    total_nodes = num_nodes_x * num_nodes_y
    total_element_groups = 1
    load_cases = 1
    solve_mode = 1

    nodes = []
    for j in range(num_nodes_y):
        for i in range(num_nodes_x):
            node_number = j * num_nodes_x + i + 1
            x_coord = i / nx
            y_coord = j / ny

            if j == 0:
                bc_x, bc_y = 1, 1
            else:
                bc_x, bc_y = 0, 0
            bc_z = 1

            nodes.append((node_number, bc_x, bc_y, bc_z, x_coord, y_coord, 0.0))

    elements = create_q4_element_data(nx, ny)

    # 修改载荷部分
    load_control_line = (1, 1)
    P = -1000.0  # 边界节点上的总力
    top_node_ids = [num_nodes_x * (num_nodes_y - 1) + i + 1 for i in range(num_nodes_x)]
    loads = []
    
    # 两端节点分配一半，其余节点分配全量（常规均布压力等效节点力分配）
    for idx, node_number in enumerate(top_node_ids):
        if idx == 0 or idx == num_nodes_x - 1:
            force = P / (2 * (num_nodes_x - 1))
        else:
            force = P / (num_nodes_x - 1)
        loads.append((node_number, 2, force))
    
    # 更新载荷控制行
    load_control_line = (1, len(loads))  # 更新载荷数量

    with open(filename, 'w') as f:
        # 标题行
        f.write(f"{title[:80]}\n")

        # 控制行
        f.write(f"{total_nodes:5d}{total_element_groups:5d}{load_cases:5d}{solve_mode:5d}\n")

        # 节点数据
        for node in nodes:
            f.write(f"{node[0]:5d}{node[1]:5d}{node[2]:5d}{node[3]:5d}")
            f.write(f"{node[4]:10.5f}{node[5]:10.5f}{node[6]:10.5f}\n")

        # 载荷数据
        f.write(f"{load_control_line[0]:5d}{load_control_line[1]:5d}\n")
        for load in loads:
            f.write(f"{load[0]:5d}{load[1]:5d}{load[2]:10.1f}\n")

        # 单元组控制数据
        total_elements = len(elements)
        num_material_groups = len(materials)
        f.write(f"{2:5d}{total_elements:5d}{num_material_groups:5d}\n")  # Q4单元类型为2

        # 材料/截面性质数据
        for mat in materials:
            f.write(f"{mat[0]:5d}{mat[1]:10.3e}{mat[2]:10.3f}\n")

        # 单元数据
        for elem in elements:
            f.write(f"{elem[0]:5d}{elem[1]:5d}{elem[2]:5d}{elem[3]:5d}{elem[4]:5d}{elem[5]:5d}\n")

        # 结束标记
        f.write("stop\n")

def main():
    # 读取配置文件
    with open('config.json', 'r') as f:
        config = json.load(f)

    materials = [
        (mat['id'], mat['E'], mat['nu']) for mat in config['materials']
    ]
    mesh_sizes = [
        (size['nx'], size['ny']) for size in config['mesh_sizes']
    ]

    exact_solution = config.get('exact_solution', {})

    output_dir = 'convergence'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for nx, ny in mesh_sizes:
        filename = os.path.join(output_dir, f'model_{nx}x{ny}.dat')
        generate_q4_model_file(filename, nx, ny, materials)

    # 可选：处理精确解表达式
    if exact_solution.get('enabled', False):
        expr_dx = exact_solution.get('dx', '')
        expr_dy = exact_solution.get('dy', '')
        print(f"[INFO] 精确解表达式已提供：dx = {expr_dx}, dy = {expr_dy}")
        # 后续可用于误差分析、收敛性对比等

if __name__ == '__main__':
    main()