# 绘制分片试验图像

import matplotlib.pyplot as plt
import numpy as np

# 节点位移数据（节点号: [dx, dy]）
displacement = {
    1: [0.00000e+00,  0.00000e+00],
    2: [4.98733e-18, -3.00000e-03],
    3: [1.00000e-02, -3.00000e-03],
    4: [1.00000e-02,  0.00000e+00],
    5: [3.00000e-03, -1.50000e-03],
    6: [8.00000e-03, -1.80000e-03],
}

node_coords = {}
elements = []

with open(r'c:\Users\ASUS\Desktop\ccode\STAPpp_git\STAPpp\src\data\T3\T3patch.dat', encoding='utf-8') as f:
    lines = f.readlines()

# 读取节点
for line in lines:
    parts = line.strip().split()
    if len(parts) == 7 and parts[0].isdigit():
        node_num = int(parts[0])
        x = float(parts[4])
        y = float(parts[5])
        node_coords[node_num] = (x, y)

# 读取单元
for line in lines:
    parts = line.strip().split()
    if len(parts) == 5 and parts[0].isdigit():
        n1 = int(parts[1])
        n2 = int(parts[2])
        n3 = int(parts[3])
        elements.append((n1, n2, n3))

plt.figure(figsize=(8,8))

# 绘制变形前单元
for elem in elements:
    pts = [node_coords[n] for n in elem] + [node_coords[elem[0]]]
    xs, ys = zip(*pts)
    plt.plot(xs, ys, 'k-', label='Original' if elem==elements[0] else "")
    # 标注单元号
    cx = sum(xs[:-1])/3
    cy = sum(ys[:-1])/3
    plt.text(cx, cy, str(elements.index(elem)+1), color='red', fontsize=12, ha='center', va='center')

# 绘制变形后单元（放大变形10倍以便观察）
scale = 20
for elem in elements:
    pts = []
    for n in elem:
        x0, y0 = node_coords[n]
        dx, dy = displacement[n]
        pts.append((x0 + dx*scale, y0 + dy*scale))
    pts.append(pts[0])
    xs, ys = zip(*pts)
    plt.plot(xs, ys, 'b--', label='Deformed (x{} exaggerated)'.format(scale) if elem==elements[0] else "")

# 标注节点
for num, (x, y) in node_coords.items():
    plt.plot(x, y, 'ko')
    plt.text(x, y, str(num), color='black', fontsize=9, ha='right', va='bottom')
    # 变形后节点
    dx, dy = displacement[num]
    plt.plot(x+dx*scale, y+dy*scale, 'bo')
    plt.text(x+dx*scale, y+dy*scale, str(num), color='blue', fontsize=9, ha='left', va='top')

plt.axis('equal')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('T3 Patch Elements: Original (black) & Deformed (blue, x{} exaggerated)'.format(scale))
plt.legend()
plt.show()