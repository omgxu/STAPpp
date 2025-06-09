import matplotlib.pyplot as plt

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

# 绘制
plt.figure(figsize=(6,6))
for elem in elements:
    pts = [node_coords[n] for n in elem] + [node_coords[elem[0]]]  # 闭合三角形
    xs, ys = zip(*pts)
    plt.plot(xs, ys, 'k-')
    plt.fill(xs, ys, edgecolor='k', fill=False)
    # 标注单元号
    cx = sum(xs[:-1])/3
    cy = sum(ys[:-1])/3
    plt.text(cx, cy, str(elements.index(elem)+1), color='red', fontsize=15, ha='center', va='center')

# 标注节点
for num, (x, y) in node_coords.items():
    plt.plot(x, y, 'bo')
    plt.text(x, y, str(num), color='blue', fontsize=9, ha='right', va='bottom')

plt.axis('equal')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('T3 Patch Elements')
plt.show()