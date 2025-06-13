# 绘制T3单元的误差范数曲线

import re
import matplotlib.pyplot as plt
import numpy as np

def gauss(ngp):
	"""
	Get Gauss points in the parent element domain [-1, 1] and
	the corresponding weights.

	Args:
		ngp : (int) number of Gauss points.

	Returns: w,gp
		w  : weights.
		gp : Gauss points in the parent element domain.
	"""
	gp = None
	w = None
	if ngp == 1:
		gp = [0]
		w = [2]
	elif ngp == 2:
		gp = [-0.57735027, 0.57735027]
		w = [1, 1]
	elif ngp == 3:
		gp = [-0.7745966692, 0.7745966692, 0.0]
		w = [0.5555555556, 0.5555555556, 0.8888888889]
	else:
		raise ValueError("The given number (ngp = {}) of Gauss points is too large and not implemented".format(ngp))
	return w, gp

def NmatElastT3(eta, psi):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	N1 = 0.25*(1-psi)*(1-eta)
	N2 = 0.25*(1+psi)*(1-eta)
	N3 = 0.5*(1+eta)

	return np.array([[N1, 0, N2, 0, N3, 0],
					 [0, N1, 0, N2, 0, N3]])

def BmatElastT3(eta, psi, C):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate
		C   : The physical coordinates

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	#Calculate the Grad(N) matrix
	GN = 0.25*np.array([[eta-1, 1-eta, 0],
						[psi-1, -psi-1,2]])

	# Compute Jacobian matrix
	J = GN@C
	detJ = np.linalg.det(J)

	BB = np.linalg.solve(J, GN)
	B1x = BB[0, 0]
	B2x = BB[0, 1]
	B3x = BB[0, 2]
	B1y = BB[1, 0]
	B2y = BB[1, 1]
	B3y = BB[1, 2]

	B = np.array([[B1x, 0, B2x, 0, B3x, 0],
				  [0, B1y, 0, B2y, 0, B3y],
				  [B1y, B1x, B2y, B2x, B3y, B3x]])

	return B, detJ

def extract_coordinate(filename):
    coordinate = []
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 找到 DISPLACEMENTS 标记
    for idx, line in enumerate(lines):
        if 'NODALPOINTDATA' in line.replace(" ", "").upper():
            print(f"Found N O D A L   P O I N T   D A T A at line {idx}")
            break
    else:
        return coordinate

    # 向下找表头 "NODE..." 所在行
    for start in range(idx + 1, len(lines)):
        if 'NUMBER' in lines[start] and 'CONDITION' in lines[start]:
            start += 1  # 数据从下一行开始
            break
    else:
        return coordinate

    # 读取数据
    for line in lines[start:]:
        print(f"DEBUG LINE: {repr(line)}")
        if not line.strip():
            print("Blank line — stopping.")
            break
        if not re.match(r'\s*\d+', line):
            print("Not a data line — stopping.")
            break
        parts = line.split()
        if len(parts) < 7:
            print(f"Line skipped — too few parts: {parts}")
            continue
        try:
            x = np.double(parts[4])
            y = np.double(parts[5])
            z = np.double(parts[6])
            coordinate.append((x, y, z))
        except Exception as e:
            print(f"Conversion error: {e} — line: {parts}")
            continue
    print(f"Extracted coordinate: {coordinate}")
    return coordinate

def extract_displacements(filename):
    displacements = []
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 找到 DISPLACEMENTS 标记
    for idx, line in enumerate(lines):
        if 'DISPLACEMENTS' in line.replace(" ", "").upper():
            print(f"Found D I S P L A C E M E N T S at line {idx}")
            break
    else:
        return displacements

    # 向下找表头 "NODE..." 所在行
    for start in range(idx + 1, len(lines)):
        if 'NODE' in lines[start] and 'DISPLACEMENT' in lines[start]:
            start += 1  # 数据从下一行开始
            break
    else:
        return displacements

    # 读取数据
    for line in lines[start:]:
        print(f"DEBUG LINE: {repr(line)}")
        if not line.strip():
            print("Blank line — stopping.")
            break
        if not re.match(r'\s*\d+', line):
            print("Not a data line — stopping.")
            break
        parts = line.split()
        if len(parts) < 4:
            print(f"Line skipped — too few parts: {parts}")
            continue
        try:
            node = int(parts[0])
            ux = float(parts[1])
            uy = float(parts[2])
            uz = float(parts[3])
            displacements.append((node, ux, uy, uz))
        except Exception as e:
            print(f"Conversion error: {e} — line: {parts}")
            continue

    return displacements

def extract_IEM(filename):
    IEN = []
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 找到 DISPLACEMENTS 标记
    for idx, line in enumerate(lines):
        if 'ELEMENTINFORMATION' in line.replace(" ", "").upper():
            print(f"Found E L E M E N T   I N F O R M A T I O N at line {idx}")
            break
    else:
        return IEN

    # 向下找表头 "NODE..." 所在行
    for start in range(idx + 1, len(lines)):
        if 'NUMBER-N' in lines[start] and 'I' in lines[start]:
            start += 1  # 数据从下一行开始
            break
    else:
        return IEN

    # 读取数据
    for line in lines[start:]:
        print(f"DEBUG LINE: {repr(line)}")
        if not line.strip():
            print("Blank line — stopping.")
            break
        if not re.match(r'\s*\d+', line):
            print("Not a data line — stopping.")
            break
        parts = line.split()
        if len(parts) < 5:
            print(f"Line skipped — too few parts: {parts}")
            continue
        try:
            nnod1 = int(parts[1])
            nnod2 = int(parts[2])
            nnod3 = int(parts[3])
            IEN.append((nnod1, nnod2, nnod3))
        except Exception as e:
            print(f"Conversion error: {e} — line: {parts}")
            continue
    print(f"Extracted IEN: {IEN}")
    return IEN

def ErrorNorm_T3(filename,i):
    n = 2**i # 纵向的单元个数，整个分为n×4n的矩形
    nel= 8*(n**2)

    ngp = 2 # 高斯积分点个数
    [w,gp] = gauss(ngp) # 高斯积分点和权重

    L2Norm = 0
    EnNorm = 0

    L2NormEx = 0
    EnNormEx = 0

    F = 1 # 载荷
    E = 1000 # 弹性模量
    I = 1/12 # 惯性矩
    l = 6 # 单元长度
    
    coordinate = extract_coordinate(filename)
    if not coordinate:
        print("No coordinate data found.")
    disp = extract_displacements(filename)
    if not disp:
            print("No displacement data found.")
    IEN = extract_IEM(filename)
    if not IEN:
            print("No IEN data found.")
    disp_x = [ux for _, ux, _, _ in disp]
    disp_y = [uy for _, _, uy, _ in disp]
    
    # 整个位移场的误差范数，每个单元计算
    for e in range(nel):
         
        C = np.zeros((3, 2)) # 单元3节点的物理坐标矩阵
        C[0, 0] = coordinate[IEN[e][0]-1][0]
        C[0, 1] = coordinate[IEN[e][0]-1][1]
        C[1, 0] = coordinate[IEN[e][1]-1][0]
        C[1, 1] = coordinate[IEN[e][1]-1][1]
        C[2, 0] = coordinate[IEN[e][2]-1][0]
        C[2, 1] = coordinate[IEN[e][2]-1][1]

        xye = np.zeros(6) # 单元3节点的物理坐标列阵
        xye[0] = C[0, 0]
        xye[1] = C[0, 1]
        xye[2] = C[1, 0]
        xye[3] = C[1, 1]
        xye[4] = C[2, 0]
        xye[5] = C[2, 1]
        
        de = np.zeros(6) # 单元3节点的有限元位移
        de[0] = disp_x[IEN[e][0]-1]
        de[1] = disp_y[IEN[e][0]-1]
        de[2] = disp_x[IEN[e][1]-1]
        de[3] = disp_y[IEN[e][1]-1]
        de[4] = disp_x[IEN[e][2]-1]
        de[5] = disp_y[IEN[e][2]-1]

        for i in range(ngp):
            for j in range(ngp):
                eta = gp[i]
                psi = gp[j]

            # derivative of the shape functions
            N= NmatElastT3(eta, psi)
            B, detJ = BmatElastT3(eta, psi, C)

            uh = N @ de # 高斯点位移
            xyt = N @ xye # 高斯点物理坐标

            uxex = 0 # 理论x位移
            uyex = F*(xyt[0]**2)/(6*E*I)*(3*l - xyt[0]) # 理论y位移
            uex = np.array([uxex, uyex]) # 理论位移

            L2Norm += detJ*w[i]*w[j]*(np.dot(uex-uh,uex-uh))

            sh = B@ de # 高斯点应变
            sex = np.array([0,0,F*l/(2*E*I)*(2*xyt[0]-xyt[0]**2)]) # 应变列阵
            
            EnNorm += 0.5*detJ*w[i]*w[j]*E*(np.dot(sex-sh,sex-sh))

    L2Norm = np.sqrt(L2Norm)
    EnNorm = np.sqrt(EnNorm)
    return 1.0/n,L2Norm, EnNorm


if __name__ == "__main__":
    files = ("T3mesh-1.out","T3mesh-2.out","T3mesh-4.out","T3mesh-8.out")
    num_files = len(files)
    h = np.zeros(num_files)
    L2Norm = np.zeros(num_files)
    EnNorm = np.zeros(num_files)
    i = 0
    if num_files == 0:
        print("No T3mesh-*.out files found.")
        exit(1)
    for file in files:
        print(f"File: {file}")
        full_file_path = "data/T3/" + file
        h[i], L2Norm[i], EnNorm[i] = ErrorNorm_T3(full_file_path,i) 
        # for node, ux, uy, uz in disp:
        #     print(f"Node {node:3d}: Ux={ux:.6e}, Uy={uy:.6e}, Uz={uz:.6e}")
        print("-" * 40)
        i += 1

    logh=np.log10(h)
    logL2Norm=np.log10(L2Norm)
    logEnNorm=np.log10(EnNorm)

    plt.plot(logh, logL2Norm, marker='o', linestyle='-')
    plt.xlabel('Log of Element Length')
    plt.ylabel('Log of L2-norm Error')
    plt.title('Log-Log Plot of L2-norm Error & Element Length')
    plt.grid(True)
    plt.show()

    plt.plot(logh, logEnNorm, marker='o', linestyle='-')
    plt.xlabel('Log of Element Length')
    plt.ylabel('Log of Energy-norm Error')
    plt.title('Log-Log Plot of Energy-norm Error & Element Length')
    plt.grid(True)
    plt.show()

    # plot the element length - error norm curve in logarithmic scale
    fig,(axs) = plt.subplots(1,2)
    plt.tight_layout()

    axs[0].set_title('Q4 element', fontsize=9); 
    axs[0].set_ylabel('L_2 error', fontsize=8)
    axs[0].xaxis.set_tick_params(labelsize=7)
    axs[0].yaxis.set_tick_params(labelsize=7)
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].plot(h,L2Norm)


    axs[1].set_xlabel('Element length (m)', fontsize=8); 
    axs[1].set_ylabel('Energy error', fontsize=8)
    axs[1].xaxis.set_tick_params(labelsize=7)
    axs[1].yaxis.set_tick_params(labelsize=7)
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].plot(h,EnNorm)



    # Print error norms obtained by the linear element and quadratic element
    #    with different element size
    print("\nError norms of linear elements")
    print('%13s %13s %13s' %('h','L2Norm','EnNorm'))
    for i in range(len(h)):
        print('%13.6E %13.6E %13.6E' %(h[i], L2Norm[i], EnNorm[i]))

    # print("\nSlope of L2-norm error:")
    # print('%13s %13s %13s' %('Slope1','Slope2','Slope3'))
    # slope = np.zeros(len(h)-1)
    # for i in range(len(h)-1):
    #     slope[i] = (np.log(L2Norm[i+1]) - np.log(L2Norm[i])) / (np.log(h[i+1]) - np.log(h[i]))
    # print('%13.6E %13.6E %13.6E' %(slope[0], slope[1], slope[2]))

    # print("\nSlope of energy-norm error:")
    # print('%13s %13s %13s' %('Slope1','Slope2','Slope3'))
    # for i in range(len(h)-1):
    #     slope[i] = (np.log(EnNorm[i+1]) - np.log(EnNorm[i])) / (np.log(h[i+1]) - np.log(h[i]))
    # print('%13.6E %13.6E %13.6E' %(slope[0], slope[1], slope[2]))

    print("\nSlope of L2-norm error:")
    header = ''.join([f"{'Slope'+str(i+1):>13s}" for i in range(len(h)-1)])
    print(header)
    slope = np.zeros(len(h)-1)
    for i in range(len(h)-1):
        slope[i] = (np.log(L2Norm[i+1]) - np.log(L2Norm[i])) / (np.log(h[i+1]) - np.log(h[i]))
    print(''.join([f"{s:13.6E}" for s in slope]))

    print("\nSlope of energy-norm error:")
    print(header)
    for i in range(len(h)-1):
        slope[i] = (np.log(EnNorm[i+1]) - np.log(EnNorm[i])) / (np.log(h[i+1]) - np.log(h[i]))
    print(''.join([f"{s:13.6E}" for s in slope]))