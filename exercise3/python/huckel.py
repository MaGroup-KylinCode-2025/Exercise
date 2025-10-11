import numpy as np

def huckel_mo(H, name="Molecule"):
    """
    Huckel分子轨道计算
    
    参数:
        H: Hamiltonian矩阵 (numpy array)
        name: 分子名称
    """
    print(f"\n{'='*40}")
    print(f"Huckel MO Calculation: {name}")
    print(f"{'='*40}")
    
    print("\nHamiltonian Matrix H:")
    print(H)
    
    # 对角化
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    print("\nEigenvalues (Orbital Energies):")
    print(eigenvalues)
    
    print("\nEigenvectors (MO Coefficients):")
    print(eigenvectors)
    
    return eigenvalues, eigenvectors


# ========================================
# 示例1: 丁二烯 (Butadiene)
# ========================================
H_butadiene = np.array([
    [ 0.0, -1.0,  0.0,  0.0],
    [-1.0,  0.0, -1.0,  0.0],
    [ 0.0, -1.0,  0.0, -1.0],
    [ 0.0,  0.0, -1.0,  0.0]
])

e, c = huckel_mo(H_butadiene, "Butadiene")


# ========================================
# 示例2: 苯 (Benzene)
# ========================================
H_benzene = np.array([
    [ 0.0, -1.0,  0.0,  0.0,  0.0, -1.0],
    [-1.0,  0.0, -1.0,  0.0,  0.0,  0.0],
    [ 0.0, -1.0,  0.0, -1.0,  0.0,  0.0],
    [ 0.0,  0.0, -1.0,  0.0, -1.0,  0.0],
    [ 0.0,  0.0,  0.0, -1.0,  0.0, -1.0],
    [-1.0,  0.0,  0.0,  0.0, -1.0,  0.0]
])

e, c = huckel_mo(H_benzene, "Benzene")
