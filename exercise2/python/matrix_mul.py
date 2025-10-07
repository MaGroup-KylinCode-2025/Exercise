import numpy as np


def matrix_multiply(a, b):
    row = a.shape[0]
    col = b.shape[1]
    ld = a.shape[1]
    assert ld == b.shape[0]
    c = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            for k in range(ld):
                # c_ij += a_ik * b_kj
                c[i, j] += a[i, k] * b[k, j]
    return c


if __name__ == "__main__":
    a = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    b = a.copy()
    print(matrix_multiply(a, b))
    print(np.dot(a, b))
