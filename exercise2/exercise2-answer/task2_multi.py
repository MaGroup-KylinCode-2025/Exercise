import numpy as np
import time

def matrix_multiply(a, b):
    row = a.shape[0]
    col = b.shape[1]
    ld = a.shape[1]
    assert ld == b.shape[0], "Matrix dimensions do not match for multiplication"
    c = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            for k in range(ld):
                c[i, j] += a[i, k] * b[k, j]
    return c

if __name__ == "__main__":
    # Define matrices
    a = np.array([[1, 2], [3, 4]])
    b = np.array([[5, 6], [7, 8]])

    # Measure time for nested loop multiplication
    start_time = time.time()
    result_nested = matrix_multiply(a, b)
    nested_time = time.time() - start_time
    print("Nested loop result:")
    print(result_nested)
    print(f"Nested loop time: {nested_time:.6f} seconds")

    # Measure time for numpy.dot multiplication
    start_time = time.time()
    result_numpy = np.dot(a, b)
    numpy_time = time.time() - start_time
    print("NumPy dot result:")
    print(result_numpy)
    print(f"NumPy dot time: {numpy_time:.6f} seconds")