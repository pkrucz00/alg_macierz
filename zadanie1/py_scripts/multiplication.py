import numpy as np

def mat_mul(A, B, mul_func):
    C = np.zeros(np.size(A, 0), np.size(B, 1))
    mul_func(A, B, C)
    return C

def multi_ijp(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for i in range(m):
        for j in range(n):
            for p in range(k):
                C[i][j] = A[i][p] * B[p][j]