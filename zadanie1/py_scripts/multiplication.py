import numpy as np


def multi_ijp(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for i in range(m):
        for j in range(n):
            for p in range(k):
                C[i][j] += A[i][p] * B[p][j]


def multi_ipj(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for i in range(m):
        for p in range(k):
            for j in range(n):
                C[i][j] += A[i][p] * B[p][j]


def multi_jip(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for j in range(n):
        for i in range(m):
            for p in range(k):
                C[i][j] += A[i][p] * B[p][j]


def multi_jpi(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for j in range(n):
        for p in range(k):
            for i in range(m):
                C[i][j] += A[i][p] * B[p][j]


def multi_pij(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for p in range(k):
        for i in range(m):
            for j in range(n):
                C[i][j] += A[i][p] * B[p][j]


def multi_pji(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for p in range(k):
        for j in range(n):
            for i in range(m):
                C[i][j] += A[i][p] * B[p][j]


def mat_mul(A, B, mul_func=multi_ijp):
    C = np.zeros((np.size(A, 0), np.size(B, 1)))
    mul_func(A, B, C)
    return C


A = np.array(
    [[1, 2, 2],
     [2, 4, 3]])
B = np.array(
    [[2, 2, 2, 1],
     [2, 3, 1, 4],
     [2, 1, 3, 7]])
C = mat_mul(A, B)
print(C)

C_expected = A @ B
print(C_expected)
