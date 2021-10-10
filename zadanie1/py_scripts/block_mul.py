import numpy as np

def multi_ijp(A, B, C):
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)
    for i in range(m):
        for j in range(n):
            for p in range(k):
                C[i][j] += A[i][p] * B[p][j]


def block_mul(A, B, m_bs, n_bs, k_bs):
    C = np.zeros((np.size(A, 0), np.size(B, 1)))
    m, n, k = np.size(C, 0), np.size(C, 1), np.size(A, 1)

    for i in range(0, m, m_bs):
        for j in range(0, n, n_bs):
            for p in range(0, k, k_bs):
                A_slice = A[i:i+m_bs, p:p+k_bs]
                B_slice = B[p:p+k_bs, j:j+n_bs]
                C_slice = C[i:i+m_bs, j:j+n_bs]
                multi_ijp(A_slice, B_slice, C_slice)
    return C

A = np.array(
    [[1, 2, 2],
     [2, 4, 3]])
B = np.array(
    [[2, 2, 2, 1],
     [2, 3, 1, 4],
     [2, 1, 3, 7]])
C = block_mul(A, B, 2, 2, 2)
print(C)


C_expected = A @ B
print(C_expected)

