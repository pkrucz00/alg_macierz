import numpy as np
from scipy.linalg import lu


def row_coord_to_matrix(A):
    result = np.zeros((len(A), len(A)))
    for y, row_dict in enumerate(A):
        for x, val in row_dict.items():
            result[x][y] = val
    return result


def matrix_to_coordinates(A):
    x_coords, y_coords = A.nonzero()
    vals = A[x_coords, y_coords]
    return list(zip(vals, x_coords, y_coords))


def coordinates_to_row_coord(coord_matrix, n):
    row_lists = [{} for _ in range(n)]
    for v, x, y in coord_matrix:
        row_lists[y][x] = v

    return row_lists


def matrix_to_row_coord(A):
    A_coord = matrix_to_coordinates(A)
    return coordinates_to_row_coord(A_coord, A.shape[0])


def dense_elimination(A):
    n = len(A)
    for k in range(n - 1):
        akk = A[k][k]
        for j in range(k+1, n):
            ajk = A[j][k]
            for i in range(k, n):
                A[j][i] -= (A[k][i]/akk)*ajk
    return A


def func(A, i, ajk, k, v):
    return v - A[k].get(i, 0)*ajk


def sparse_elimination(A):
    A = matrix_to_row_coord(A)
    n = len(A)
    for k, k_row_vals in enumerate(A):
        akk = k_row_vals[k]
        A[k] = {x: (v/akk if x >= k else v) for x, v in k_row_vals.items()}
        for j in range(k+1, n):
            j_row_vals = A[j]
            if k in j_row_vals:
                ajk = j_row_vals[k]
                A[j] = {i: (func(A, i, ajk, k, v) if i >= k else v) for i, v in j_row_vals.items()}
    return row_coord_to_matrix(A)


def get_matrix_from_csv(csv_file):
    return np.loadtxt(open(csv_file, "rb"), delimiter=",", skiprows=0)


def coord_to_graph(A_coord, n):
    graph = {i: set() for i in range(n)}
    for _, x, y in A_coord:
        if x != y:
            graph[x].add(y)

    return graph


def matrix_to_graph(A):
    return coord_to_graph(matrix_to_coordinates(A), A.shape[0])


def subtract_graphs(G, H):
    return {v: (G[v] - H[v] if v in H else G[v]) for v in G}


def get_fill_in_edges(A):
    def get_graph_without_edges_to_v(graph, v):
        return {u: adj_to_u - {v} for u, adj_to_u in graph.items()}

    n = A.shape[0]
    original_graph = matrix_to_graph(A)
    graph_copy = original_graph.copy()

    fill_ins = {v: set() for v in range(n)}
    for v in range(n):
        adj_to_v = graph_copy.pop(v)
        graph_copy = get_graph_without_edges_to_v(graph_copy, v)
        for u in adj_to_v:
            new_edges_adj_to_u = adj_to_v - {u}
            graph_copy[u] = graph_copy[u] | new_edges_adj_to_u
            fill_ins[u] = new_edges_adj_to_u

    return subtract_graphs(fill_ins, original_graph)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    A = get_matrix_from_csv("../matrices/matrix_0_18_2_0.csv")[:40, :40]
    print(get_fill_in_edges(A))


