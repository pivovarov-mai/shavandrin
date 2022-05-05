import numpy as np
from numpy.linalg import norm
from tabulate import tabulate


def householder(matrix):
    n = matrix.shape[0]
    Q = np.identity(n)
    R = np.copy(matrix)
    for i in range(n - 1):
        x = R[i:, i]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        Q_i = np.identity(n)
        Q_i[i:, i:] -= 2.0 * np.outer(v, v)
        R = np.dot(Q_i, R)
        Q = np.dot(Q, Q_i)
    return Q, R


def find_roots(A, i):
    n = len(A)
    a11 = A[i][i]
    a12 = A[i][i + 1] if i + 1 < n else 0
    a21 = A[i + 1][i] if i + 1 < n else 0
    a22 = A[i + 1][i + 1] if i + 1 < n else 0
    return np.roots((1, -a11 - a22, a11 * a22 - a12 * a21))


def check_complex(A, eps, i):
    Q, R = householder(A)
    A_next = np.dot(R, Q)
    lambda1 = find_roots(A, i)
    lambda2 = find_roots(A_next, i)
    return True if abs(lambda1[0] - lambda2[0]) <= eps and \
                   abs(lambda1[1] - lambda2[1]) <= eps else False


def get_eigenvalue(a, eps, i):
    A = np.copy(a)
    while True:
        Q, R = householder(A)
        A = np.dot(R, Q)
        if norm(A[i + 1:, i]) <= eps:
            return A[i][i], "real", A
        elif norm(A[i + 2:, i]) <= eps and check_complex(A, eps, i):
            return find_roots(A, i), "complex", A


def QR_algorithm(a, eps):
    res = []
    i = 0
    a = np.copy(a)
    while i < len(a):
        res_i, type_, a = get_eigenvalue(a, eps, i)
        if type_ == "complex":
            for r in res_i:
                res.append(complex(r))
            i += 2
        elif type_ == "real":
            res.append(res_i)
            i += 1
    return np.array(res)


def main():
    np.set_printoptions(precision=4)
    n = int(input())
    e = float(input())
    a = np.array([list(map(float, input().split())) for _ in range(n)])

    print("Matrix A:")
    print(tabulate(a), '\n')
    
    q, r = householder(a)

    res = QR_algorithm(a, e)
    print("Eigen values (QR):")
    print(res, '\n')
    
    print("Eigen values (numpy):")
    v, vectors = np.linalg.eig(a)
    print(v)


if __name__ == "__main__":
    main()