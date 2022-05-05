import numpy as np
import math
import copy
from tabulate import tabulate


def finish_iterations(a, x_old, x_new, eps):
    a_norm = np.linalg.norm(a, ord=np.inf) 
    v_norm = np.linalg.norm(x_new - x_old, ord=np.inf)

    if a_norm >= 1:
        ans = v_norm
    else:
        ans = a_norm / (1 - a_norm) * v_norm
    return ans <= eps


def iterative(a, b, eps=0.001):
    n = len(a)
    count = len(b)
    x = np.array([0. for k in range(count)])
    alpha = a.copy()
    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -alpha[i][j] / a[i][i]
    
    print("Iterative method:")
    it = 0
    while True:
        x_prev = copy.deepcopy(x)
        for i in range(count):
            s = 0
            for j in range(count):
                if j != i:
                    s = s + a[i][j] * x_prev[j] 
            x[i] = b[i] / a[i][i] - s / a[i][i]
        print(f'{it + 1}:', x)
        if finish_iterations(alpha, x_prev, x, eps):
            break
        it += 1
    return x    
            

def seidel(a, b, eps=0.001):
    n = len(a)
    count = len(b)
    x = np.array([0. for k in range(count)])
    alpha = a.copy()
    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -alpha[i][j] / a[i][i]
    B = np.tril(alpha, -1)
    C = alpha - B
    alpha = np.dot(np.linalg.inv(np.eye(n) - B), C)

    print("Gauss-Seidel method:")
    it = 0
    while True:
        x_prev = copy.deepcopy(x)
        for i in range(count):
            s = 0
            for j in range(count):
                if j < i:
                    s = s + a[i][j] * x[j] 
                elif j > i:
                    s = s + a[i][j] * x_prev[j] 
            x[i] = b[i] / a[i][i] - s / a[i][i]
        print(f'{it + 1}:', x)
        if finish_iterations(alpha, x_prev, x, eps):
            break
        it += 1
    return x    


def main():
    # input
    np.set_printoptions(precision=8)
    n = int(input())
    eps = float(input())
    a = np.array([list(map(float, input().split())) for _ in range(n)])
    b = np.array([float(input()) for _ in range(n)])
    print("Matrix A:\n", tabulate(a), '\n')
    print("Matrix B:\n", b, '\n')

    x = iterative(a, b, eps)
    print("x =", x)
    print()
    x = seidel(a, b, eps) 
    print("x =", x)

    print("\nNumpy:")
    print("x =", np.linalg.solve(a, b))


if __name__ == "__main__":
    main()