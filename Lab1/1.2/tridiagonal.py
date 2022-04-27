import numpy as np


def solution(a, b):
    n = len(a)
    x = np.array([0. for _ in range(n)])
    v = [0 for _ in range(n)]
    u = [0 for _ in range(n)]
    # Прямой ход
    v[0] = a[0][1] / -a[0][0]
    u[0] = b[0] / a[0][0]
    for i in range(1, n - 1):
        v[i] = a[i][i + 1] / (-a[i][i] - a[i][i - 1] * v[i - 1])
        u[i] = (a[i][i - 1] * u[i - 1] - b[i]) / (-a[i][i] - a[i][i - 1] * v[i - 1])
    v[n - 1] = 0
    u[n - 1] = (a[n - 1][n - 2] * u[n - 2] - b[n - 1]) / (-a[n - 1][n - 1] - a[n - 1][n - 2] * v[n - 2])
    
    # Обратный ход
    x[n - 1] = u[n - 1]
    for i in range(n - 1, 0, -1):
        x[i - 1] = v[i - 1] * x[i] + u[i - 1]
    return x    
                

def main():
    # input
    np.set_printoptions(precision=5)
    n = int(input("Enter the matrix size: "))
    print("Enter matrix A:")
    a = np.array([list(map(float, input().split())) for _ in range(n)])
    print("Enter matrix B:")
    b = np.array([float(input()) for _ in range(n)])
    print("Matrix A:\n", a, '\n')
    print("Matrix B:", b, '\n')

    # solution
    x = solution(a, b)
    print('x =', x, '\n')

    print("A * x:", a.dot(x))


if __name__ == "__main__":
    main()