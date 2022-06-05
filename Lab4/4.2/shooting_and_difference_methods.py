import numpy as np
import matplotlib.pyplot as plt

# epsilon
eps = 0.001

# Functions from previous labs 4.1 and 1.2

def Runge_Kutta_method(f, g, y0, z0, interval, h, return_Z=False):
    start, finish = interval
    X = [np.around(j, 5) for j in np.arange(start, finish + h, h)]
    Y = [y0]
    Z = [z0]
    for i in range(len(X) - 1):
        K1 = h * g(X[i], Y[i], Z[i])
        L1 = h * f(X[i], Y[i], Z[i])
        K2 = h * g(X[i] + h / 2, Y[i] + K1 / 2, Z[i] + L1 / 2)
        L2 = h * f(X[i] + h / 2, Y[i] + K1 / 2, Z[i] + L1 / 2)
        K3 = h * g(X[i] + h / 2, Y[i] + K2 / 2, Z[i] + L2 / 2)
        L3 = h * f(X[i] + h / 2, Y[i] + K2 / 2, Z[i] + L2 / 2)
        K4 = h * g(X[i] + h, Y[i] + K3, Z[i] + L3)
        L4 = h * f(X[i] + h, Y[i] + K3, Z[i] + L3)
        delta_y = (K1 + 2 * K2 + 2 * K3 + K4) / 6
        delta_z = (L1 + 2 * L2 + 2 * L3 + L4) / 6
        Y.append(np.around(Y[i] + delta_y, 5))
        Z.append(np.around(Z[i] + delta_z, 5))
    if not return_Z:
        return X, Y
    return X, Y, Z


def tridiagonal(a, b):
    n = len(a)
    x = np.array([0. for _ in range(n)])
    v = [0 for _ in range(n)]
    u = [0 for _ in range(n)]
    # Forward
    v[0] = a[0][1] / -a[0][0]
    u[0] = b[0] / a[0][0]
    for i in range(1, n - 1):
        v[i] = a[i][i + 1] / (-a[i][i] - a[i][i - 1] * v[i - 1])
        u[i] = (a[i][i - 1] * u[i - 1] - b[i]) / (-a[i][i] - a[i][i - 1] * v[i - 1])
    v[n - 1] = 0
    u[n - 1] = (a[n - 1][n - 2] * u[n - 2] - b[n - 1]) / (-a[n - 1][n - 1] - a[n - 1][n - 2] * v[n - 2])
    
    # Backward
    x[n - 1] = u[n - 1]
    for i in range(n - 1, 0, -1):
        x[i - 1] = v[i - 1] * x[i] + u[i - 1]
    return x  


# Source functions

def f(x, y, z):
    global eps
    return ((2 * x + 4) * z - 2 * y) /  (eps + x * (x + 4))


def g(x, y, z):
    return z


# Functions for finite difference method
# y'' + p_fd(x)y' + q_fd(x)y = f_fd(x)

def p_fd(x):
    global eps
    return -(2 * x + 4) / (x * (x + 4))


def q_fd(x):
    return 2 / (x * (x + 4))


def f_fd(x):
    return 0


# Exact solution for differential equation

def exact_solution(x):
    return x ** 2 + x + 2


# Function for calculating true values of function

def count_true_values(X):
    Y_true = []
    for i in range(len(X)):
        Y_true.append(np.around(exact_solution(X[i]), 5))
    return Y_true


# Shooting method

def shooting_method(f, g, y0, yn, interval, h, eps):
    n_prev = 1.0
    n = 0.8
    iterations = 0
    while True:
        iterations += 1
        x_prev, y_prev = Runge_Kutta_method(f, g, y0, n_prev, interval, h)
        x, y = Runge_Kutta_method(f, g, y0, n, interval, h)
        if abs(y[-1] - yn) < eps:
            break
        n_prev, n = n, n - (y[-1] - yn) * (n - n_prev) / ((y[-1] - yn) - (y_prev[-1] - yn))
    return x, y, iterations


def finite_difference_method(p, q, f, y0, yn, interval, h):
    A = []
    B = []
    rows = []
    a, b = interval
    x = np.arange(a, b + h, h)
    n = len(x)

    for i in range(n):
        if i == 0:
            rows.append(1)
        else:
            rows.append(0)
    A.append(rows)
    B.append(y0)

    for i in range(1, n - 1):
        rows = []
        B.append(f(x[i]))
        for j in range(n):
            if j == i - 1:
                rows.append(1 / h ** 2 - p(x[i]) / (2 * h))
            elif j == i:
                rows.append(-2 / h ** 2 + q(x[i]))
            elif j == i + 1:
                rows.append(1 / h ** 2 + p(x[i]) / (2 * h))
            else:
                rows.append(0)
        A.append(rows)

    rows = []
    B.append(yn)
    for i in range(n):
        if i == n - 1:
            rows.append(1)
        else:
            rows.append(0)

    A.append(rows)
    y = tridiagonal(A, B)
    x = np.ndarray.tolist(x)
    y = np.ndarray.tolist(y)
    return x, y


# Functions for errors

def runge_romberg_method(v1, v2, p, h1, h2):
    # the method only works for such steps, one of which is twice as large as the other
    # Example: h1 = 1, h2 = 0.5; (h1 = 2 * h2)
    assert h1 == h2 * 2
    error = 0
    for i in range(len(v1)):
        error += (v1[i] - v2[i * 2]) ** 2
    return error ** 0.5 / (2 ** p + 1)


def error_with_exact_solution(Y, Y_true):
    assert len(Y) == len(Y_true)
    error = 0
    for i in range(len(Y)):
        error += abs(Y_true[i] - Y[i])
    return error / len(Y)


if __name__ == '__main__':
    # initial values
    interval = (0, 2)
    y0 = 2
    y1 = 8
    h1 = 0.1
    h2 = 0.05

    # Exact solution
    x_exact = [np.around(j, 5) for j in np.arange(interval[0], interval[1] + h1, h1)]
    y_exact = [np.around(exact_solution(x), 5) for x in x_exact]
    print('-------------------------------------------------------------------------------------------------------------------')
    print(f'True x: {x_exact}')
    print(f'True y: {y_exact}')
    print('-------------------------------------------------------------------------------------------------------------------')
    plt.plot(x_exact, y_exact, label='exact solution')

    # Shooting method
    x_shooting, y_shooting, iters_shooting = shooting_method(f, g, y0, y1, interval, h1, eps)
    x_shooting2, y_shooting2, iters_shooting2 = shooting_method(f, g, y0, y1, interval, h2, eps)

    print(f'Shooting method solution: {y_shooting}')    # errors for shooting method
    error_shooting = runge_romberg_method(y_shooting, y_shooting2, 1, h1, h2)
    error_exact_shooting = error_with_exact_solution(y_shooting, y_exact)
    print('-------------------------------------------------------------------------------------------------------------------')
    print(f'Error for shooting method using Runge-Romberg method: {error_shooting}')
    print(f'Error |y_true - y_shooting|: {error_exact_shooting}')
    print('-------------------------------------------------------------------------------------------------------------------')

    # plots
    plt.plot(x_shooting, y_shooting, label=f'Shooting method, step={h1}')
    plt.plot(x_shooting2, y_shooting2, label=f'Shooting method, step={h2}')

    # Finite difference method
    x_fd, y_fd = finite_difference_method(p_fd, q_fd, f_fd, y0, y1, interval, h1)
    x_fd2, y_fd2 = finite_difference_method(p_fd, q_fd, f_fd, y0, y1, interval, h2)

    print(f'Finite difference method solution: {y_fd}')    # errors for shooting method
    error_fd= runge_romberg_method(y_fd, y_fd2, 4, h1, h2)
    error_exact_fd = error_with_exact_solution(y_fd, y_exact)
    print('-------------------------------------------------------------------------------------------------------------------')
    print(f'Error for finite difference method using Runge-Romberg method: {error_fd}')
    print(f'Error |y_true - y_fd|: {error_exact_fd}')
    print('-------------------------------------------------------------------------------------------------------------------')

    # plots
    plt.plot(x_fd, y_fd, label=f'finite difference method, step={h1}')
    plt.plot(x_fd2, y_fd2, label=f'finite difference method, step={h2}')

    plt.legend()
    plt.show()