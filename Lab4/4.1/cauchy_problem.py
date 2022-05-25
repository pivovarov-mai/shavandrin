import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x):
    return np.sin(x - 1) + np.cos(x - 1) / x


def f(x, y, z):
    # y'' = z' = f(x, y, z)
    return -(x ** 2 - 2) * y / x ** 2


def g(x, y, z):
    # y' = z = g(x, y, z)
    return z


def count_true_values(X):
    Y_true = []
    for i in range(len(X)):
        Y_true.append(np.around(exact_solution(X[i]), 5))
    return Y_true


def Euler_method(interval, y0, z0, h, return_delta_Y=False):
    start, finish = interval
    X = [np.around(j, 2) for j in np.arange(start, finish + h, h)]
    Y = [y0]
    delta_Y = []
    z = z0
    Y_true = count_true_values(X)
    for i in range(len(X) - 1):
        z += h * f(X[i], Y[i], z)
        delta_Y.append(h * f(X[i], Y[i], z))
        Y.append(np.around(Y[i] + h * g(X[i], Y[i], z), 5))
    if not return_delta_Y:
        return X, Y
    return X, Y, delta_Y


def Runge_Kutta_method(interval, y0, z0, h, return_Z=False):
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


def Adams_method(interval, y0, z0, h, return_Z=False):
    x_rk, y_rk, z_rk = Runge_Kutta_method(interval, y0, z0, h, return_Z=True)
    X, Y, Z = x_rk, y_rk[:4], z_rk[:4]
    for i in range(3, len(X) - 1):
        z = Z[i] + h * (55 * f(X[i], Y[i], Z[i])
                        - 59 * f(X[i - 1], Y[i - 1], Z[i - 1])
                        + 37 * f(X[i - 2], Y[i - 2], Z[i - 2])
                        - 9 * f(X[i - 3], Y[i - 3], Z[i - 3])) / 24
        Z.append(np.around(z, 5))
        y = Y[i] + h * (55 * g(X[i], Y[i], Z[i])
                        - 59 * g(X[i - 1], Y[i - 1], Z[i - 1])
                        + 37 * g(X[i - 2], Y[i - 2], Z[i - 2])
                        - 9 * g(X[i - 3], Y[i - 3], Z[i - 3])) / 24
        Y.append(np.around(y, 5))
    if not return_Z:
        return X, Y
    return X, Y, Z


def runge_romberg_method(v1, v2, p, h1, h2, real_value):
    # the method only works for such steps, one of which is twice as large as the other
    # Example: h1 = 1, h2 = 0.5; (h1 = 2 * h2)
    assert h1 == 2 * h2
    k = h1 / h2
    result = []
    errors = []
    for i in range(len(v1)):
        result.append(np.around(v2[i * 2] + (v2[i * 2] - v1[i]) / (k ** p - 1), 5))
        errors.append(np.around(real_value[i] - result[i], 5))
    return result, sum(errors) / len(v1)


def error_with_exact_solution(Y, Y_true):
    errors = []
    assert len(Y) == len(Y_true)
    for i in range(len(Y)):
        errors.append(np.around(np.abs(Y_true[i] - Y[i]), 5))
    return sum(errors) / len(Y)


if __name__ == "__main__":
    # initializating
    interval = (1, 2)
    h = 0.1
    x0, y0, z0 = 1, 1, 0
    # printing
    # Euler method
    X_euler, Y_euler = Euler_method(interval, y0, z0, h)
    Y_true = count_true_values(X_euler)
    error_euler = error_with_exact_solution(Y_euler, Y_true)
    print('--------------------------------------------------------------------------------------------')
    print(f'True function values: {Y_true}')
    print('--------------------------------------------------------------------------------------------')
    print('Euler method:')
    print(f'X: {X_euler}')
    print(f'Y: {Y_euler}')
    print('--------------------------------------------------------------------------------------------')
    print(f'Error |y_true - y_euler|: {error_euler}')
    print('--------------------------------------------------------------------------------------------')
    # Runge_Kutta_method
    X_rk, Y_rk = Runge_Kutta_method(interval, y0, z0, h)
    error_rk = error_with_exact_solution(Y_rk, Y_true)
    print('--------------------------------------------------------------------------------------------')
    print('Runge-Kutta method:')
    print(f'X: {X_rk}')
    print(f'Y: {Y_rk}')
    print('--------------------------------------------------------------------------------------------')
    print(f'Error |y_true - y_runge_kutta|: {error_rk}')
    print('--------------------------------------------------------------------------------------------')
    # Adams method
    X_adams, Y_adams = Adams_method(interval, y0, z0, h)
    error_adams = error_with_exact_solution(Y_adams, Y_true)
    print('--------------------------------------------------------------------------------------------')
    print('Adams method:')
    print(f'X: {X_adams}')
    print(f'Y: {Y_adams}')
    print('--------------------------------------------------------------------------------------------')
    print(f'Error |y_true - y_adams|: {error_adams}')
    print('--------------------------------------------------------------------------------------------')
    # Runge-Romberg method and errors
    X_euler1, Y_euler1 = Euler_method(interval, y0, z0, 0.05, return_delta_Y=False)
    X_rk1, Y_rk1 = Runge_Kutta_method(interval, y0, z0, 0.05)
    X_adams1, Y_adams1 = Adams_method(interval, y0, z0, 0.05)
    Y_runge, error_euler = runge_romberg_method(Y_euler, Y_euler1, 1, h, 0.05, Y_true)
    Y_runge, error_kutta = runge_romberg_method(Y_rk, Y_rk1, 4, h, 0.05, Y_true)
    Y_runge, error_adams = runge_romberg_method(Y_adams, Y_adams1, 4, h, 0.05, Y_true)
    print('--------------------------------------------------------------------------------------------')
    print('Runge-Romberg method:')
    print(f'X: {X_euler1}')
    print(f'Y: {Y_runge}')
    print('--------------------------------------------------------------------------------------------')
    print(f'Error calculated by Runge-Romberg method for Euler method: {error_euler}')
    print(f'Error calculated by Runge-Romberg method for Runge-Kutta method: {error_kutta}')
    print(f'Error calculated by Runge-Romberg method for Adams method: {error_adams}')
    print('--------------------------------------------------------------------------------------------')
    # plots
    plt.plot(X_euler, Y_euler, color='red', label='Euler method, h = 0.1')
    plt.plot(X_euler1, Y_euler1, color='orange', label='Euler method, h = 0.05')
    plt.plot(X_rk, Y_rk, '--', color='green', label='Runge-Kutta method')
    plt.plot(X_adams, Y_adams, '-.', color='pink', label='Adams method')
    plt.plot(X_euler, Y_true, color='blue', label='Exact solution')
    plt.legend()
    plt.show()