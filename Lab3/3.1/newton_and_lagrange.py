import numpy as np

def f(x):
    return np.pi / 2 - np.arctan(x) + x


def w(x, X, index):
    mult = 1
    for i in range(len(X)):
        if i != index:
            mult *= x - X[i]
    return mult


def calc_lagrange(A, X, x_point):
    value = 0
    for i in range(len(A)):
        value += A[i] * w(x_point, X, i)
    return value


def lagrange_interpolate(X, x_point):
    try:
        Y = [f(x) for x in X]
        f_div_w_array = []
        for i in range(len(X)):
            w_i = w(X[i], X, i)
            f_div_w_array.append(Y[i] / w_i)
        return calc_lagrange(f_div_w_array, X, x_point)
    except:
        print("Error")


def divided_difference(X):
    if len(X) == 1:
        return f(X[0])
    return (divided_difference(X[:-1]) - divided_difference(X[1:])) / (X[0] - X[-1])


def calc_newton(all_diffs, X, x_point):
    newton = f(X[0])
    multiplier = 1
    for i in range(len(all_diffs)):
        multiplier *= (x_point - X[i])
        newton += multiplier * all_diffs[i]
    return newton


def newton_interpolate(X, x_point):
    try:
        Y = [f(X) for x in X]
        j = 0
        div_diffs = []
        for j in range(1, len(X)):
            #div_diffs.append(divided_difference(X[i:(i + j + 1)]) for i in range(len(X) - j))
            div_diffs.append(divided_difference(X[:(j + 1)]))
        return calc_newton(div_diffs, X, x_point)
    except:
        print("Error")


if __name__ == "__main__":
    x_point = -0.5
    X_a = [-3, -1, 1, 3]
    X_b = [-3, 0, 1, 3]
    l1 = lagrange_interpolate(X_a, x_point)
    l2 = lagrange_interpolate(X_b, x_point)
    true_y = f(x_point)
    print('-----------------------------------------------------------------------------------------')
    print(f'True function value: {true_y}')
    print('-----------------------------------------------------------------------------------------')
    print(f'Lagrange interpolation for first array of points: {l1}')
    print(f'Lagrange interpolation for second array of points: {l2}')
    print(f'Absolute error of Lagrange interpolation for first array of points: {abs(l1 - true_y)}')
    print(f'Absolute error of Lagrange interpolation for second array of points: {abs(l2 - true_y)}')
    print('-----------------------------------------------------------------------------------------')
    n1 = newton_interpolate(X_a, x_point)
    n2 = newton_interpolate(X_b, x_point)
    print(f'Newton interpolation for first array of points: {n1}')
    print(f'Newton interpolation for second array of points: {n2}')
    print(f'Absolute error of Newton interpolation for first array of points: {abs(n1 - true_y)}')
    print(f'Absolute error of Newton interpolation for second array of points: {abs(n2 - true_y)}')
    print('-----------------------------------------------------------------------------------------')