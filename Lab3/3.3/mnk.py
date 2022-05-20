import numpy as np
import matplotlib.pyplot as plt

def find_coeffs(X, Y, n):
    A = []
    B = []
    for k in range(n + 1):
        A.append([sum(map(lambda X: X ** (k + i), X)) for i in range(n + 1)])
        B.append(sum(map(lambda X: X[0] * X[1] ** k, zip(Y, X))))
    A = np.array(A)
    B = np.array(B)
    return np.linalg.solve(A, B)


def mnk(x, A):
    function_value = 0
    for i in range(len(A)):
        function_value += A[i] * x ** i
    return function_value


def calc_error(calc_Y, true_Y):
    error = 0
    for i in range(len(calc_Y)):
        error += (calc_Y[i] - true_Y[i]) ** 2
    return error


if __name__ == "__main__":
    X = [-5, -3, -1, 1, 3, 5]
    Y = [-2.0558, -0.18016, 1.3562, 1.7854, 3.3218, 5.1974]
    A1 = find_coeffs(X, Y, 1)
    A2 = find_coeffs(X, Y, 2)
    approximating_polynomial1 = [mnk(x, A1) for x in X]
    approximating_polynomial2 = [mnk(x, A2) for x in X]
    error1 = calc_error(approximating_polynomial1, Y)
    error2 = calc_error(approximating_polynomial2, Y)
    # printing
    print('----------------------------------------------------------------------------------------------')
    print(f'Coefficients in approximating polynomial of the first degree: {A1}')
    print(f'Sum of quadratic errors for approximating polynomial of the first degree: {error1}')
    print('----------------------------------------------------------------------------------------------')
    print(f'Coefficients in approximating polynomial of the second degree: {A2}')
    print(f'Sum of quadratic errors for approximating polynomial of the second degree: {error2}')
    print('----------------------------------------------------------------------------------------------')
    # making plots
    A3 = find_coeffs(X, Y, 3)
    approximating_polynomial3 = [mnk(x, A3) for x in X]
    for i in range(len(X)):
        plt.plot(X[i], Y[i], 'ro')
    plt.plot(X, approximating_polynomial1, color='black', label='1')
    plt.plot(X, approximating_polynomial2, '--', color='green', label='2')
    plt.plot(X, approximating_polynomial3, '-.', color='orange', label='3')
    plt.legend()
    plt.show()