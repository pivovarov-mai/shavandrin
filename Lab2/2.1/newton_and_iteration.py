import numpy as np

def f(x):
    return 10 ** x - 5 * x - 2


def derivative(x):
    return np.log(10) * 10 ** x - 5


def phi(x):
    return np.log10(5 * x + 2)


def phi_derivative(x):
    return 5 / (5 * x + 2) / np.log(10)


def newton(a, b, eps):
    try:
        iterations = 1
        x_k = (a + b) / 2
        x_k_next = x_k - f(x_k) / derivative(x_k)
        while abs(x_k_next - x_k) > eps:
            x_k = x_k_next
            x_k_next = x_k - f(x_k) / derivative(x_k)
            iterations += 1
        return x_k_next, iterations
    except:
        print("Error")


def simple_iteration(a, b, eps):
    try:
        x_k = (a + b) / 2
        x_k_next = phi(x_k)
        q = abs(phi(b))
        iterations = 1
        while q * (1 - q) * abs(x_k_next - x_k) > eps:
            x_k = x_k_next
            x_k_next = phi(x_k)
            iterations += 1
        return x_k_next, iterations
    except:
        print("Error")


if __name__ == "__main__":
    # It was defined graphically that root is in [0, 1] range.
    eps = float(input("Enter epsilon: "))
    a = float(input("Enter a: "))
    b = float(input("Enter b: "))
    print('-----------------------------------------------------------------------------------------')
    print("Newton method")
    root, iterations = newton(a, b, eps)
    print(f"The root of equation: {root}")
    print(f"The number of iterations: {iterations}")
    print('-----------------------------------------------------------------------------------------')
    print("Simple iteration method")
    root, iterations = simple_iteration(a, b, eps)
    print(f"The root of equation: {root}")
    print(f"The number of iterations: {iterations}")
    print('-----------------------------------------------------------------------------------------')