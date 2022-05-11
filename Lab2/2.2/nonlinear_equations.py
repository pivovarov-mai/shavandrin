import numpy as np

def f1(X):
    return X[0] ** 2 - X[0] + X[1] ** 2 - 1


def f2(X):
    return X[1] - np.tan(X[0])


def f1_dx1(x1):
    return 2 * x1 - 1


def f1_dx2(x2):
    return 2 * x2


def f2_dx1(x1):
    return -1 / np.cos(x1) ** 2


def f2_dx2(x1):
    return 1


def phi1(X):
    return np.sqrt(X[0] - X[1] ** 2 + 1)


def phi2(X):
    return np.tan(X[0])


def phi1_dx1(X):
    return 1 / (2 * np.sqrt(X[0] - X[1] ** 2 + 1))


def phi1_dx2(X):
    return -X[1] / np.sqrt(X[0] - X[1] ** 2 + 1)


def phi2_dx1(X):
    return 1 / np.cos(X[0]) ** 2


def phi2_dx2(X):
    return 0


def jacobi(X):
    jacobi = []
    jacobi.append([f1_dx1(X[0]), f1_dx2(X[1])])
    jacobi.append([f2_dx1(X[0]), f2_dx2(X[1])])
    return np.array(jacobi)


def calc_norm(v1, v2):
    return max(v2[0] - v1[0], v2[1] - v1[1])


def calc_q1(a1, b1, a2, b2):
    mid1 = (a1 + b1) / 2
    mid2 = (a2 + b2) / 2
    val1 = mid1 + abs(b1 - a1)
    val2 = mid2 + abs(b2 - a2)
    max1 = abs(phi1_dx1([val1, val2])) + abs(phi1_dx2([val1, val2]))
    max2 = abs(phi2_dx1([val1, val2])) + abs(phi2_dx2([val1, val2]))
    return max(max1, max2)


def calc_q(a1, b1, a2, b2):
    q = None
    for x1 in [a1, b1]:
        for x2 in [a2, b2]:
            max1 = abs(phi1_dx1([x1, x2])) + abs(phi1_dx2([x1, x2]))
            max2 = abs(phi2_dx1([x1, x2])) + abs(phi2_dx2([x1, x2]))
            q_cur = max(max1, max2)
        if q is None or q_cur > q:
            q = q_cur
    return q


def newton(a1, b1, a2, b2, eps):
    try:
        x_k = np.array([(a1 + b1) / 2, (a2 + b2) / 2])
        inverse_jacobi = np.linalg.inv(jacobi(x_k))
        x_k_next = x_k - np.dot(inverse_jacobi, np.array([f1(x_k), f2(x_k)]))
        iterations = 1
        while calc_norm(x_k, x_k_next) > eps:
            x_k = x_k_next
            x_k_next = x_k - np.dot(inverse_jacobi, np.array([f1(x_k), f2(x_k)]))
            iterations += 1
        return x_k_next, iterations
    except:
        print("Error")
    return


def simple_iteration(a1, b1, a2, b2, eps):
    try:
        x_k = np.array([(a1 + b1) / 2, (a2 + b2) / 2])
        x_k_next = np.array([phi1(x_k), phi2(x_k)])
        iterations = 1
        q = calc_q(a1, b1, a2, b2)
        while q / (1 - q) * calc_norm(x_k, x_k_next) > eps:
            x_k = x_k_next
            x_k_next = np.array([phi1(x_k), phi2(x_k)])
            iterations += 1
        return x_k_next, iterations
    except:
        print("Error")


if __name__ == "__main__":
    # It was defined graphically that x1 is in [0.7, 1] range, x2 is in [1, 1.2] range.
    eps = float(input("Enter epsilon: "))
    a1 = float(input("Enter a1: "))
    b1 = float(input("Enter b1: "))
    a2 = float(input("Enter a2: "))
    b2 = float(input("Enter b2: "))
    print('-----------------------------------------------------------------------------------------')
    print("Newton method")
    solution = newton(a1, b1, a2, b2, eps)
    print(f"The root of equation: {solution[0]}")
    print(f"The number of iterations: {solution[1]}")
    print('-----------------------------------------------------------------------------------------')
    print("Simple iteration method")
    solution = simple_iteration(a1, b1, a2, b2, eps)
    print(f"The root of equation: {solution[0]}")
    print(f"The number of iterations: {solution[1]}")
    print('-----------------------------------------------------------------------------------------')