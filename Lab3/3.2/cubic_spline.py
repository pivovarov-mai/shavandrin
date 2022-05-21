import numpy as np
import matplotlib.pyplot as plt

def calc_h(X):
    h = []
    for i in range(1, len(X)):
        h.append(X[i] - X[i - 1])
    return h


def calc_c(X, Y):
    all_h = calc_h(X)
    n = len(X)
    C = [[2 * (all_h[0] + all_h[1]), all_h[1], 0]]
    H = [[3 * ((Y[2] - Y[1]) / all_h[1] - (Y[1] - Y[0]) / all_h[0])]]
    for i in range(3, n - 1):
        C.append([all_h[i - 2], 2 * (all_h[i - 2] + all_h[i - 1]), all_h[i - 1]])
        H.append([3 * ((Y[i] - Y[i - 1]) / all_h[i - 1] - (Y[i - 1] - Y[i - 2]) / all_h[i - 2])])
    C.append([0, all_h[n - 3], 2 * (all_h[n - 3] + all_h[n - 2])])
    H.append([3 * ((Y[n - 1] - Y[n - 2]) / all_h[n - 2] - (Y[n - 2] - Y[n - 3]) / all_h[n - 3])])
    C = np.array(C)
    H = np.array(H)
    return np.linalg.solve(C, H)


def find_range(X, x_point):
    ind = 0
    for i in range(len(X) - 1):
        if X[i] <= x_point < X[i + 1]:
            ind = i
            break
    if ind == 0 or ind == len(X) - 1:
        print("This is boundary point")
        return
    return ind


def calc_spline(A, B, C, D, X, x_point, i):
    return A[i - 1] + B[i - 1] * (x_point - X[i - 1]) + C[i - 1] * (x_point - X[i - 1]) ** 2 + D[i - 1] * (x_point - X[i - 1]) ** 3


def calc_koeffs(X, Y):
    try:
        all_c = [0]
        c_from_equations = calc_c(X, Y)
        tmp = np.ndarray.flatten(c_from_equations)
        all_c.extend(tmp)
        all_h = calc_h(X)
        a = []
        b = []
        d = []
        n = len(X)
        for i in range(1, n):
            a.append(Y[i - 1])
        for i in range(1, n - 1):
            b.append((Y[i] - Y[i - 1]) / all_h[i - 1] - all_h[i - 1] * (all_c[i] + 2 * all_c[i - 1]) / 3)
            d.append((all_c[i] - all_c[i - 1]) / (3 * all_h[i - 1]))
        b.append((Y[n - 1] - Y[n - 2]) / all_h[n - 2] - 2 * all_h[n - 2] * all_c[n - 2] / 3)
        d.append(-all_c[n - 2] / (3 * all_h[n - 2]))
        return a, b, all_c, d
    except:
        print("Error")

def cubic_spline(X, Y, x_point):
    try:
        a, b, c, d = calc_koeffs(X, Y)
        index = find_range(X, x_point)
        return calc_spline(a, b, c, d, X, x_point, index)
    except:
        print("Error")

if __name__ == "__main__":
    x_point = -0.5
    X = [-3, -1, 1, 3, 5]
    Y = [-0.18016, 1.3562, 1.7814, 3.3218, 5.1974]
    spline_value = cubic_spline(X, Y, x_point)
    print(f"Spline value in point {x_point}: {spline_value}")
    # graphic of cubic spline for check
    a, b, c, d = calc_koeffs(X, Y)
    for i in range(len(X) - 1):
        x_spline = [np.around(j, 2) for j in np.arange(X[i], X[i + 1] + 0.2, 0.2)]
        y_spline = [calc_spline(a, b, c, d, X, elem, i + 1) for elem in x_spline]
        plt.plot(x_spline, y_spline, color='blue')
        plt.scatter(X[i], Y[i], color='red')
    plt.scatter(X[len(X) - 1], Y[len(X) - 1], color='red')
    plt.show()