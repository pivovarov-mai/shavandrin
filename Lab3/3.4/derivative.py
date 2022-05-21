def find_x_in_points(X, val):
    for i in range(len(X)):
        if X[i] == val:
            return i
    return -1


def calc_left_derivative(x, X, Y):
    i = find_x_in_points(X, x)
    if i == -1:
        print("There is no such point!")
        return
    elif i == 0:
        print("There is no left point!")
        return 
    return (Y[i] - Y[i - 1]) / (X[i] - X[i - 1])


def calc_right_derivative(x, X, Y):
    i = find_x_in_points(X, x)
    if i == -1:
        print("There is no such point!")
        return
    elif i == len(X) - 1:
        print("There is no right point!")
        return
    return (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])


def calc_first_derivative(x, X, Y):
    left_deriv = calc_left_derivative(x, X, Y)
    right_deriv = calc_right_derivative(x, X, Y)
    if left_deriv == -1 or right_deriv == -1:
        return -1
    i = find_x_in_points(X, x)
    return left_deriv + (right_deriv - left_deriv) / (X[i + 1] - X[i - 1]) * (2 * x - X[i - 1] - X[i])


def calc_second_derivative(x, X, Y):
    left_deriv = calc_left_derivative(x, X, Y)
    right_deriv = calc_right_derivative(x, X, Y)
    if left_deriv == -1 or right_deriv == -1:
        return -1
    i = find_x_in_points(X, x)
    return 2 * (right_deriv - left_deriv) / (X[i + 1] - X[i - 1])


if __name__ == "__main__":
    X = [1, 1.2, 1.4, 1.6, 1.8]
    Y = [1, 0.69444, 0.5102, 0.39062, 0.30864]
    x_p = 1.4
    print(f'First derivative of tabular function is {calc_first_derivative(x_p, X, Y)}')
    print(f'Second derivative of tabular function is {calc_second_derivative(x_p, X, Y)}')