import math

def f(x):
    return x * math.sqrt(2 * x + 3)


def generate_X(X0, X1, h):
    X = []
    start = X0
    end = X1
    while start <= end:
        X.append(start)
        start += h
    return X


def rectangle_method(X0, X1, h):
    X = generate_X(X0, X1, h)
    integral = 0
    for i in range(0, len(X) - 1):
        integral += f((X[i] + X[i + 1]) / 2)
    return h * integral


def trapezoid_method(X0, X1, h):
    X = generate_X(X0, X1, h)
    integral = f(X[0]) / 2
    for i in range(1, len(X) - 1):
        integral += f(X[i])
    integral += f(X[len(X) - 1]) / 2
    return h * integral


def simpson_method(X0, X1, h):
    X = generate_X(X0, X1, h)
    integral = f(X[0])
    for i in range(1, len(X) - 1):
        integral += 4 * f(X[i]) if i % 2 != 0 else 2 * f(X[i])
    integral += f(X[len(X) - 1])
    return h * integral / 3


def runge_romberg_richardson(int1, int2, p, h1, h2, real_value):
    k = h1 / h2
    result = int2 + (int2 - int1) / (k ** p - 1)
    return result, real_value - result


if __name__ == "__main__":
    X0 = -1
    X1 = 1
    h1 = 0.5
    h2 = 0.25
    real_val = 0.4
    p = 2
    r1, r2 = rectangle_method(X0, X1, h1), rectangle_method(X0, X1, h2)
    t1, t2 = trapezoid_method(X0, X1, h1), trapezoid_method(X0, X1, h2)
    s1, s2 = simpson_method(X0, X1, h1), simpson_method(X0, X1, h2)
    rrr_r, err_r = runge_romberg_richardson(r1, r2, p, h1, h2, real_val)
    rrr_t, err_t = runge_romberg_richardson(t1, t2, p, h1, h2, real_val)
    rrr_s, err_s = runge_romberg_richardson(s1, s2, p, h1, h2, real_val)
    print('-------------------------------------------------------')
    print(f'Real value of integral: {real_val}')
    print('-------------------------------------------------------')
    print(f'Rectange method with step {h1}: {r1}')
    print(f'Rectange method with step {h2}: {r2}')
    print(f'Runge-Romberg-Richardson method for rectangle method: {rrr_r}')
    print(f'Error for rectangle method: {err_r}')
    print('-------------------------------------------------------')
    print(f'Trapezoid method with step {h1}: {t1}')
    print(f'Trapezoid method with step {h2}: {t2}')
    print(f'Runge-Romberg-Richardson method for trapezoid method: {rrr_t}')
    print(f'Error for trapezoid method: {err_t}')
    print('-------------------------------------------------------')
    print(f'Simpson method with step {h1}: {s1}')
    print(f'Simpson method with step {h2}: {s2}')
    print(f'Runge-Romberg-Richardson method for Simpson method: {rrr_s}')
    print(f'Error for Simpson method: {err_s}')
    print('-------------------------------------------------------')