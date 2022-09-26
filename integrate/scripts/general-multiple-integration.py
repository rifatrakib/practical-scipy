from scipy.integrate import nquad, dblquad
import numpy as np


def I(n):
    return dblquad(lambda t, x: np.exp(-x * t) / t ** n, 0, np.inf, lambda x: 1, lambda x: np.inf)


print(f"{I(4) = }")
print(f"{I(3) = }")
print(f"{I(2) = }")

area = dblquad(lambda x, y: x * y, 0, 0.5, lambda x: 0, lambda x: 1 - 2 * x)
print(f"{area = }")

n = 5


def f(t, x):
    return np.exp(-x * t) / t ** n


print(f"{nquad(f, [[1, np.inf], [0, np.inf]]) = }")


def f_n(x, y):
    return x * y


def bounds_y():
    return [0, 0.5]


def bounds_x(y):
    return [0, 1 - 2 * y]


print(f"{nquad(f_n, [bounds_x, bounds_y]) = }")
