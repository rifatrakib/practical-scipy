import numpy as np
from scipy import integrate


def f1(x):
    return x ** 2


def f2(x):
    return x ** 3


x = np.array([1, 3, 4])
y1 = f1(x)
I1 = integrate.simpson(y1, x)
print(f"{I1 = }")

y2 = f2(x)
I2 = integrate.simpson(y2, x)
print(f"{I2 = }")
