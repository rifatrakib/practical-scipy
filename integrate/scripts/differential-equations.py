from scipy import integrate
from scipy import special
import numpy as np


def func(t, y):
    return [t * y[1], y[0]]


def gradient(t, y):
    return [[0, t], [1, 0]]


y1_0 = 1 / 3 ** (2 / 3) / special.gamma(2 / 3)
y0_0 = -1 / 3 ** (1 / 3) / special.gamma(1 / 3)
y0 = [y0_0, y1_0]
t_span = [0, 4]

sol_1 = integrate.solve_ivp(func, t_span, y0)
print(f"{sol_1.t = }")
print(f"{sol_1.t[1] = }")
print(f"{special.airy(sol_1.t)[0] = }")

rtol, atol = (1e-8, 1e-8)
sol_2 = integrate.solve_ivp(func, t_span, y0, rtol=rtol, atol=atol)
print(f"{sol_2.y[1][0::6] = }")
print(f"{special.airy(sol_2.t)[0][0::6] = }")

t = np.linspace(0, 4, 100)
sol_3 = integrate.solve_ivp(func, t_span, y0, t_eval=t)
print(f"{sol_3.y[1] = }")

sol_4 = integrate.solve_ivp(func, t_span, y0, method="Radau", jac=gradient)
print(f"{sol_4.y[1] = }")
