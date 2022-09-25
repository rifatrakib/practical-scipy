from scipy import integrate
from scipy import special
import numpy as np


result = integrate.quad(lambda x: special.jv(2.5, x), 0, 4.5)
print(f"{result = }")

I = np.sqrt(
    2 / np.pi) * (18.0 / 27 * np.sqrt(2)\
    * np.cos(4.5) - 4.0 / 27 * np.sqrt(2)\
    * np.sin(4.5) + np.sqrt(2 * np.pi)\
    * special.fresnel(3 / np.sqrt(np.pi))[0]
)
print(f"{I = }")
print(f"{abs(result[0] - I) = }")


def simple_integrand(x, a, b):
    return a * x ** 2 + b


a = 2
b = 1
I = integrate.quad(simple_integrand, 0, 1, args=(a, b))
print(f"{I = }")


def integrand(t, n, x):
    return np.exp(-x * t) / t ** n


def expint(n, x):
    return integrate.quad(integrand, 1, np.inf, args=(n, x))[0]


vec_expint = np.vectorize(expint)
result = vec_expint(3, np.arange(1.0, 4.0, 0.5))
print(f"{vec_expint(3, np.arange(1.0, 4.0, 0.5)) = }")

print(f"{special.expn(3, np.arange(1.0, 4.0, 0.5)) = }")

result = integrate.quad(lambda x: expint(3, x), 0, np.inf)
I3 = 1.0 / 3.0
print(f"{result = }")
print(f"{I3 - result[0] = }")
