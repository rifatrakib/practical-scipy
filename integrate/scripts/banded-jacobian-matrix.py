from scipy import integrate
from time import time
import numpy as np


def G(u, v, f, k):
    return f * (1 - u) - u * v ** 2


def H(u, v, f, k):
    return -(f + k) * v + u * v ** 2


def grayscott1d(y, t, f, k, du, dv, dx):
    """
    Differential equations for the 1-D Gray-Scott equations.
    
    The ODEs are derived using the method of lines.
    """
    # the vectors u and v are interleaved in y
    # we define views of u and v by slicing y
    u = y[::2]
    v = y[1::2]
    
    # dydt is the return value of this function
    dydt = np.empty_like(y)
    
    # just like u and v are views of the interleaved vectors in y
    # dudt and dvdt are views of the interleaved output vectors in dydt
    dudt = dydt[::2]
    dvdt = dydt[1::2]
    
    # compute du/dt and dv/dt
    # the end points and the interior points are handled separately
    dudt[0] = G(u[0], v[0], f, k) + du * (-2.0 * u[0] + 2.0 * u[1]) / dx ** 2
    dudt[1:-1] = G(u[1:-1], v[1:-1], f, k) + du * np.diff(u, 2) / dx ** 2
    dudt[-1] = G(u[-1], v[-1], f, k) + du * (-2.0 * u[-1] + 2.0 * u[-2]) / dx ** 2
    
    dvdt[0] = H(u[0], v[0], f, k) + dv * (-2.0 * v[0] + 2.0 * v[1]) / dx ** 2
    dvdt[1:-1] = H(u[1:-1], v[1:-1], f, k) + dv * np.diff(v, 2) / dx ** 2
    dvdt[-1] = H(u[-1], v[-1], f, k) + dv * (-2.0 * v[-1] + 2.0 * v[-2]) / dx ** 2
    
    return dydt


rng = np.random.default_rng()
y0 = rng.standard_normal(5000)
t = np.linspace(0, 50, 11)
f = 0.024
k = 0.055
du = 0.01
dv = 0.005
dx = 0.025

start = time()
sol_a = integrate.odeint(grayscott1d, y0, t, args=(f, k, du, dv, dx))
print(f"{sol_a = } was computed in {time() - start} s")

start = time()
sol_b = integrate.odeint(grayscott1d, y0, t, args=(f, k, du, dv, dx), ml=2, mu=2)
print(f"{sol_b = } was computed in {time() - start} s")

print(f"{np.allclose(sol_a, sol_b) = }")
