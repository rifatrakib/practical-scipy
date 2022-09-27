## Integration (`scipy.integrate`)

The `scipy.integrate` sub-package provides _several integration techniques_ including an _ordinary differential equation integrator_. An overview of the module is provided by the help command.


#### General integration (`quad`)

The function `quad` is provided to integrate a function of one variable between two points. The points can be $\pm\infty$ ($\pm inf$) to indicate infinite limits. For example, suppose you wish to integrate a bessel function `jv(2.5, x)` along the interval `[0, 4.5]`.

$I=\int_{0}^{4.5}J_{2.5}\left(x\right)\, dx.$

This could be computed using `quad`.

The first argument to `quad` is a `"callable" Python object` (_i.e., a function, method, or class instance_). Notice the use of a `lambda function` in this case as the argument. The next _two arguments_ are the __limits of integration__. The return value is a `tuple`, with the `first element` holding the _estimated value of the integral_ and the `second element` holding an _upper bound on the error_. Notice, that in this case, the true value of this integral is

$I=\sqrt{\frac{2}{\pi}}\left(\frac{18}{27}\sqrt{2}\cos\left(4.5\right)-\frac{4}{27}\sqrt{2}\sin\left(4.5\right)+\sqrt{2\pi}\textrm{Si}\left(\frac{3}{\sqrt{\pi}}\right)\right)$,

where

$\textrm{Si}\left(x\right)=\int_{0}^{x}\sin\left(\frac{\pi}{2}t^{2}\right)\, dt$.

is the `Fresnel sine integral`. Note that the `numerically-computed integral` is within $1.04\times10^{-11}$  of the exact result — well below the reported error bound.

If the function to integrate takes additional parameters, they can be provided in the args argument. Suppose that the following integral shall be calculated:

$I(a,b)=\int_{0}^{1} ax^2+b \, dx$.

`Infinite` inputs are also allowed in `quad` by using $\pm inf$ as one of the arguments. For example, suppose that a numerical value for the exponential integral:

$E_{n}\left(x\right)=\int_{1}^{\infty}\frac{e^{-xt}}{t^{n}}\, dt$.

is desired (and the fact that _this integral can be computed_ as `special.expn(n,x)` is forgotten). The functionality of the function `special.expn` can be _replicated_ by defining a new function `vec_expint` based on the routine `quad`.

The function which is integrated can even use the `quad` argument (though the `error bound` may _underestimate_ the error due to possible `numerical error` in the _integrand_ from the use of `quad`). The integral in this case is

$I_{n}=\int_{0}^{\infty}\int_{1}^{\infty}\frac{e^{-xt}}{t^{n}}\, dt\, dx=\frac{1}{n}$.


## General multiple integration (`dblquad`, `tplquad`, `nquad`)

The mechanics for _double and triple integration_ have been wrapped up into the functions `dblquad` and `tplquad`. These functions take the `function to integrate` and `four, or six arguments`, respectively. _The limits of all inner integrals need to be defined as functions_.

As example for `non-constant limits` consider the integral

$I=\int_{y=0}^{1/2}\int_{x=0}^{1-2y} x y \, dx\, dy=\frac{1}{96}$.

This `integral` can be evaluated using the expression below (Note the use of the `non-constant lambda functions` for the `upper limit` of the `inner integral`).

For `n-fold integration`, `scipy` provides the function `nquad`. The `integration bounds` are an `iterable` object: either _a list of constant bounds_, or _a list of functions for the non-constant integration bounds_. The `order of integration` (and therefore the `bounds`) is from the _innermost integral_ to the _outermost_ one.

Note that the _order of arguments for `f` **must match** the order of the integration bounds; `i`_.e., the `inner integral` with respect to $t$ is on the interval $[1, \infty]$ and the `outer integral` with respect to $x$ is on the interval $[0, \infty]$. `Non-constant integration bounds` can be treated in a similar manner.


## Integrating using Samples

If the samples are _equally-spaced_ and the _number of samples_ available is $2^{k}+1$ for some integer $k$, then __Romberg `romb` integration__ can be used to obtain _high-precision estimates of the integral_ using the available samples. _`Romberg integration`_ uses the `trapezoid rule` at `step-sizes` related by a `power of two` and then performs _`Richardson extrapolation`_ on these estimates to approximate the integral with a higher degree of accuracy.

In case of _arbitrary spaced samples_, the two functions `trapezoid` and `simpson` are available. They are using _`Newton-Coates formulas` of order 1 and 2 respectively_ to perform integration. The `trapezoidal rule` approximates the function as a _straight line between adjacent points_, while `Simpson`'s rule approximates the function between _three adjacent points as a parabola_.

For an _odd number of samples_ that are `equally spaced Simpson's rule` is __exact__ if the function is a _polynomial of order 3 or less_. If the samples are _not equally spaced_, then the result is _**exact** only if the function is a `polynomial of order 2 or less`_.


## Ordinary differential equations (solve_ivp)

Integrating a set of `ordinary differential equations (ODEs)` given _initial conditions_ is another useful example. The function `solve_ivp` is available in `SciPy` for integrating a `first-order vector differential equation`:

$\frac{d\mathbf{y}}{dt}=\mathbf{f}\left(\mathbf{y},t\right)$,

given initial conditions $\mathbf{y}\left(0\right)=y_{0}$, where $\mathbf{y}$ is a length $N$ vector and $\mathbf{f}$ is a mapping from $\mathcal{R}^{N}$ to $\mathcal{R}^{N}$. A `higher-order ordinary differential equation` can always be __reduced__ to a _differential equation_ of this type by introducing __`intermediate derivatives`__ into the $\mathbf{y}$ vector.

For example, suppose it is desired to find the solution to the following second-order differential equation:

$\frac{d^{2}w}{dz^{2}}-zw(z)=0$

with _initial conditions_ $w\left(0\right)=\frac{1}{\sqrt[3]{3^{2}}\Gamma\left(\frac{2}{3}\right)}$ and $\left.\frac{dw}{dz}\right|_{z=0}=-\frac{1}{\sqrt[3]{3}\Gamma\left(\frac{1}{3}\right)}$. It is known that the solution to this _differential equation_ with these `boundary conditions` is the __`Airy function`__

$w=\textrm{Ai}\left(z\right)$,

which gives a means to check the _integrator_ using `special.airy`.

First, convert this `ODE` into _standard form_ by setting $\mathbf{y}=\left[\frac{dw}{dz},w\right]$ and $t=z$. Thus, the _differential equation_ becomes

$\frac{d\mathbf{y}}{dt}=\left[\begin{array}{c} ty_{1}\\ y_{0}\end{array}\right]=\left[\begin{array}{cc} 0 & t\\ 1 & 0\end{array}\right]\left[\begin{array}{c} y_{0}\\ y_{1}\end{array}\right]=\left[\begin{array}{cc} 0 & t\\ 1 & 0\end{array}\right]\mathbf{y}$.

In other words,

$\mathbf{f}\left(\mathbf{y},t\right)=\mathbf{A}\left(t\right)\mathbf{y}$.

As an interesting reminder, if $\mathbf{A}\left(t\right)$ commutes with $\int_{0}^{t}\mathbf{A}\left(\tau\right)\, d\tau$ under __`matrix multiplication`__, then this `linear differential equation` has an exact solution using the __matrix exponential__:

$\mathbf{y}\left(t\right)=\exp\left(\int_{0}^{t}\mathbf{A}\left(\tau\right)d\tau\right)\mathbf{y}\left(0\right)$,

However, in this case, $\mathbf{A}\left(t\right)$ and its integral _do not commute_.

This `differential equation` can be solved using the function `solve_ivp`. It requires the _`derivative`, `fprime`, the `time span [t_start, t_end]` and the initial conditions vector, `y0`_, as input arguments and returns an object whose `y` field is _an array_ with consecutive solution values as columns. The initial conditions are therefore given in the first output column.

As it can be seen `solve_ivp` _determines its time steps_ **automatically** if not specified otherwise. To compare the solution of `solve_ivp` with the `airy` function the _time vector_ created by `solve_ivp` is passed to the `airy` function.

The solution of `solve_ivp` with its _standard parameters_ shows a _big deviation_ to the `airy` function. To __minimize__ this deviation, _relative and absolute tolerances_ can be used.

To specify _user defined time points_ for the solution of `solve_ivp`, `solve_ivp` offers __two possibilities__ that can also be used _complementarily_. By passing the `t_eval` option to the function call `solve_ivp` returns the solutions of these time points of `t_eval` in its output.

If the `jacobian matrix` of function is known, it can be passed to the `solve_ivp` to _achieve better results_. Please be aware however that the _default integration method_ __`RK45`__ _**does not support** jacobian matrices_ and thereby another integration method has to be chosen. One of the integration methods that support a jacobian matrix is the `Radau` method.


## Solving a system with a `banded Jacobian matrix`

`odeint` can be told that the `Jacobian` is __banded__. For a _large system of differential equations_ that are known to be __stiff__, this can `improve performance significantly`.

As an example, we'll solve the `1-D Gray-Scott partial differential equations` using the _method of lines [`MOL`]_. The `Gray-Scott equations` for the functions $u(x, t)$ and $v(x, t)$ on the interval $x \in [0, L]$ are

$\frac{\partial u}{\partial t} = D_u \frac{\partial^2 u}{\partial x^2} - uv^2 + f(1-u) \\$
$\frac{\partial v}{\partial t} = D_v \frac{\partial^2 v}{\partial x^2} + uv^2 - (f + k)v \\$

where $D_u$ and $D_v$ are the `diffusion coefficients` of the components $u$ and $v$, respectively, and $f$ and $k$ are `constants`. (For more information about the system, see [http://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/](http://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/))

We'll assume `Neumann boundary conditions` (i.e., `"no flux"`):

$\frac{\partial u}{\partial x}(0,t) = 0, \quad$
$\frac{\partial v}{\partial x}(0,t) = 0, \quad$
$\frac{\partial u}{\partial x}(L,t) = 0, \quad$
$\frac{\partial v}{\partial x}(L,t) = 0$

To apply the `method of lines`, we __discretize__ the $x$ variable by defining the `uniformly spaced grid` of $N$ points $\left\{x_0, x_1, \ldots, x_{N-1}\right\}$, with $x_0 = 0$ and $x_{N-1} = L$. We define $u_j(t) \equiv u(x_k, t)$ and $v_j(t) \equiv v(x_k, t)$, and replace the $x$ derivatives with _finite differences_. That is,

$\frac{\partial^2 u}{\partial x^2}(x_j, t) \rightarrow
    \frac{u_{j-1}(t) - 2 u_{j}(t) + u_{j+1}(t)}{(\Delta x)^2}$

We then have a `system` of $2N$ `ordinary differential equations`:

$\frac{du_j}{dt} = \frac{D_u}{(\Delta x)^2} \left(u_{j-1} - 2 u_{j} + u_{j+1}\right)
          -u_jv_j^2 + f(1 - u_j) \\
    \frac{dv_j}{dt} = \frac{D_v}{(\Delta x)^2} \left(v_{j-1} - 2 v_{j} + v_{j+1}\right)
          + u_jv_j^2 - (f + k)v_j$

For convenience, the $(t)$ arguments have been dropped.

To __enforce__ the `boundary conditions`, we introduce `"ghost" points` $x_{-1}$ and $x_N$, and define $u_{-1}(t) \equiv u_1(t)$, $u_N(t) \equiv u_{N-2}(t)$; $v_{-1}(t)$ and $v_N(t)$
are defined analogously.

Then,

$\frac{du_0}{dt} = \frac{D_u}{(\Delta x)^2} \left(2u_{1} - 2 u_{0}\right)
          -u_0v_0^2 + f(1 - u_0) \\
    \frac{dv_0}{dt} = \frac{D_v}{(\Delta x)^2} \left(2v_{1} - 2 v_{0}\right)
          + u_0v_0^2 - (f + k)v_0$

and

$\frac{du_{N-1}}{dt} = \frac{D_u}{(\Delta x)^2} \left(2u_{N-2} - 2 u_{N-1}\right)
          -u_{N-1}v_{N-1}^2 + f(1 - u_{N-1}) \\
    \frac{dv_{N-1}}{dt} = \frac{D_v}{(\Delta x)^2} \left(2v_{N-2} - 2 v_{N-1}\right)
          + u_{N-1}v_{N-1}^2 - (f + k)v_{N-1}$

We can now starting implementing this system in code. We must combine $\{u_k\}$ and $\{v_k\}$ into a single vector of length $2N$. The two obvious choices are $\{u_0, u_1, \ldots, u_{N-1}$, $v_0, v_1, \ldots, v_{N-1}\}$ and $\{u_0, v_0, u_1, v_1, \ldots, u_{N-1}$, $v_{N-1}\}$. Mathematically, it does not matter, but the choice affects how efficiently `odeint` can solve the system. The reason is in how the `order` affects the _pattern of the nonzero elements_ of the `Jacobian matrix`.

When the variables are __interleaved__, the _`bandwidth` is much smaller_. The `main diagonal` and the `two diagonals immediately above` and the `two immediately below` the main diagonal are the __nonzero diagonals__. This is important, because the inputs `mu` and `ml` of `odeint` are the _upper and lower bandwidths_ of the `Jacobian matrix`. When the variables are __interleaved__, `mu` and `ml` are 2. When the variables are stacked with $\{v_k\}$ following $\{u_k\}$, the _upper and lower bandwidths_ are N.

With that decision made, we can write the function that implements the system of differential equations.

* First, we define the functions for the `source` and `reaction` terms of the system.

* Next, we define the function that computes the right-hand side of the system of differential equations.

We won't implement a function to compute the `Jacobian`, but we will tell `odeint` that the _Jacobian matrix is banded_. This allows the `underlying solver (LSODA)` to _avoid computing values that it knows are zero_. For a _large system_, this `improves the performance significantly`.

* First, we define the required inputs.

* Time the computation without taking advantage of the banded structure of the Jacobian matrix.

* Now set `ml=2` and `mu=2`, so `odeint` knows that the `Jacobian matrix` is __banded__.
