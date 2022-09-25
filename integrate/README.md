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
