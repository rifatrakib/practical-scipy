# Special functions ([`scipy.special`](https://docs.scipy.org/doc/scipy/reference/special.html#module-scipy.special))

The main feature of the [`scipy.special`](https://docs.scipy.org/doc/scipy/reference/special.html#module-scipy.special) package is the definition of numerous special functions of `mathematical physics`. Available functions include `airy, elliptic, bessel, gamma, beta, hypergeometric, parabolic cylinder, mathieu, spheroidal wave, struve, and kelvin`. There are also some `low-level stats functions` that are not intended for general use as an easier interface to these functions is provided by the [stats](https://docs.scipy.org/doc/scipy/tutorial/stats.html) module. Most of these functions can take array arguments and return array results following the same broadcasting rules as other math functions in `Numerical Python`. Many of these functions also accept complex numbers as input. For a complete list of the available functions with a one-line description type `>>> help(special)`. Each function also has its own documentation accessible using help.

## Bessel functions of real order(`jv`, `jn_zeros`)

Bessel functions are a family of solutions to `Bessel's differential equation` with real or complex order `alpha`:

$x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + (x^2 - \alpha^2)y = 0$

Among other uses, these functions arise in wave propagation problems, such as the vibrational modes of a thin drum head.
