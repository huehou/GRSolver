# GR Solver

This is a Python program that translates the *Mathematica* notebook **Curvature and the Einstein Equation** to a Python program. As of now, the plan is, from a given metric !equation(https://latex.codecogs.com/gif.latex?g_%7B%5Cmu%20%5Cnu%7D), it computes the components of the following:
- The inverse metric, $g^{\alpha \beta}$
- The Christoffel symbols, $\Gamma^{\lambda}_{\mu \nu}$
- The Riemann tensor, $R^{\lambda}_{\mu \nu \sigma}$
- The Ricci tensor, $R_{\mu \nu}$
- The scalar curvature, $R$
- The Einstein tensor, $G_{\mu \nu}$

You must input the covariant components of the metric tensor providing the relevant input. You may also wish to change the names of the coordinates. All the components computed are in the *coordinate basis* in which the metric was specified.

## Defining a list of coordinates
As an example, let's consider the Schwarzschild metric. We first setup the coordinate system by specifying the name. This is done using the `sympy` module.
```python
import sympy as sym
t, r, theta, phi, m = sym.symbols('t r theta phi m')
coord = [t, r, theta, phi]
```
You can change the names of the coordinates by choosing different inputs for `coord`, for example, to `[t, x, y, z]`. when another set of coordinate names is more appropriate. In this program indices range over 0 to $n$, where $n$ is the number of dimensions.

## Defining the metric
Input the metric as a list of lists, i.e., as a matrix. You can input the components of any metric here, but you must specify them as explicit functions of the coordinates.
```python
metric = [[-1 + 2*m/r, 0, 0, 0], [0, (1 - 2*m/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sym.sin(theta)**2]]
```
Once the metric and coordinates are defined, we are ready to use the GR solver. Using the solver is as simple as 
```python
test = GR(coord, metric)
```
and we are ready to go.

## Calculating the inverse metric
The contravariant version of the metric $g^{\alpha \beta}$ is the inverse of the metric $g_{\alpha \beta}$. This is obtained through matrix inversion, which can be called by
```python
print(np.array(test.inversemetric))
```

## Calculating the Christoffel symbols
The Christoffel symbols are defined as
> $$\Gamma^{\lambda}_{\mu \nu} = \frac{1}{2} g^{\lambda \sigma} (\partial_\mu g_{\sigma \nu} + \partial_\nu g_{\sigma \mu} - \partial_\sigma g_{\mu \nu}) \;.$$

To calculate the Christoffel symbols, one needs only to call
```python
print(np.array(test.Christ()))
```
If one is interested in only the non-zero Christoffel symbols, then one can call
```python
test.Christ(True)
```
In the output the symbol `Gamma[1,2,3]` stands for $\Gamma^1_{23}$. Because the Christoffel symbols are symmetric under interchange of the last two indices, only the independent components are displayed.