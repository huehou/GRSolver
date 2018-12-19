# GR Solver

## Disclaimer
This program translates the codes in **Curvature and the Einstein Equation** from *Mathematica* to *Python*. As such, most text in the original file has been copied into the README file. Those who are interested can download the original *Mathematica* file from [here](http://web.physics.ucsb.edu/~gravitybook/mathematica.html).

## Introduction

This is a Python program that translates the *Mathematica* notebook **Curvature and the Einstein Equation** to a Python program. As of now, the plan is, from a given metric $g_{\alpha \beta}$, it computes the components of the following:
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
In the output the symbol `Gamma[1,2,3]` stands for $\Gamma^{1}_{23}$. Because the Christoffel symbols are symmetric under interchange of the last two indices, only the independent components are displayed.

## Calculating Riemann tensor

The components of the Riemann tensor, $R^{\lambda}_{\mu \nu \sigma}$ are calculated using the definition

> $$R^{\lambda}_{\mu \nu \sigma} = \partial_\nu \Gamma^{\lambda}_{\mu \sigma} - \partial_\sigma \Gamma^{\lambda}_{\mu \nu} + \Gamma^{\eta}_{\mu \sigma} \Gamma^{\lambda}_{\eta \nu} - \Gamma^{\eta}_{\mu \nu} \Gamma^{\lambda}_{\eta \sigma} \;.$$

To calculate the Riemann tensor, one needs only to call

```python
test.Riem()
```

If one is interested in only the non-zero components of the Riemann tensor, then one can call

```python
test.Riem(True)
```
In the output, the symbol `R[1,2,1,3]` stands for $R^{1}_{213}$, and similarly for the other components. You can obtain `R[1,2,3,1]` from `R[1,2,1,3]` using the antisymmetry of the Riemann tensor under exchange of the last two indices. The antisymmetry under exchange of the first two indices of $R_{\lambda \mu \nu \sigma}$ is not evident in the output because the components of $R^{\lambda}_{\mu \nu \sigma}$ are displayed.

## Calculating the Ricci tensor

The Ricci tensor $R_{\mu \nu}$ is obtained by contracting the first and third indices of the Riemann tensor,

> $$R_{\mu \nu} = R^{\lambda}_{\mu \lambda \nu}$$

This can be calculated by calling
```python
test.Ric()
```
To display the non-zero components of the Ricci tensor, use
```python
test.Ric(True)
```
Since Ricci tensor is symmetric, only the independent components are displayed. A vanishing table (as with the Schwarzschild metric example) means that the vacuum Einstein equation is satisfied.s

## Calculating Ricci Scalar

The Ricci scalar $R$ is calculated by contracting the Ricci tensor. This is done with the help of the contravariant metric

> $$R = g^{\mu \nu} R_{\mu \nu} \;.$$

This can be calculated by calling

```python
print(test.RicSca())
```

## Calculating the Einstein tensor

The Einstein tensor is

> $$ G_{\mu \nu} = R_{\mu \nu} - \frac{1}{2} g_{\mu \nu} R \;,

which is the left hand side of the Einstein field equation

> $$G_{\mu \nu} = 8 \pi T_{\mu \nu}$$

This can be calculated by calling
```python
test.Eins()
```
To display only non-zero components of the Einstein tensor, one can call
```python
test.Eins(True)
```
to print the non-zero components of the Einstein tensor. A vanishing table meanst that the vacuum Einstein equatoin is satisfied.