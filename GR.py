import numpy as np
import sympy as sym

class GR():
    '''
    From a given metric, it computes the components of the following:
    - The inverse metric
    - The Christoffel symbols
    - The Riemann tensor
    - The Ricci tensor
    - The scalar curvature
    - The Einstein tensor
    '''

    def __init__(self, coord, metric):
        self.coord = coord # Defining a list of coordinates
        self.metric = metric # Defining the metric
        self.inversemetric = metric**(-1)
        self.dim = len(coord)
        self.Christoffel = None

    def Christ(self):
        pass 

if __name__ == "__main__":
    t, r, theta, phi = sym.symbols('t r theta phi')
    coord = [t, r, theta, phi]
    m = sym.symbols('m')
    metric = sym.Matrix([[-1 + 2*m/r, 0, 0, 0], [0, (1 - 2*m/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sym.sin(theta)**2]])
    test = GR(coord, metric)
    print(test.metric)
    print(test.inversemetric)
    print(test.dim)