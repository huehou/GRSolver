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
        self.Riemann = None

    def Christ(self, prnt = False):
        '''
        Calculate the Christoffel symbols of a given metric.
        Only independent non-zero components are displayed. The Christoffel symbols are symmetric under interchange of the last two indices.
        '''
        if type(self.Christoffel) == type(None):
            Christoffel = np.ones([self.dim, self.dim, self.dim])*coord[0]
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(self.dim):
                        res = 0
                        for s in range(self.dim):
                            res += self.inversemetric[i,s]*(sym.diff(self.metric[s,j], self.coord[k]) + sym.diff(self.metric[s,k], self.coord[j]) - sym.diff(self.metric[j,k], self.coord[s]))
                        res /= 2
                        Christoffel[i,j,k] = sym.simplify(res)
            self.Christoffel = Christoffel
        Christoffel = self.Christoffel
        if prnt:
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(j+1):
                        if Christoffel[i,j,k] != 0:
                            print('Gamma[{},{},{}] = {}'.format(i,j,k,Christoffel[i,j,k]))
        return self.Christoffel

    def Riem(self):
        Riemann = np.ones([self.dim, self.dim, self.dim, self.dim])*coord[0]
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        res = 0
                        for s in range(self.dim):
                            res += self.Christ()[s,j,l] * self.Christ()[i,k,s] - self.Christ()[s,j,k] * self.Christ()[i,l,s]
                        res += sym.diff(self.Christ()[i,j,l], coord[k])
                        res -= sym.diff(self.Christ()[i,j,k], coord[l])
                        Riemann[i,j,k,l] = res
        self.Riemann = Riemann
        return Riemann

if __name__ == "__main__":
    t, r, theta, phi = sym.symbols('t r theta phi')
    coord = [t, r, theta, phi]
    m = sym.symbols('m')
    metric = sym.Matrix([[-1 + 2*m/r, 0, 0, 0], [0, (1 - 2*m/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sym.sin(theta)**2]])
    test = GR(coord, metric)
    test.Riem()