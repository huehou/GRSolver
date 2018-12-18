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
        self.Ricci = None

    def Christ(self, display = False):
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
        if display:
            Christoffel = self.Christoffel
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(j+1):
                        if Christoffel[i,j,k] != 0:
                            print('Gamma[{},{},{}] = {}'.format(i,j,k,Christoffel[i,j,k]))
        return self.Christoffel

    def Riem(self, display = False):
        if type(self.Riemann) == type(None):
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
                            Riemann[i,j,k,l] = sym.simplify(res)
            self.Riemann = Riemann
        if display:
            Riemann = self.Riemann
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(self.dim):
                        for l in range(k):
                            if Riemann[i,j,k,l] != 0:
                                print('Riemann[{},{},{},{}] = {}'.format(i,j,k,l,Riemann[i,j,k,l]))
        return self.Riemann

    def Ric(self, display = False):
        if type(self.Ricci) == type(None):
            Ricci = np.ones([self.dim, self.dim])*coord[0]
            for j in range(self.dim):
                for l in range(self.dim):
                    res = 0
                    for i in range(self.dim):
                        res += self.Riem()[i,j,i,l]
                    Ricci[j,l] = sym.simplify(res)
            self.Ricci = Ricci
        if display:
            Ricci = self.Ricci
            for j in range(self.dim):
                for l in range(j+1):
                    if Ricci[j,l] != 0:
                        print('Ricci[{},{}] = {}'.format(j, l, Ricci[j,l]))
        return self.Ricci

if __name__ == "__main__":
    t, r, theta, phi = sym.symbols('t r theta phi')
    coord = [t, r, theta, phi]
    m = sym.symbols('m')
    metric = sym.Matrix([[-1 + 2*m/r, 0, 0, 0], [0, (1 - 2*m/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sym.sin(theta)**2]])
    test = GR(coord, metric)
    test.Ric(True)