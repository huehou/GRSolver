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
        '''
        Input the metric as a sympy Matrix object. You can input the components of any metric here, but you must specify them as explicit functions of the coordinates.
        '''
        self.coord = coord # Defining a list of coordinates
        self.metric = metric # Defining the metric
        self.inversemetric = metric**(-1) # The inverse metric
        self.dim = len(coord) # The dimension of the system
        self.Christoffel = None # Christoffel symbols
        self.Riemann = None # Riemann curvature tensor
        self.Ricci = None # Ricci tensor
        self.RicciScalar = None # Ricci scalar
        self.Einstein = None # Einstein tensor

    def Christ(self, display = False):
        '''
        Calculate the Christoffel symbols of a given metric.
        Only independent non-zero components are displayed. The first index is the upper index while the last 2 are the lower indices. The Christoffel symbols are symmetric under interchange of the last two indices.
        '''
        if type(self.Christoffel) == type(None):
            # Calculate only when it was not calculated before
            Christoffel = np.ones([self.dim, self.dim, self.dim])*self.coord[0] # Initialise the array with sympy symbols
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
            # To display only independent non-zero components
            Christoffel = self.Christoffel
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(j+1):
                        if Christoffel[i,j,k] != 0:
                            print('Gamma[{},{},{}] = {}'.format(i,j,k,Christoffel[i,j,k]))
        return self.Christoffel

    def Riem(self, display = False):
        '''
        Calculate the Riemann curvature tensor. Only the independent non-zero components are displayed. In the output, the first index is contravariant while the 3 other indices are covariant indices. Note that the Riemann tensor is antisymmetric under the exchange of the last two indices. The antisymmetry under exchange of the first two indices is not evident in the output because of the contravariant nature of the first index.
        '''
        if type(self.Riemann) == type(None):
            # Calculate only when it is not calculated before 
            Riemann = np.ones([self.dim, self.dim, self.dim, self.dim])*self.coord[0] # Initialise array with sympy symbols
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
            # To display only independent non-zero components
            Riemann = self.Riemann
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(self.dim):
                        for l in range(k):
                            if Riemann[i,j,k,l] != 0:
                                print('Riemann[{},{},{},{}] = {}'.format(i,j,k,l,Riemann[i,j,k,l]))
        return self.Riemann

    def Ric(self, display = False):
        '''
        Calculate the Ricci tensor by contracting the first and third indices. Only the non-zero components are displayed. The Ricci tensor is symmetric.
        '''
        if type(self.Ricci) == type(None):
            # Calculate only when it is not calculated before
            Ricci = np.ones([self.dim, self.dim])*self.coord[0] # Initialise array with sympy symbols
            for j in range(self.dim):
                for l in range(self.dim):
                    res = 0
                    for i in range(self.dim):
                        res += self.Riem()[i,j,i,l]
                    Ricci[j,l] = sym.simplify(res)
            self.Ricci = Ricci
        if display:
            # To display independent non-zero components
            Ricci = self.Ricci
            for j in range(self.dim):
                for l in range(j+1):
                    if Ricci[j,l] != 0:
                        print('Ricci[{},{}] = {}'.format(j, l, Ricci[j,l]))
        return self.Ricci

    def RicSca(self, display = False):
        '''
        Calculate the Ricci scalar by contracting the Ricci tensor. 
        '''
        if type(self.RicciScalar) == type(None):
            # Calculate only when it is not calculated before
            RicciScalar = 0
            for i in range(self.dim):
                for j in range(self.dim):
                    RicciScalar += self.inversemetric[i,j] * self.Ric()[i,j]
            self.RicciScalar = RicciScalar
        if display:
            print(self.RicciScalar)
        return self.RicciScalar

    def Eins(self, display = False):
        '''
        Calculate the Einstein tensor, which is the expression at the LHS of the Einstein field equation. A vanishing Einstein tensor means that the vacuum Einstein equation is satisfied.
        '''
        if type(self.Einstein) == type(None):
            # Calculate only when it is not calculated before
            Einstein = sym.Matrix(self.Ric()) - 1/2 * self.RicSca() * self.metric
            self.Einstein = Einstein
        if display:
            # To display only independent non-zero components. No output means that it satisfies vacuum Einstein field equation
            Einstein = self.Einstein
            for j in range(self.dim):
                for l in range(j+1):
                    if Einstein[j,l] != 0:
                        print('Einstein[{},{}] = {}'.format(j, l, Einstein[j,l]))
        return self.Einstein


if __name__ == "__main__":
    t, r, theta, phi = sym.symbols('t r theta phi')
    coord = [t, r, theta, phi]
    m = sym.symbols('m')
    metric = sym.Matrix([[-1 + 2*m/r, 0, 0, 0], [0, (1 - 2*m/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2 * sym.sin(theta)**2]])
    test = GR(coord, metric)
    test.Eins(True)