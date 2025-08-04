#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize
# ~ import ternary


def row(a,b):
    r=np.zeros((1,2))
    r[0,0]=a
    r[0,1]=b
    return r

def column(a,b):
    c=np.zeros((2,1))
    c[0,0]=a
    c[1,0]=b
    return c
def mat(a,b,c):
    return np.array([[a, c],[c, b]])
class ThermoK:

    
    def __init__(self,K0,ceqA0,ceqB0,K1,ceqA1,ceqB1):

        if not np.all(np.linalg.eigvals(K0) > 0):
            raise RuntimeError("Matrix K0 should be positive definite")
        
        if not K0[1,0] == K0[0,1]:
            raise RuntimeError("Matrix K0 should be symetric")
            
        if not np.all(np.linalg.eigvals(K1) > 0):
            raise RuntimeError("Matrix K1 should be positive definite")
            
        if not K1[1,0] == K1[0,1]:
            raise RuntimeError("Matrix K1 should be symetric")
        
        self.K0=K0
        self.ceq0=column(ceqA0,ceqB0)
        self.K1=K1
        self.ceq1=column(ceqA1,ceqB1)
        
        self.init_computed_params()



    def init_computed_params(self):
        self.K0inv = np.linalg.inv(self.K0) 
        self.K1inv = np.linalg.inv(self.K1) 
        
        
        print("Thermo params:")
        rm=5
        print("ceq0".center(rm)+ " | " + "K0".center(2*rm) + " | " + "K0inv".center(2*rm))
        print( str(self.ceq0[0,0]).center(rm) + " | " + str(self.K0[0,0]).center(rm) + str(self.K0[0,1]).center(rm) + " | " + str(self.K0inv[0,0]).center(rm) + str(self.K0inv[0,1]).center(rm))
        print( str(self.ceq0[1,0]).center(rm) + " | " + str(self.K0[1,0]).center(rm) + str(self.K0[1,1]).center(rm) + " | " + str(self.K0inv[1,0]).center(rm) + str(self.K0inv[1,1]).center(rm))
        print("ceq1".center(rm)+ " | " + "K1".center(2*rm) + " | " + "K1inv".center(2*rm))
        print( str(self.ceq1[0,0]).center(rm) + " | " + str(self.K1[0,0]).center(rm) + str(self.K1[0,1]).center(rm) + " | " + str(self.K1inv[0,0]).center(rm) + str(self.K1inv[0,1]).center(rm))
        print( str(self.ceq1[1,0]).center(rm) + " | " + str(self.K1[1,0]).center(rm) + str(self.K1[1,1]).center(rm) + " | " + str(self.K1inv[1,0]).center(rm) + str(self.K1inv[1,1]).center(rm))

        # print the eigen props of deltaKinv
        res=np.linalg.eig(self.K1inv - self.K0inv)
        print("eigenvalues of deltaKinv")
        print(res.eigenvalues)
        print("eigenvectors of deltaKinv")
        print(res.eigenvectors)

        # nature of the phase diagram:
        # the function dw is a quadratic function of 2 variables whose hessian is H=(Kinv1 - Kinv0)
        # the equation dw = 0 takes the isocontour at 0 of this function
        # Then, the sign of  det(H) gives the kind of graph of the quadratic function and its isocontours
        # see https://en.wikipedia.org/wiki/Quadratic_function#Bivariate_(two_variable)_quadratic_function
        # this is the same as taking a conic by instersecting a plane with a cone
        # det(H) > 0 : dw is an elliptic paraboloid and the frontiers will be ellipses
        # det(H) < 0 : dw is an hyperbolic paraboloid and the frontiers will be hyperbolas
        # det(H) = 0 : two possibilities:
        #          - if H = 0: the curve is simply a line (mu orthogonal to ceq1 - ceq0)
        #          - if H !=0: dw is a parabolic cylinder function, and the isocontour is a parabola. An example is given below, where H has a single nonzero eigenvalue.
        
        # more details in https://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections
        # actually, there should also be the case of two intersecting lines possible (degenerated hyberbola)
        # and single point (degenerate ellipse), possible here if delta_ceq=0
        
        # actually, the diagram may be computable analytically:
        # see https://en.wikipedia.org/wiki/Conic_section#General_Cartesian_form
        
        H=self.K1inv - self.K0inv
        rhs =  -2 * (self.ceq1 - self.ceq0)
        
        A = H[0,0]
        B = H[1,0] + H[0,1]
        C = H[1,1]
        D = -rhs[0,0]
        E = -rhs[1,0]
        F=0
        Aq = [ [ A  , B/2, D/2 ],
               [ B/2, C  , E/2 ],
               [ D/2, E/2, F   ] ]
        A33 = [ [ A  , B/2 ],
                [ B/2, C   ] ]
        
        
        print("Matrix of the quadratic function describing the diagram (Aq):")
        print(Aq[0])
        print(Aq[1])
        print(Aq[2])
        detAq = np.linalg.det(Aq)
        detA33 = np.linalg.det(A33)
        eps=1e-10
        if abs(detAq)>eps:
            if abs(detA33) < eps:
                print("diagram should be made of parabolas")
            elif detA33 < 0:
                print("diagram should be made of hyperbolas")
            elif detA33 > 0:
                print("diagram should be made of ellipses")
            else:
                raise RuntimeError("wrong case of detAq ??")
        else:
            print("degenerate conic section")
            if abs(detA33) < eps:
                print("diagram should be made of straight lines")
            elif detA33 < 0:
                print("diagram should be made of intersecting lines")
            elif detA33 > 0:
                print("diagram should be made of a single point (should not happen)")
            else:
                raise RuntimeError("wrong case of detAq ??")
                
        # old way (that i came up myself...)
        # ~ det = res.eigenvalues[0] * res.eigenvalues[1]
        
        # ~ if det>0:
            # ~ print("deltaKinv has positive determinant, phase diagram will be composed of ellipses")
        # ~ if det<0:
            # ~ print("deltaKinv has negative determinant, phase diagram will be composed of hyperbolas")
        # ~ if det==0 and np.all(res.eigenvalues == 0 ):
            # ~ print("deltaKinv is zero, phase diagram will be made of straigth lines")
        # ~ if det==0 and np.any(res.eigenvalues != 0 ):
            # ~ print("deltaKinv has null determinant and is non zero, phase diagram will be made of parabolas")

    def fs(self,cA, cB):
        return 0 * cA
        
    # vector of mu
    def mu0(self,cA, cB):
        return np.matmul(self.K0, column(cA, cB) - self.ceq0)
    
    # for external interface
    def mus1(self,cA, cB):
        return cA * self.K0[0,0] + cB * self.K0[0,1]
        
    def mus2(self,cA, cB):
        return cB * self.K0[1,1] + cA * self.K0[0,1]

    # vector of C
    def c0(self,muA, muB):
        return self.ceq0 + np.matmul(self.K0inv, column(muA, muB) )
        
    def Cs1(self, mu1, mu2):
        return self.ceq0[0,0] + mu1 * self.K0inv[0,0] + mu2 * self.K0inv[0,1]
        
    def Cs2(self, mu1, mu2):
        return self.ceq0[1,0] + mu2 * self.K0inv[1,1] + mu1 * self.K0inv[0,1]
        
    def fl(self,c1, c2):
        return 0 * c1
        
    # vector of mu
    def mu1(self,cA, cB):
        return np.matmul(self.K1, column(cA, cB) - self.ceq1)
        
    def mul1(self,cA, cB):
        return cA * self.K1[0,0] + cB * self.K1[0,1]
        
    def mul2(self,cA, cB):
        return cB * self.K1[1,1] + cA * self.K1[0,1]
        

    # vector of C
    def c1(self,muA, muB):
        return self.ceq1 + np.matmul(self.K1inv, column(muA, muB) )
        
    def Cl1(self,mu1, mu2):
        return self.ceq1[0,0] + mu1 * self.K1inv[0,0] + mu2 * self.K1inv[0,1]
        
    def Cl2(self,mu1, mu2):
        return self.ceq1[1,0] + mu2 * self.K1inv[1,1] + mu1 * self.K1inv[0,1]
        
    def ws(self,mu1, mu2):
        v1 =   mu1 * self.K0inv[0,0] + mu2 * self.K0inv[0,1]
        v2 =   mu2 * self.K0inv[1,1] + mu1 * self.K0inv[0,1]
        w =  - (mu1 * self.ceq0[0,0] + mu2 * self.ceq0[1,0] ) - 0.5* (mu1 * v1 + mu2 * v2)
        return w
        
    def wl(self,mu1, mu2): 
        v1 =   mu1 * self.K1inv[0,0] + mu2 * self.K1inv[0,1]
        v2 =   mu2 * self.K1inv[1,1] + mu1 * self.K1inv[0,1]
        w =   - (mu1 * self.ceq1[0,0] + mu2 * self.ceq1[1,0] ) - 0.5* (mu1 * v1 + mu2 * v2)
        return w
        
    def dw(self, mu1, mu2):
        return self.wl(mu1, mu2) - self.ws(mu1, mu2)
        
    def df(self, c1, c2):
        return self.fl(c1, c2) - self.fs(c1, c2)
        
    def init_diagram(self):
        
        self.tielines=[]
        self.front0=[[], []]
        self.front1=[[], []]
        
        
    def add_point(self,muA, muB):
        c0, c1 = self.c0(muA, muB), self.c1(muA, muB)
        c0A, c0B = c0[0,0], c0[1,0]
        c1A, c1B = c1[0,0], c1[1,0]
        
        if c0A < 0 or c0B < 0 or c1A <0 or c1B <0:
            return
        if c0A + c0B > 1 or c1A + c1B > 1:
            return
        self.front0[0].append(c0A)
        self.front0[1].append(c0B)
        self.front1[0].append(c1A)
        self.front1[1].append(c1B)
        cA = [c0A, c1A]
        cB = [c0B, c1B]
        self.tielines.append([ cA, cB ])
        
    def compute_diagram(self):
        self.init_diagram()
        
        muA_max=-np.infty
        muA_min=np.infty
        
        muB_max=-np.infty
        muB_min=np.infty
        n=500
        
        for cA in np.linspace(0,1,n):
            for cB in np.linspace(0,1-cA,int(n * (1-cA))):
                muAs = [ self.mul1(cA, cB) ,self.mus1(cA, cB) ]
                muBs = [ self.mul2(cA, cB) ,self.mus2(cA, cB) ]
                muA_max = max(muA_max, *muAs)
                muA_min = min(muA_min, *muAs)
                muB_max = max(muB_max, *muBs)
                muB_min = min(muB_min, *muBs)
        
        print("\nmuA range:", muA_min, muA_max)
        print("\nmuB range:", muB_min, muB_max)
        
        
        nA=100
        nB=0
        for muA in np.linspace(muA_min,muA_max,nA):
            muB = None
            found=False
            fun = lambda mu : self.dw(muA, mu)
            
            
            sol = scipy.optimize.root_scalar(fun, x0 = muA, method='newton')
            if sol.converged:
                muB = sol.root
                print("found root from newton")
                found=True
                self.add_point(muA, muB)
            # ~ else:
            
            # must bruteforce like this to find multiple roots
            # may be better to do this on coarse space then refine with newton
            for muBloc in np.linspace(muB_min,muB_max, nB):
                dw = self.dw(muA, muBloc)
                
                if abs(dw) < 1e-3:
                    muB = muBloc
                    found=True
                    sol = scipy.optimize.root_scalar(fun, x0 = muBloc, method='newton')
                    if sol.converged:
                        muB = sol.root
                        
                    # ~ fun2 = lambda mu : self.dw(mu[0], mu[1])
                    # ~ sol = scipy.optimize.root_scalar(fun2, x0 = [muA, muBloc], method='newton')
                    # ~ if sol.converged:
                        # ~ muB = sol.root
                if found:
                    self.add_point(muA, muB)
        
        return self.tielines, self.front0, self.front1

    def compute_diagram_from_c(self):
        self.init_diagram()
        
        n=200
        
        for c0A in np.linspace(0,1,n):
            for c0B in np.linspace(0,1-c0A,int(n * (1-c0A))):
                mu0 = self.mu0(c0A, c0B)
                muA, muB = mu0[:,0]
                # find c1A, c1B with same mu
                fun = lambda c: ((mu0 - self.mu1(c[0], c[1]))[:,0])
                sol = scipy.optimize.root(fun, x0 = [c0A, c0B])
                if sol.success:
                    c1A, c1B = sol.x
                    dw = self.dw( muA, muB)
                    if abs(dw)<1e-2:
                        fun2 = lambda mu: (self.dw(mu[0], mu[1]))**2
                        sol = scipy.optimize.minimize(fun2, x0 = [muA, muB], method='Nelder-Mead')
                        if sol.success:
                            # ~ print("may have improved")
                            self.add_point(*sol.x)
                        # ~ else:
                            # ~ self.add_point(muA, muB)
                        


def fit_ellipse(XX, YY):
    # see https://stackoverflow.com/questions/47873759/how-to-fit-a-2d-ellipse-to-given-points
    # Extract x coords and y coords of the data as column vectors
    X = np.array([[vv] for vv in XX])
    Y = np.array([[vv] for vv in YY])
    
    # Formulate and solve the least squares problem ||Ax - b ||^2
    A = np.hstack([X**2, X * Y, Y**2, X, Y])
    b = np.ones_like(X)
    x = np.linalg.lstsq(A, b)[0].squeeze()
    
    # Print the equation of the ellipse in standard form
    print('The ellipse is given by {0:.3}x^2 + {1:.3}xy + {2:.3}y^2 + {3:.3}x + {4:.3}y = 1'.format(x[0], x[1],x[2],x[3],x[4]))
    return x, [min(XX), max(XX)], [min(YY), max(YY)]

if __name__ == "__main__":
    
    
    # diagram with parabolas: deltaKinv has exaclty 1 non zero eigenvalue
    K0 = np.linalg.inv(mat(1,2,1))
    K1 = np.linalg.inv(mat(1,3,1))
    K0 = mat(2,1,0)
    K1 = mat(2,4,0)
    ceqA0 = 0.3
    ceqB0 = 0.3
    ceqA1 = 0.4
    ceqB1 = 0.4
    
    # ~ # diagram with straigth lines: deltaKinv = 0, degenerate parabolas
    # ~ K0 = mat(1,1,0)
    # ~ K1 = mat(1,1,0)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # diagram with circles deltaKinv != 0, eigenvalues are equals
    # ~ K0 = mat(1,1,0)
    # ~ K1 = 2*mat(1,1,0)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # diagram with ellipses: det(deltaKinv) > 0
    # ~ K0 = mat(1,1,0)
    # ~ K1 = mat(2,2,0.2)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # diagram with a single point: degenerate ellipses. Numerical resolution is very unstable
    # ~ K0 = mat(1,1,0)
    # ~ K1 = 2*mat(1,1,0)
    # ~ ceqA0 = 0.2
    # ~ ceqB0 = 0.2
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # diagram with hyperbolas: det(deltaKinv) < 0. managed to get both section in the gibbs triangle, yay !
    # ~ K0 = mat(2,1,0)
    # ~ K1 = mat(1,2,1)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # diagram with intersecting lines (degenerate hyperbola)
    # ~ K0 = mat(2,1,0)
    # ~ K1 = mat(1,2,0)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    
    # ~ from ThermoK import ThermoK, mat
    # params capu MoO3 = 1.5 %mol
    # ~ K0 = mat(1.4207, 1.9431, 1.6123)
    # ~ K1 = mat(6.4636, 22.7486, -9.8138)
    # ~ ceqA0 = 0.600056
    # ~ ceqB0 = 0.392136
    # ~ ceqA1 = 0.023910
    # ~ ceqB1 = 0.510693
    # ~ ciniA = 0.591
    # ~ ciniB = 0.394
    
    # params capu MoO3 = 2 %mol
    # ~ K0 = mat(1.4705,1.9490,1.6757)
    # ~ K1 = mat(6.9739,30.2869,-12.6910)
    # ~ ceqA0 = 0.603498
    # ~ ceqB0 = 0.38876
    # ~ ceqA1 = 0.022844
    # ~ ceqB1 = 0.510151
    # ~ ciniA = 0.588
    # ~ ciniB = 0.392
    
    # params capu MoO3 = 3 %mol
    # ~ K0 = mat(2.3345,3.0121,2.6463)
    # ~ K1 = mat(10.4042,44.4530,-21.0754)
    # ~ ceqA0 = 0.610548
    # ~ ceqB0 = 0.381839
    # ~ ceqA1 = 0.020803
    # ~ ceqB1 = 0.509110
    ciniA = 0.582
    ciniB = 0.388
    
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    TK.compute_diagram_from_c()
    
    plt.plot([0,1,0,0], [0,0,1,0], color="k")
    
    Na2O=1
    plt.scatter(TK.front0[Na2O], TK.front0[1-Na2O], color="b", label="Frontiere phase 0")
    plt.scatter(TK.front1[Na2O], TK.front1[1-Na2O], color="r", label="Frontiere phase 1")
    

    # ~ x, rx, ry = fit_ellipse(TK.front0[Na2O], TK.front0[1-Na2O])
    
    # ~ rr = 10000
    # ~ x_coord = np.linspace(rx[0],rx[1],rr)
    # ~ y_coord = np.linspace(ry[0],ry[1],rr)
    # ~ X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
    # ~ Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
    # ~ plt.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('red'), linewidths=2)

    # ~ x, rx, ry = fit_ellipse(TK.front1[Na2O], TK.front1[1-Na2O])
    
    # ~ rr = 10000
    # ~ x_coord = np.linspace(rx[0],rx[1],rr)
    # ~ y_coord = np.linspace(ry[0],ry[1],rr)
    # ~ X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
    # ~ Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
    # ~ plt.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('blue'), linewidths=2)
    
    plt.scatter([ciniB], [ciniA], color="purple", label="Cini")
    plt.plot([ceqB0, ceqB1], [ceqA0, ceqA1], color="green", label="Ceq", marker="o")
    
    for tl in TK.tielines:
        plt.plot(tl[Na2O], tl[1-Na2O], color = "gray", zorder= -10)
    plt.xlabel("Na2O")
    plt.ylabel("SiO2")
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.legend()
    plt.savefig("diagramK.pdf")
    plt.show()