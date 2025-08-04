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


        
        self.K0=K0
        self.ceq0=column(ceqA0,ceqB0)
        self.K1=K1
        self.ceq1=column(ceqA1,ceqB1)
        
        self.init_computed_params()



    def init_computed_params(self):
        self.K0inv = np.linalg.inv(self.K0) 
        self.K1inv = np.linalg.inv(self.K1) 
        
        
        
        print("Thermo params:")
        print("ceq0 =", self.ceq0)
        print("K0 = ", self.K0)
        print("K0inv = ", self.K0inv)
        print("ceq1 =", self.ceq1)
        print("K1 = ", self.K1)
        print("K1inv = ", self.K1inv)

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
        w =  - (mu1 * self.ceq0[0,0] + mu2 * self.ceq0[1,0] )
        v1 =   mu1 * self.K0inv[0,0] + mu2 * self.K0inv[0,1]
        v2 =   mu2 * self.K0inv[1,1] + mu1 * self.K0inv[0,1]
        w += - 0.5* (mu1 * v1 + mu2 * v2)
        return w
        
    def wl(self,mu1, mu2):
        w =  - (mu1 * self.ceq1[0,0] + mu2 * self.ceq1[1,0] )
        v1 =   mu1 * self.K1inv[0,0] + mu2 * self.K1inv[0,1]
        v2 =   mu2 * self.K1inv[1,1] + mu1 * self.K1inv[0,1]
        w += - 0.5* (mu1 * v1 + mu2 * v2)
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
        # ~ else:
            # ~ print("found valid point")
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
                        else:
                            self.add_point(muA, muB)
                        


if __name__ == "__main__":
    
    # cas test de base
    K0 = np.array([[1, 0],[0, 1]])
    # peut enlever le 2* ci dessous, devrait donner deux droites
    # (je crois que ca correspond au cas test basique que werner utilisait dans sa NT)
    # avec le 2*, deux cercles concentriques
    K1 = 2*np.array([[1, 0],[0, 1]]) 
    ceqA0 = 0.1
    ceqB0 = 0.1
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    
    # ~ from ThermoK import ThermoK, mat
    # params capu MoO3 = 1.5 %mol
    K0 = mat(1.4207, 1.9431, 1.6123)
    K1 = mat(6.4636, 22.7486, -9.8138)
    ceqA0 = 0.600056
    ceqB0 = 0.392136
    ceqA1 = 0.023910
    ceqB1 = 0.510693
    ciniA = 0.591
    ciniB = 0.394
    
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
    K0 = mat(2.3345,3.0121,2.6463)
    K1 = mat(10.4042,44.4530,-21.0754)
    ceqA0 = 0.610548
    ceqB0 = 0.381839
    ceqA1 = 0.020803
    ceqB1 = 0.509110
    ciniA = 0.582
    ciniB = 0.388
    
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    TK.compute_diagram_from_c()
    
    plt.plot([0,1,0,0], [0,0,1,0], color="k")
    
    Na2O=1
    plt.plot(TK.front0[Na2O], TK.front0[1-Na2O], color="b", label="Frontiere phase 0")
    plt.plot(TK.front1[Na2O], TK.front1[1-Na2O], color="r", label="Frontiere phase 1")
    
    
    plt.scatter([ciniB], [ciniA], color="purple", label="Cini")
    plt.scatter([ceqB0, ceqB1], [ceqA0, ceqA1], color="green", label="Ceq")
    
    for tl in TK.tielines:
        plt.plot(tl[Na2O], tl[1-Na2O], color = "gray", zorder= -10)
    plt.xlabel("Na2O")
    plt.ylabel("SiO2")
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.legend()
    plt.show()