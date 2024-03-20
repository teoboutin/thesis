#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize
import ternary


class ThermoDilute:

    
    def __init__(self,es1,es2,es3,el1,el2,el3):

        self.thermoParams={}
        self.thermoParams[0]={}
        self.thermoParams[1]={}
        self.thermoParams[0]["A"]=es1
        self.thermoParams[0]["B"]=es2
        self.thermoParams[0]["C"]=es3
        self.thermoParams[1]["A"]=el1
        self.thermoParams[1]["B"]=el2
        self.thermoParams[1]["C"]=el3
        
        self.es1=es1
        self.es2=es2
        self.es3=es3
        self.el1=el1
        self.el2=el2
        self.el3=el3
        
        self.init_computed_params()



    def init_computed_params(self):
        self.de3=self.el3-self.es3
        
        self.k1=np.exp(self.es1-self.el1)
        self.k2=np.exp(self.es2-self.el2)
        
        
        self.Cl0_1=self.k1*self.de3/(self.k1-1)
        self.Cs0_1=self.de3/(self.k1-1)
        self.Cl0_2=self.k2*self.de3/(self.k2-1)
        self.Cs0_2=self.de3/(self.k2-1)
        
        
        print("Thermo params:")
        print('el1={:2.5}, el2={:2.5}, el3={:2.5}'.format(float(self.el1),float(self.el2),float(self.el3)))
        print('es1={:2.5}, es2={:2.5}, es3={:2.5}'.format(float(self.es1),float(self.es2),float(self.es3)))
        print('Cl0_1={:2.3}, Cs0_1={:2.3}, Cl0_2={:2.3}, Cs0_2={:2.3}'.format(self.Cl0_1,self.Cs0_1,self.Cl0_2,self.Cs0_2))

    def fs(self,c1, c2):
        return c1*(self.es1-1+np.log(c1)) + c2*(self.es2-1+np.log(c2)) + self.es3
    def mus1(self,c1,c2):
        return self.es1+np.log(c1)
    def mus2(self,c1,c2):
        return self.es2+np.log(c2)
        
    def Cs1(self, mu1, mu2):
        return np.exp(mu1-self.es1)
    def Cs2(self, mu1, mu2):
        return np.exp(mu2-self.es2)
        
    def fl(self,c1, c2):
        return c1*(self.el1-1+np.log(c1)) + c2*(self.el2-1+np.log(c2)) + self.el3
        
    def mul1(self,c1,c2):
        return self.el1+np.log(c1)
        
    def mul2(self,c1,c2):
        return self.el2+np.log(c2)
        

    def Cl1(self,mu1, mu2):
        return np.exp(mu1-self.el1)
    def Cl2(self,mu1, mu2):
        return np.exp(mu2-self.el2)
        
    def ws(self,mu1, mu2):
        #  return self.fs(self.Cs1(mu1, mu2), self.Cs2(mu1, mu2)) - mu1*self.Cs1(mu1, mu2) - mu2*self.Cs2(mu1, mu2)
        return self.es3 - self.Cs1(mu1, mu2) - self.Cs2(mu1, mu2)
    def wl(self,mu1, mu2):
        #  return self.fl(self.Cl1(mu1, mu2), self.Cl2(mu1, mu2)) - mu1*self.Cl1(mu1, mu2) - mu2*self.Cl2(mu1, mu2)
        return self.el3 - self.Cl1(mu1, mu2) - self.Cl2(mu1, mu2)
        
    def dw(self, mu1, mu2):
        return self.wl(mu1, mu2) - self.ws(mu1, mu2)
        
    def df(self, c1, c2):
        return self.fl(c1, c2) - self.fs(c1, c2)
        
