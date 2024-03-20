#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize
import scienceplots

from solvers.ThermoDilute import ThermoDilute

DiffusionSpec = namedtuple("DiffusionSpec", "Dl1 Dl2 Ds1 Ds2")
InitSpec = namedtuple("InitSpec", "Clinf1 Clinf2 Csinf1 Csinf2")

plt.style.use('science')
fsize=plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = (fsize[0]*2,fsize[1]*2)

def logspace(start=0, end=1, npoints=100, base=10):
    length=abs(end-start)
    sign= (end>=start) * 1 + (end<start) * -1
    ls=np.logspace(np.log(1)/np.log(base),np.log(length+1)/np.log(base),num=npoints, base=base, endpoint=True)
    space=start +sign* (ls-1)
    space.sort()
    return space
    
class Stefan2p2cSolution:

    
    def __init__(self,Thermo, Diff, Init,xi,muco1,muco2):
        self.Thermo=Thermo
        self.Diff=Diff
        self.Init=Init
        self.xi=xi
        self.Xi=-xi*self.Thermo.de3
        self.muco1=muco1
        self.muco2=muco2
        
        
        
        self.Clco1=Thermo.Cl1(muco1,muco2)
        self.Clco2=Thermo.Cl2(muco1,muco2)
        self.Csco1=Thermo.Cs1(muco1,muco2)
        self.Csco2=Thermo.Cs2(muco1,muco2)
        
        self.zeta=self.Clco1/Thermo.Cl0_1
        #  self.Clco1=(1-zeta)*Thermo.Cl0_1
        #  self.Clco2=zeta*Thermo.Cl0_2
        #  self.Csco1=(1-zeta)*Thermo.Cl0_1
        #  self.Csco2=zeta*Thermo.Cs0_2
        print('Clco1={:2.3}, Clco2={:2.3}, Csco1={:2.3}, Csco2={:2.3}'.format(self.Clco1,self.Clco2,self.Csco1,self.Csco2))
        print('muco1={:2.3}, muco2={:2.3}'.format(muco1,muco2))
        

        rad=5
        D=max(max(Diff.Ds1,Diff.Ds2),max(Diff.Dl1,Diff.Dl2))
        #  minf=xi-rad*(max(Diff.Ds1,Diff.Ds2))**0.5
        #  pinf=xi+rad*(max(Diff.Dl1,Diff.Dl2))**0.5
        minf=xi-rad*D**0.5
        pinf=xi+rad*D**0.5
        nptxi=10000
        
        #  self.tls=np.linspace(minf, xi, nptxi)
        self.tls=logspace(start=xi, end=minf, npoints=nptxi)
        self.tCs1=self.fC(-self.tls, -self.xi, Init.Csinf1, self.Csco1, Diff.Ds1)
        self.tCs2=self.fC(-self.tls, -self.xi, Init.Csinf2, self.Csco2, Diff.Ds2)
    
        #  self.tll=np.linspace(xi, pinf, nptxi)
        self.tll=logspace(start=xi, end=pinf, npoints=nptxi)
        self.tCl1=self.fC(self.tll, xi, Init.Clinf1, self.Clco1, Diff.Dl1)
        self.tCl2=self.fC(self.tll, xi, Init.Clinf2, self.Clco2, Diff.Dl2)
        
        self.tC1=np.append(self.tCs1, self.tCl1)
        self.tC2=np.append(self.tCs2, self.tCl2)
        
        mulinf1=Thermo.mul1(Init.Clinf1,Init.Clinf2)
        mulinf2=Thermo.mul2(Init.Clinf1,Init.Clinf2)
        musinf1=Thermo.mus1(Init.Csinf1,Init.Csinf2)
        musinf2=Thermo.mus2(Init.Csinf1,Init.Csinf2)
        
        self.dwlinf=Thermo.dw(mulinf1,mulinf2)
        self.dwsinf=Thermo.dw(musinf1,musinf2)
        
        tmul1=Thermo.mul1(self.tCl1,self.tCl2)
        tmul2=Thermo.mul2(self.tCl1,self.tCl2)
        tmus1=Thermo.mus1(self.tCs1,self.tCs2)
        tmus2=Thermo.mus2(self.tCs1,self.tCs2)
        

        tmu1=np.append(tmus1, tmul1)
        tmu2=np.append(tmus2, tmul2)
        
        self.tmu1=tmu1
        self.tmu2=tmu2
        self.tl=np.append(self.tls, self.tll)
        
        self.dw=Thermo.dw(tmu1,tmu2)
        
        self.t0s=np.linspace(minf, 0, nptxi)
        self.t0l=np.linspace(0, pinf, nptxi)
        
        self.unstable_dws=abs((self.dw-self.dwsinf) * (self.tl < self.xi) * (self.dw < 0))
        self.unstable_dwl=abs((self.dw-self.dwlinf) * (self.tl > self.xi) * (self.dw > 0))
        
        self.excess_dws=np.trapz(self.unstable_dws, x=self.tl)
        self.excess_dwl=np.trapz(self.unstable_dwl, x=self.tl)
        
        self.total_excess_dw=abs(self.excess_dws)+abs(self.excess_dwl)
        #  print("Excess dws : ", self.excess_dws)
        #  print("Excess dwl : ", self.excess_dwl)
        
        wlinf=Thermo.wl(mulinf1, mulinf2)
        wsinf=Thermo.ws(musinf1, musinf2)
        self.twl=Thermo.wl(tmul1, tmul2)
        self.tws=Thermo.ws(tmus1, tmus2)
        self.tw=np.append(self.tws, self.twl)
        self.Xi_Omega=self.calc_Xi(self.tws, wsinf, self.twl,wlinf)
        #  self.Xi_Omega=np.trapz(self.tls, self.tws)+np.trapz(self.tll, self.twl) - np.trapz(self.tls-self.xi, wsinf*self.tls**0) - np.trapz(self.tll, wlinf * self.tll**0)
        #  self.Xi_Omega=np.trapz(self.tl, self.tw) - self.xi * (wlinf-wsinf)
        print("Xi = ",self.Xi)
        print("Xi_Omega = ",self.Xi_Omega)
        #  plt.plot(self.tl, self.tw)
        
        flinf=Thermo.fl(Init.Clinf1, Init.Clinf2)
        fsinf=Thermo.fs(Init.Csinf1, Init.Csinf2)
        self.tfl=Thermo.fl(self.tCl1, self.tCl2)
        self.tfs=Thermo.fs(self.tCs1, self.tCs2)
        self.tf=np.append(self.tfs, self.tfl)
        self.Xi_F=self.calc_Xi(self.tfs, fsinf, self.tfl, flinf)
        
        self.dfl=Thermo.fl(self.tCl1, self.tCl2)-Thermo.fs(self.tCl1, self.tCl2)
        self.dfs=Thermo.fl(self.tCs1, self.tCs2)-Thermo.fs(self.tCs1, self.tCs2)
        self.df=np.append(self.dfs, self.dfl)
        
        print("Xi_F = ",self.Xi_F)
        
        #  Xi_C1=self.calc_Xi(self.tCs1, Init.Csinf1, self.tCl1, Init.Clinf1)
        #  print("Xi_C1 = ",Xi_C1)
        #  print(self.dwd.shape)
    
    def calc_Xi(self, vs, ls, vl, ll):
        #  return (np.trapz( vs-ls, x=self.tls)+np.trapz(vl-ll, x=self.tll))-self.xi * (ll-ls)
        return (np.trapz( vs, x=self.tls)+np.trapz(vl, x=self.tll))- np.trapz(self.t0s**0*ls, x=self.t0s)- np.trapz(self.t0l**0*ll, x=self.t0l)
    # analytical composition profiles
    def fC(self,l, li, Cinf, Cco, D):
        return Cinf + (Cco-Cinf)*special.erfc(l/2/D**0.5)/special.erfc(li/2/D**0.5)
            

    def get_diffuse_profile(self, W, t):
        Wl=W/t**0.5
        tphi= (1+np.tanh(2*(self.tl-self.xi)/Wl))/2
        c1l=self.Thermo.Cl1(self.tmu1, self.tmu2)
        c1s=self.Thermo.Cs1(self.tmu1, self.tmu2)
        c1=tphi*c1l + (1-tphi)*c1s
        c2l=self.Thermo.Cl2(self.tmu1, self.tmu2)
        c2s=self.Thermo.Cs2(self.tmu1, self.tmu2)
        c2=tphi*c2l + (1-tphi)*c2s
        
        return c1,c2
        

    def show(self):
        print("Printing all solution parameters")
        plist=[("xi", self.xi),
                ("zeta", self.zeta),
                ("Csco1", self.Csco1),
                ("Csco2", self.Csco2),
                ("Clco1", self.Clco1),
                ("Clco2", self.Clco2),
                ("Xi_Omega", self.Xi_Omega),
                ("Xi_F", self.Xi_F),
                ("excess_dws", self.excess_dws),
                ("excess_dwl", self.excess_dwl),
                ("total_excess_dw", self.total_excess_dw),
            ]
        for name, val in plist:
            print(name.ljust(10," "), "=", "%7.4f" % val)


class SolverStefan2p2c:

    
    def __init__(self,Thermo, Diff, Init):
        self.Thermo=Thermo
        self.Diff=Diff
        self.Init=Init

# intermediary functions
    def u(self,x):
        return (2/np.pi**0.5)*np.exp(-x**2)/special.erfc(x)
    def ul1(self,x):
        D=self.Diff.Dl1
        return (D)**0.5 * self.u(x/2/D**0.5)
    def us1(self,x):
        D=self.Diff.Ds1
        return (D)**0.5 * self.u(-x/2/D**0.5)
    
    def ul2(self,x):
        D=self.Diff.Dl2
        return (D)**0.5 * self.u(x/2/D**0.5)
    def us2(self,x):
        D=self.Diff.Ds2
        return (D)**0.5 * self.u(-x/2/D**0.5)

    # individual equations to solve
    def f1(self,xi,muco1,muco2):
        lhs=xi*(self.Thermo.Cl1(muco1,muco2)-self.Thermo.Cs1(muco1,muco2))
        rhs=self.ul1(xi)*(self.Thermo.Cl1(muco1,muco2)-self.Init.Clinf1) + self.us1(xi)*(self.Thermo.Cs1(muco1,muco2)-self.Init.Csinf1)
        return rhs-lhs
        
    def f2(self,xi,muco1,muco2):
        lhs=xi*(self.Thermo.Cl2(muco1,muco2)-self.Thermo.Cs2(muco1,muco2))
        rhs=self.ul2(xi)*(self.Thermo.Cl2(muco1,muco2)-self.Init.Clinf2) + self.us2(xi)*(self.Thermo.Cs2(muco1,muco2)-self.Init.Csinf2)
        return rhs-lhs
        
    def f3(self,xi,muco1,muco2):
        lhs=self.Thermo.dw(muco1,muco2)
        rhs=0
        return rhs-lhs
        
        
    # transcendental function to solve for 0
    def f(self,vect):
        xi,muco1,muco2=vect
        
        res1=self.f1(xi,muco1,muco2)
        res2=self.f2(xi,muco1,muco2)
        res3=self.f3(xi,muco1,muco2)

        res=[res1,res2,res3]
        return res
        

    
    def solve(self,start, maxfev):
        print("Searching xi starting from x0 = ",start)
        sol=scipy.optimize.root(self.f,start,method='hybr')
        
        print(sol.success)
        xi,muco1,muco2=sol.x
        
        print("Found solution xi,muco1,muco2= ",sol.x)

        return Stefan2p2cSolution(self.Thermo, self.Diff, self.Init,xi,muco1,muco2)
    
#  def plotdw(sol):
    #  fig_dw, ax_dw= plt.subplots()
    
    #  dw=sol.dwd
    
    #  minl=min(sol.tl)
    #  maxl=max(sol.tl)
    #  ax_dw.set_xlim([minl,maxl])
    #  ax_dw.set_ylim([min(dw), max(dw)])
    
    #  ax_dw.plot(sol.tl,dw)

    #  ax_dw.plot([sol.xi, sol.xi], [min(dw), max(dw)], color="k")
    #  ax_dw.axvline(x=sol.xi, color="k", linewidth=1.4)
    #  ax_dw.axhline(y=0, color="k", linewidth=1.4)
    #  alp=0.1
    #  ax_dw.fill_between([minl, sol.xi],[max(dw), max(dw)], color="green", alpha=alp)
    #  ax_dw.fill_between([sol.xi, maxl],[max(dw), max(dw)], color="red", alpha=alp)
    #  ax_dw.fill_between([minl, sol.xi],[min(dw), min(dw)], color="red", alpha=alp)
    #  ax_dw.fill_between([sol.xi, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    #  ax_dw.set_xlabel("$\\lambda$")
    #  ax_dw.set_ylabel("$\\Delta^\\Phi\\omega$")
    #  plt.savefig("maugis_test_dw.pdf")
    
if __name__=="__main__":
    

    
# input parameters
    
    scale=100
    el1=0.0
    el2=0.0
    el3=5.0/scale
    
    es1=el1+np.log(25/20)
    es2=el2+np.log(20/15)
    es3=0

    Thermo=ThermoDilute(es1,es2,0,el1,el2,el3)
    Diff=DiffusionSpec(Dl1=36,Ds1=16,Dl2=4,Ds2=1)
    
    Init=InitSpec(Clinf1=23/scale,Clinf2=3/scale,Csinf1=3/scale,Csinf2=12/scale)
# print the params

    #  print('DlA={:2.3}, DlB={:2.3}, DsA={:2.3}, DsB={:2.3}'.format(float(DlA),float(DlB),float(DsA),float(DsB)))
    #  print('ClinfA={:2.3}, ClinfB={:2.3}, CsinfA={:2.3}, CsinfB={:2.3}'.format(float(ClinfA),float(ClinfB),float(CsinfA),float(CsinfB)))

    Solver = SolverStefan2p2c(Thermo, Diff, Init)
    sol1 = Solver.solve([-4.7,2,2], 1000)
    sol2 = Solver.solve([-1,2,2], 1000)
    sol3 = Solver.solve([3.1,2,2], 1000)

    fig_diag, ax_diag = plt.subplots()

    ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    ax_diag.set_xlim([0,Thermo.Cl0_1])
    ax_diag.set_ylim([0,Thermo.Cl0_2])
    ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

    ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*scale) % 5 == 0])
    ax_diag.set_xlabel("$C_1$")
    ax_diag.set_ylabel("$C_2$")
    
    p1color = ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f" % (sol1.xi)))[0].get_color()
    ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f" % (sol2.xi)), marker="x")
    ax_diag.plot(sol3.tC1, sol3.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f" % (sol3.xi)), marker="x")

    ax_diag.legend()

   
    #  plotdw(sol1)
    plt.show()





    







    
    
