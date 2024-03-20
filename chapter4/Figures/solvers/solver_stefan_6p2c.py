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
plt.rcParams["figure.raise_window"] = False


def linlogspace(start, end, npoints, base=10):
    return np.logspace(np.log(start)/np.log(base),np.log(end)/np.log(base),num=npoints, base=base)
class Stefan6p2cSolution:

    
    def __init__(self,Thermo, Diff, Init,sol):
        
        self.Np=6
        self.Thermo=Thermo
        self.Diff=Diff
        self.Init=Init
        
        xi12,xi23,xi34,xi45,xi56,zeta12,zeta23,zeta34,zeta45,zeta56=sol.x
    
    
        D1_1=Diff.Ds1
        D1_2=Diff.Ds2
        
        D2_1=Diff.Dl1
        D2_2=Diff.Dl2
        
        D3_1=Diff.Ds1
        D3_2=Diff.Ds2
        
        D4_1=Diff.Dl1
        D4_2=Diff.Dl2
        
        D5_1=Diff.Ds1
        D5_2=Diff.Ds2
        
        D6_1=Diff.Dl1
        D6_2=Diff.Dl2
        CompCoeffs=dict()
        
        CompCoeffs["D1A"]=Diff.Ds1
        CompCoeffs["D1B"]=Diff.Ds2
        
        CompCoeffs["D2A"]=Diff.Dl1
        CompCoeffs["D2B"]=Diff.Dl2
        
        CompCoeffs["D3A"]=Diff.Ds1
        CompCoeffs["D3B"]=Diff.Ds2
        
        CompCoeffs["D4A"]=Diff.Dl1
        CompCoeffs["D4B"]=Diff.Dl2
        
        CompCoeffs["D5A"]=Diff.Ds1
        CompCoeffs["D5B"]=Diff.Ds2
        
        CompCoeffs["D6A"]=Diff.Dl1
        CompCoeffs["D6B"]=Diff.Dl2
        
        self.xi12=xi12
        self.xi23=xi23
        self.xi34=xi34
        self.xi45=xi45
        self.xi56=xi56
        self.Xi=-self.Thermo.de3*(xi12+xi34-xi23)
        
        
        self.zeta12=zeta12
        self.zeta23=zeta23
        self.zeta34=zeta34
        self.zeta45=zeta45
        self.zeta56=zeta56
        
        self.xi_mean=(xi12+xi34-xi23)
        
        CompCoeffs["C112A"]=(1-zeta12)*self.Thermo.Cs0_1
        CompCoeffs["C112B"]=zeta12*self.Thermo.Cs0_2
        
        CompCoeffs["C212A"]=(1-zeta12)*self.Thermo.Cl0_1
        CompCoeffs["C212B"]=zeta12*self.Thermo.Cl0_2
        
        CompCoeffs["C223A"]=(1-zeta23)*self.Thermo.Cl0_1
        CompCoeffs["C223B"]=zeta23*self.Thermo.Cl0_2
        
        CompCoeffs["C323A"]=(1-zeta23)*self.Thermo.Cs0_1
        CompCoeffs["C323B"]=zeta23*self.Thermo.Cs0_2
        
        CompCoeffs["C334A"]=(1-zeta34)*self.Thermo.Cs0_1
        CompCoeffs["C334B"]=zeta34*self.Thermo.Cs0_2
        
        CompCoeffs["C434A"]=(1-zeta34)*self.Thermo.Cl0_1
        CompCoeffs["C434B"]=zeta34*self.Thermo.Cl0_2
        
        CompCoeffs["C445A"]=(1-zeta45)*self.Thermo.Cl0_1
        CompCoeffs["C445B"]=zeta45*self.Thermo.Cl0_2
        
        CompCoeffs["C545A"]=(1-zeta45)*self.Thermo.Cs0_1
        CompCoeffs["C545B"]=zeta45*self.Thermo.Cs0_2
        
        CompCoeffs["C556A"]=(1-zeta56)*self.Thermo.Cs0_1
        CompCoeffs["C556B"]=zeta56*self.Thermo.Cs0_2
        
        CompCoeffs["C656A"]=(1-zeta56)*self.Thermo.Cl0_1
        CompCoeffs["C656B"]=zeta56*self.Thermo.Cl0_2
        
    
        # first phase
        CompCoeffs["A1A"]=Init.Csinf1
        CompCoeffs["B1A"]=(CompCoeffs["C112A"]-Init.Csinf1)/special.erfc(-xi12/2/D1_1**0.5)
        CompCoeffs["A1B"]=Init.Csinf2
        CompCoeffs["B1B"]=(CompCoeffs["C112B"]-Init.Csinf2)/special.erfc(-xi12/2/D1_2**0.5)
        
        # second phase
        CompCoeffs["B2A"]=(CompCoeffs["C223A"]-CompCoeffs["C212A"])/(special.erfc(xi23/2/D2_1**0.5)-special.erfc(xi12/2/D2_1**0.5))
        #  CompCoeffs["B2A"]=(CompCoeffs["C223A"]-CompCoeffs["C212A"])/(special.erfc(-xi12/2/D2_1**0.5)-special.erfc(-xi23/2/D2_1**0.5))
        CompCoeffs["A2A"]=CompCoeffs["C212A"]-CompCoeffs["B2A"]*special.erfc(xi12/2/D2_1**0.5)
        CompCoeffs["A2A2"]=CompCoeffs["C223A"]-CompCoeffs["B2A"]*special.erfc(xi23/2/D2_1**0.5)
    

        CompCoeffs["B2B"]=(CompCoeffs["C223B"]-CompCoeffs["C212B"])/(special.erfc(xi23/2/D2_2**0.5)-special.erfc(xi12/2/D2_2**0.5))
        #  CompCoeffs["B2B"]=(CompCoeffs["C223B"]-CompCoeffs["C212B"])/(special.erfc(-xi12/2/D2_2**0.5)-special.erfc(-CompCoeffs["xi23"]/2/D2_2**0.5))
        CompCoeffs["A2B"]=CompCoeffs["C212B"]-CompCoeffs["B2B"]*special.erfc(xi12/2/D2_2**0.5)
        CompCoeffs["A2B2"]=CompCoeffs["C223B"]-CompCoeffs["B2B"]*special.erfc(xi23/2/D2_2**0.5)
    
        
        # third phase
        CompCoeffs["B3A"]=(CompCoeffs["C323A"]-CompCoeffs["C334A"])/(special.erfc(-xi23/2/D3_1**0.5)-special.erfc(-xi34/2/D3_1**0.5))
        CompCoeffs["A3A"]=CompCoeffs["C334A"]-CompCoeffs["B3A"]*special.erfc(-xi34/2/D3_1**0.5)
        

        CompCoeffs["B3B"]=(CompCoeffs["C323B"]-CompCoeffs["C334B"])/(special.erfc(-xi23/2/D3_2**0.5)-special.erfc(-xi34/2/D3_2**0.5))
        CompCoeffs["A3B"]=CompCoeffs["C334B"]-CompCoeffs["B3B"]*special.erfc(-xi34/2/D3_2**0.5)
    
        # fourth phase
        CompCoeffs["B4A"]=(CompCoeffs["C434A"]-CompCoeffs["C445A"])/(special.erfc(xi45/2/D4_1**0.5)-special.erfc(xi34/2/D4_1**0.5))
        CompCoeffs["A4A"]=CompCoeffs["C445A"]-CompCoeffs["B4A"]*special.erfc(xi45/2/D4_1**0.5)
        

        CompCoeffs["B4B"]=(CompCoeffs["C434B"]-CompCoeffs["C445B"])/(special.erfc(xi45/2/D3_2**0.5)-special.erfc(xi34/2/D3_2**0.5))
        CompCoeffs["A4B"]=CompCoeffs["C445B"]-CompCoeffs["B4B"]*special.erfc(xi45/2/D3_2**0.5)
    
        # fifth phase
        CompCoeffs["B5A"]=(CompCoeffs["C545A"]-CompCoeffs["C556A"])/(special.erfc(-xi45/2/D5_1**0.5)-special.erfc(-xi56/2/D5_1**0.5))
        CompCoeffs["A5A"]=CompCoeffs["C556A"]-CompCoeffs["B5A"]*special.erfc(-xi56/2/D5_1**0.5)
        

        CompCoeffs["B5B"]=(CompCoeffs["C545B"]-CompCoeffs["C556B"])/(special.erfc(-xi45/2/D5_2**0.5)-special.erfc(-xi56/2/D5_2**0.5))
        CompCoeffs["A5B"]=CompCoeffs["C556B"]-CompCoeffs["B5B"]*special.erfc(-xi56/2/D5_2**0.5)
    
        # sixth phase
        CompCoeffs["A6A"]=Init.Clinf1
        CompCoeffs["B6A"]=(CompCoeffs["C656A"]-Init.Clinf1)/special.erfc(xi56/2/D6_1**0.5)
        CompCoeffs["A6B"]=Init.Clinf2
        CompCoeffs["B6B"]=(CompCoeffs["C656B"]-Init.Clinf2)/special.erfc(xi56/2/D6_2**0.5)
    
        self.CompCoeffs=CompCoeffs
        
      
        rad=5
        minf=xi12-rad*(max(Diff.Ds1,Diff.Ds2))**0.5
        pinf=xi34+rad*(max(Diff.Dl1,Diff.Dl2))**0.5
        nptxi=10000
        
        def flux_left(x,y,pos):
            return (y[pos]-y[pos-1])/(x[pos]-x[pos-1])
        def flux_right(x,y,pos):
            return (y[pos+1]-y[pos])/(x[pos+1]-x[pos])
            
        
        limits=[minf,xi12,xi23,xi34,xi45,xi56,pinf]
        
        self.tl=dict(tot=np.array([]))
        self.tCA=dict(tot=np.array([]))
        self.tCB=dict(tot=np.array([]))
        
        for i in range(1,self.Np+1):
            self.tl[i]=np.linspace(limits[i-1],limits[i], nptxi)
            self.tCA[i]=self.fC( self.tl[i], i, "A")
            self.tCB[i]=self.fC( self.tl[i], i, "B")
            
            self.tl["tot"]=np.append(self.tl["tot"], self.tl[i])
            self.tCA["tot"]=np.append(self.tCA["tot"], self.tCA[i])
            self.tCB["tot"]=np.append(self.tCB["tot"], self.tCB[i])
        
        # plot first phase solution
        self.tl1=np.linspace(minf, xi12, nptxi)
        tC1_1=self.fC1A(self.tl1)
        tC1_2=self.fC1B(self.tl1)
        

        tmu1_1=Thermo.mus1(tC1_1,tC1_2)
        tmu1_2=Thermo.mus2(tC1_1,tC1_2)

        # plot second phase sol
        self.tl2=np.linspace(xi12, xi23, nptxi)
        tC2_1=self.fC2A(self.tl2)
        tC2_2=self.fC2B(self.tl2)

        #  print(xi12/2 * (tC1_1[-1]-tC2_1[0]))
        
        tmu2_1=Thermo.mul1(tC2_1,tC2_2)
        tmu2_2=Thermo.mul2(tC2_1,tC2_2)
        
        # plot third phase sol
        self.tl3=np.linspace(xi23, xi34, nptxi)
        tC3_1=self.fC3A(self.tl3)
        tC3_2=self.fC3B(self.tl3)

        tmu3_1=Thermo.mus1(tC3_1,tC3_2)
        tmu3_2=Thermo.mus2(tC3_1,tC3_2)
    
        # plot fourth phase sol
        self.tl4=np.linspace(xi34, pinf, nptxi)
        tC4_1=self.fC4A(self.tl4)
        tC4_2=self.fC4B(self.tl4)
        
        tmu4_1=Thermo.mul1(tC4_1,tC4_2)
        tmu4_2=Thermo.mul2(tC4_1,tC4_2)


        print("Checking stefan condtion at xi12")
        print("C1:",D2_1 * flux_right(self.tl2, tC2_1, 0) - D1_1 * flux_left(self.tl1, tC1_1, -1)-xi12/2 * (tC1_1[-1]-tC2_1[0]))
        print("C2:",D2_2 * flux_right(self.tl2, tC2_2, 0) - D1_2 * flux_left(self.tl1, tC1_2, -1)-xi12/2 * (tC1_2[-1]-tC2_2[0]))
        
        #  print("Checking stefan condtion at xi23")
        #  print("C1:",D3_1 * flux_right(self.tl3, tC3_1, 0) - D2_1 * flux_left(self.tl2, tC2_1, -1)-xi23/2 * (tC2_1[-1]-tC3_1[0]))
        #  print("C2:",D3_2 * flux_right(self.tl3, tC3_2, 0) - D2_2 * flux_left(self.tl2, tC2_2, -1)-xi23/2 * (tC2_2[-1]-tC3_2[0]))
        
        #  print("Checking stefan condtion at xi34")
        #  print("C1:",D4_1 * flux_right(self.tl4, tC4_1, 0) - D3_1 * flux_left(self.tl3, tC3_1, -1)-xi34/2 * (tC3_1[-1]-tC4_1[0]))
        #  print("C2:",D4_2 * flux_right(self.tl4, tC4_2, 0) - D3_2 * flux_left(self.tl3, tC3_2, -1)-xi34/2 * (tC3_2[-1]-tC4_2[0]))
        
        # append solutions
        self.tl=np.append(self.tl1,[self.tl2,self.tl3,self.tl4])
        self.tC1=np.append(tC1_1,[tC2_1,tC3_1,tC4_1])
        self.tC2=np.append(tC1_2,[tC2_2,tC3_2,tC4_2])

        tmu1=np.append(tmu1_1,[tmu2_1,tmu3_1,tmu4_1])
        tmu2=np.append(tmu1_2,[tmu2_2,tmu3_2,tmu4_2])
        
        self.dw1=Thermo.dw(tmu1_1,tmu1_2)
        self.dw2=Thermo.dw(tmu2_1,tmu2_2)
        self.dw3=Thermo.dw(tmu4_1,tmu3_2)
        self.dw4=Thermo.dw(tmu4_1,tmu4_2)
        self.dw=Thermo.dw(tmu1,tmu2)
        
        mulinf1=Thermo.mul1(Init.Clinf1,Init.Clinf2)
        mulinf2=Thermo.mul2(Init.Clinf1,Init.Clinf2)
        musinf1=Thermo.mus1(Init.Csinf1,Init.Csinf2)
        musinf2=Thermo.mus2(Init.Csinf1,Init.Csinf2)
        
        wlinf=Thermo.wl(mulinf1, mulinf2)
        wsinf=Thermo.ws(musinf1, musinf2)
        
        tw1=Thermo.ws(tmu1_1, tmu1_2)
        tw2=Thermo.wl(tmu2_1, tmu2_2)
        tw3=Thermo.ws(tmu3_1, tmu3_2)
        tw4=Thermo.wl(tmu4_1, tmu4_2)
        
        self.tw=np.append(tw1,[tw2,tw3,tw4])
        
        self.t0s=np.linspace(minf, 0, nptxi)
        self.t0l=np.linspace(0, pinf, nptxi)
        
        self.Xi_Omega=self.calc_Xi(tw1, tw2, tw3, tw4, wsinf, wlinf)

        flinf=Thermo.fl(Init.Clinf1, Init.Clinf2)
        fsinf=Thermo.fs(Init.Csinf1, Init.Csinf2)
        
        tf1=Thermo.fs(tC1_1, tC1_2)
        tf2=Thermo.fl(tC2_1, tC2_2)
        tf3=Thermo.fs(tC3_1, tC3_2)
        tf4=Thermo.fl(tC4_1, tC4_2)
        self.Xi_F=self.calc_Xi(tf1, tf2, tf3, tf4, fsinf, flinf)


        df1=Thermo.df(tC1_1, tC1_2)
        df2=Thermo.df(tC2_1, tC2_2)
        df3=Thermo.df(tC3_1, tC3_2)
        df4=Thermo.df(tC4_1, tC4_2)
        self.df=np.append(df1,[df2,df3,df4])
        
        self.show()
        
    
        #  Xi_C1=self.calc_Xi(tC1_1, tC2_1, tC3_1, tC4_1, Init.Csinf1, Init.Clinf1)
        #  Xi_C2=self.calc_Xi(tC1_2, tC2_2, tC3_2, tC4_2, Init.Csinf2, Init.Clinf2)
        #  print("Xi_C1 = ",Xi_C1)
        #  print("Xi_C2 = ",Xi_C2)
        #  Xi_l=self.calc_Xi(self.tl1, self.tl2, self.tl3, self.tl4, 0,0)
        #  Xi_1=self.calc_Xi(1*self.tl1**0, 1*self.tl2**0, 1*self.tl3**0, 1*self.tl4**0, 0, 0)
        #  Xi_0=self.calc_Xi(0*self.tl1**0, 0*self.tl2**0, 0*self.tl3**0, 0*self.tl4**0, -1, -1)
        #  print("Xi_l = ",Xi_l,"th:",(pinf**2-minf**2)/2)
        #  print("Xi_1 = ",Xi_1,"th:",pinf-minf)
        #  print("Xi_0 = ",Xi_0,"th:",pinf-minf)

    def calc_Xi(self, v1, v2, v3, v4, l1, l4):
        return np.trapz(v1, x=self.tl1) + np.trapz(v2, x=self.tl2) + np.trapz(v3, x=self.tl3) + np.trapz(v4,x=self.tl4) - np.trapz(self.t0s**0*l1,x=self.t0s) - np.trapz(self.t0l**0*l4,x=self.t0l)
    # analytical composition profiles
    def fC1A(self,l):
        return self.CompCoeffs["A1A"]+self.CompCoeffs["B1A"]*special.erfc(-l/2/self.Diff.Ds1**0.5)
        
    def fC1B(self,l):
        return self.CompCoeffs["A1B"]+self.CompCoeffs["B1B"]*special.erfc(-l/2/self.Diff.Ds2**0.5)
        
    def fC2A(self,l):
        return self.CompCoeffs["A2A"]+self.CompCoeffs["B2A"]*special.erfc(l/2/self.Diff.Dl1**0.5)
        
    def fC2B(self,l):
        return self.CompCoeffs["A2B"]+self.CompCoeffs["B2B"]*special.erfc(l/2/self.Diff.Dl2**0.5)
    
    def fC3A(self,l):
        return self.CompCoeffs["A3A"]+self.CompCoeffs["B3A"]*special.erfc(-l/2/self.Diff.Ds1**0.5)
        
    def fC3B(self,l):
        return self.CompCoeffs["A3B"]+self.CompCoeffs["B3B"]*special.erfc(-l/2/self.Diff.Ds2**0.5)
        
    def fC4A(self,l):
        return self.CompCoeffs["A4A"]+self.CompCoeffs["B4A"]*special.erfc(l/2/self.Diff.Dl1**0.5)
        
    def fC4B(self,l):
        return self.CompCoeffs["A4B"]+self.CompCoeffs["B4B"]*special.erfc(l/2/self.Diff.Dl2**0.5)

    def fC(self, l, p, c):
        index=str(p)+c
        sign= 1 if (p%2==0) else -1
        return self.CompCoeffs["A"+index]+self.CompCoeffs["B"+index]*special.erfc(sign*l/2/self.CompCoeffs["D"+index]**0.5)

    def show(self):
        print("Printing all solution parameters")
        print("Interf 12")
        plist=[
                ("xi12", self.xi12),
                ("xi23", self.xi23),
                ("xi34", self.xi34),
                ("xi34", self.xi45),
                ("xi34", self.xi56),
                ("zeta12", self.zeta12),
                ("zeta23", self.zeta23),
                ("zeta34", self.zeta34),
                ("zeta34", self.zeta45),
                ("zeta34", self.zeta56),
                ("Xi_Omega", self.Xi_Omega),
                ("Xi_F", self.Xi_F),
            ]
        for name, val in plist:
            print(name.ljust(10," "), "=", "%7.4f" % val)












class SolverStefan6p2c:

    
    def __init__(self,Thermo, Diff, Init):
        self.Thermo=Thermo
        self.Diff=Diff
        self.Init=Init
        
        self.D1_1=Diff.Ds1
        self.D1_2=Diff.Ds2
        
        self.D2_1=Diff.Dl1
        self.D2_2=Diff.Dl2
        
        self.D3_1=Diff.Ds1
        self.D3_2=Diff.Ds2
        
        self.D4_1=Diff.Dl1
        self.D4_2=Diff.Dl2
        
        self.D5_1=Diff.Ds1
        self.D5_2=Diff.Ds2
        
        self.D6_1=Diff.Dl1
        self.D6_2=Diff.Dl2
        
        #  self.use_divide=True
        self.use_divide=False

    # intermediary functions
    # first interface
    def u112A(self,x12):
        D=self.D1_1
        l112=x12/2/(self.D1_1)**0.5
        return 2*(self.D1_1/np.pi)**0.5*np.exp(-l112**2)/special.erfc(-l112)
        
    def u112B(self,x12):
        l112=x12/2/(self.D1_2)**0.5
        return 2*(self.D1_2/np.pi)**0.5*np.exp(-l112**2)/special.erfc(-l112)
        
    def u212A(self,x12,x23):
        D=self.D2_1
        l212=x12/2/(D)**0.5
        l223=x23/2/(D)**0.5
        du212A=(special.erfc(l223)-special.erfc(l212))
        return 2*(D/np.pi)**0.5*np.exp(-l212**2)/du212A
        
    def u212B(self,x12,x23):
        D=self.D2_2
        l212=x12/2/(D)**0.5
        l223=x23/2/(D)**0.5
        du212B=(special.erfc(l223)-special.erfc(l212))
        return 2*(D/np.pi)**0.5*np.exp(-l212**2)/du212B
        

        
    ## 2nd interf
    def u223A(self,x12,x23):
        D=self.D2_1
        l212=x12/2/(D)**0.5
        l223=x23/2/(D)**0.5
        du223A=(special.erfc(l223)-special.erfc(l212))
        return -2*(D/np.pi)**0.5*np.exp(-l223**2)/du223A
        
    def u223B(self,x12,x23):
        D=self.D2_2
        l212=x12/2/(D)**0.5
        l223=x23/2/(D)**0.5
        du223B=(special.erfc(l223)-special.erfc(l212))
        return -2*(D/np.pi)**0.5*np.exp(-l223**2)/du223B
        
    def u323A(self,x23,x34):
        D=self.D3_1
        l323=x23/2/(D)**0.5
        l334=x34/2/(D)**0.5
        du323A=(special.erfc(-l323)-special.erfc(-l334))
        return -2*(D/np.pi)**0.5*np.exp(-l323**2)/du323A
        
    def u323B(self,x23,x34):
        D=self.D3_2
        l323=x23/2/(D)**0.5
        l334=x34/2/(D)**0.5
        du323B=(special.erfc(-l323)-special.erfc(-l334))
        return -2*(D/np.pi)**0.5*np.exp(-l323**2)/du323B
        
 
    ## third interf
    def u334A(self,x23,x34):
        D=self.D3_1
        l323=x23/2/(D)**0.5
        l334=x34/2/(D)**0.5
        du334A=(special.erfc(-l323)-special.erfc(-l334))
        return 2*(D/np.pi)**0.5*np.exp(-l334**2)/du334A
        
    def u334B(self,x23,x34):
        D=self.D3_2
        l323=x23/2/(D)**0.5
        l334=x34/2/(D)**0.5
        du334B=(special.erfc(-l323)-special.erfc(-l334))
        return 2*(D/np.pi)**0.5*np.exp(-l334**2)/du334B
        
    def u434A(self,x34,x45):
        D=self.D4_1
        l434=x34/2/(D)**0.5
        l445=x45/2/(D)**0.5
        du434A=(special.erfc(l445)-special.erfc(l434))
        return 2*(D/np.pi)**0.5*np.exp(-l434**2)/du434A
        
    def u434B(self,x34,x45):
        D=self.D4_2
        l434=x34/2/(D)**0.5
        l445=x45/2/(D)**0.5
        du434B=(special.erfc(l445)-special.erfc(l434))
        return 2*(D/np.pi)**0.5*np.exp(-l434**2)/du434B
        
    ## fourth interf
    def u445A(self,x34,x45):
        D=self.D4_1
        l434=x34/2/(D)**0.5
        l445=x45/2/(D)**0.5
        du445A=(special.erfc(l445)-special.erfc(l434))
        return -2*(D/np.pi)**0.5*np.exp(-l445**2)/du445A
        

        
    def u445B(self,x34,x45):
        D=self.D4_2
        l434=x34/2/(D)**0.5
        l445=x45/2/(D)**0.5
        du445B=(special.erfc(l445)-special.erfc(l434))
        return -2*(D/np.pi)**0.5*np.exp(-l445**2)/du445B
        
        
    def u545A(self,x45,x56):
        D=self.D5_1
        l545=x45/2/(D)**0.5
        l556=x56/2/(D)**0.5
        du545A=(special.erfc(-l545)-special.erfc(-l556))
        return -2*(D/np.pi)**0.5*np.exp(-l545**2)/du545A
        
    def u545B(self,x45,x56):
        D=self.D5_2
        l545=x45/2/(D)**0.5
        l556=x56/2/(D)**0.5
        du545B=(special.erfc(-l545)-special.erfc(-l556))
        return -2*(D/np.pi)**0.5*np.exp(-l545**2)/du545B
        
    ## fifth interf
    def u556A(self,x45,x56):
        D=self.D5_1
        l545=x45/2/(D)**0.5
        l556=x56/2/(D)**0.5
        du=(special.erfc(-l545)-special.erfc(-l556))
        return 2*(D/np.pi)**0.5*np.exp(-l556**2)/du
        

        
    def u556B(self,x45,x56):
        D=self.D5_2
        l545=x45/2/(D)**0.5
        l556=x56/2/(D)**0.5
        du=(special.erfc(-l545)-special.erfc(-l556))
        return 2*(D/np.pi)**0.5*np.exp(-l556**2)/du
        
        
    def u656A(self,x56):
        D=self.D6_1
        l656=x56/2/(D)**0.5
        return 2*(D/np.pi)**0.5*np.exp(-l656**2)/(special.erfc(l656))
        
    def u656B(self,x56):
        D=self.D6_2
        l656=x56/2/(D)**0.5
        return 2*(D/np.pi)**0.5*np.exp(-l656**2)/(special.erfc(l656))
    
    
    def DeltaC_1(self, zeta):
        return (1-zeta)*(self.Thermo.Cl0_1-self.Thermo.Cs0_1) 
    
    def DeltaC_2(self, zeta):
        return (zeta)*(self.Thermo.Cl0_2-self.Thermo.Cs0_2)

    
    # individual equations to solve
    def f1A(self,x12,zeta12,x23,zeta23):
        lhs  = x12 * self.DeltaC_1(zeta12)
        rhs  = self.u112A(x12)*((1-zeta12)*self.Thermo.Cs0_1-self.Init.Csinf1)
        rhs += self.u212A(x12,x23)*((1-zeta23)-(1-zeta12))*self.Thermo.Cl0_1
        return rhs-lhs
        
    def f2A(self,x12,zeta12,x23,zeta23,x34,zeta34):
        lhs  = x23 * -self.DeltaC_1(zeta23)
        rhs  = self.u223A(x12,x23)*((1-zeta23)-(1-zeta12))*self.Thermo.Cl0_1
        rhs += self.u323A(x23,x34)*((1-zeta23)-(1-zeta34))*self.Thermo.Cs0_1
        return rhs-lhs
        
    def f3A(self,x23,zeta23,x34,zeta34,x45,zeta45):
        lhs  = x34 * self.DeltaC_1(zeta34)
        rhs  = self.u334A(x23,x34)*((1-zeta34)-(1-zeta23))*self.Thermo.Cs0_1
        rhs += self.u434A(x34,x45)*((1-zeta34)-(1-zeta45))*self.Thermo.Cl0_1
        return rhs-lhs
        
    def f4A(self,x34,zeta34,x45,zeta45,x56,zeta56):
        lhs  = x45 * -self.DeltaC_1(zeta45)
        rhs  = self.u445A(x34,x45)*(zeta34-zeta45)*self.Thermo.Cl0_1
        rhs += self.u545A(x45,x56)*(zeta56-zeta45)*self.Thermo.Cs0_1
        return rhs-lhs
        
    def f5A(self,x45,zeta45,x56,zeta56):
        lhs  = x56 * self.DeltaC_1(zeta56)
        rhs  = self.u556A(x45,x56)*(zeta56-zeta45)*self.Thermo.Cs0_1 
        rhs += self.u656A(x56)*((1-zeta56)*self.Thermo.Cl0_1-self.Init.Clinf1)
        return rhs-lhs
        
    def f1B(self,x12,zeta12,x23,zeta23):
        lhs=x12 * self.DeltaC_2(zeta12)
        rhs  = self.u112B(x12)*((zeta12)*self.Thermo.Cs0_2-self.Init.Csinf2) 
        rhs += self.u212B(x12,x23)*(zeta23-zeta12)*self.Thermo.Cl0_2
        return rhs-lhs
        
    def f2B(self,x12,zeta12,x23,zeta23,x34,zeta34):
        lhs  = x23 * -self.DeltaC_2(zeta23)
        rhs  = self.u223B(x12,x23)*(zeta23-zeta12)*self.Thermo.Cl0_2
        rhs += self.u323B(x23,x34)*(zeta23-zeta34)*self.Thermo.Cs0_2
        return rhs-lhs
        
    def f3B(self,x23,zeta23,x34,zeta34,x45,zeta45):
        lhs  = x34 * self.DeltaC_2(zeta34)
        rhs  = self.u334B(x23,x34)*(zeta34-zeta23)*self.Thermo.Cs0_2
        rhs += self.u434B(x34,x45)*(zeta34-zeta45)*self.Thermo.Cl0_2
        return rhs-lhs
        
    def f4B(self,x34,zeta34,x45,zeta45,x56,zeta56):
        lhs  = x45 * -self.DeltaC_2(zeta45)
        rhs  = self.u445B(x34,x45)*(zeta45-zeta34)*self.Thermo.Cl0_2
        rhs += self.u545B(x45,x56)*(zeta45-zeta56)*self.Thermo.Cs0_2
        return rhs-lhs
        
    def f5B(self,x45,zeta45,x56,zeta56):
        lhs  = x56 * self.DeltaC_2(zeta56)
        rhs  = self.u556B(x45,x56)*(zeta45-zeta56)*self.Thermo.Cs0_2
        rhs += self.u656B(x56)*((zeta56)*self.Thermo.Cl0_2-self.Init.Clinf2)
        return rhs-lhs
        
    # transcendental function to solve for 0
    def f(self,vect):
        x12,x23,x34,x45,x56,zeta12,zeta23,zeta34,zeta45,zeta56=vect
        
        res1A=self.f1A(x12,zeta12,x23,zeta23)
        res2A=self.f2A(x12,zeta12,x23,zeta23,x34,zeta34)
        res3A=self.f3A(x23,zeta23,x34,zeta34,x45,zeta45)
        res4A=self.f4A(x34,zeta34,x45,zeta45,x56,zeta56)
        res5A=self.f5A(x45,zeta45,x56,zeta56)
        
        res1B=self.f1B(x12,zeta12,x23,zeta23)
        res2B=self.f2B(x12,zeta12,x23,zeta23,x34,zeta34)
        res3B=self.f3B(x23,zeta23,x34,zeta34,x45,zeta45)
        res4B=self.f4B(x34,zeta34,x45,zeta45,x56,zeta56)
        res5B=self.f5B(x45,zeta45,x56,zeta56)
        res=[res1A,res2A,res3A,res4A,res5A,res1B,res2B,res3B,res4B,res5B]
        return res

        

    
    def solve(self,start, maxfev):
        print("Searching xi starting from x0 = ",start)
        sol=scipy.optimize.root(self.f,start,method='hybr')
        
        print(sol.success)
        
        print("Found solution ",sol.x)

        return Stefan6p2cSolution(self.Thermo, self.Diff, self.Init, sol)

    





    

    

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
    #  x0=[-4.03936979 , 0.20633457 , 3.96190201 , 0.58974754 , 0.47092233 , 0.18618692]
    x0=[-10,-5,0,5,10,0.5,0.5,0.5,0.5,0.5]
    
    #  Init=InitSpec(Clinf1=23/scale,Clinf2=10/scale,Csinf1=3/scale,Csinf2=5/scale)
    #  x0=[-0.9,-0.9,-0.9,  0.5,  0.5,  0.5] # only valid one
    
    Solver = SolverStefan6p2c(Thermo, Diff, Init)
    sol1 = Solver.solve(x0, 1000)
    
    fig_diag, ax_diag = plt.subplots()

    ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    ax_diag.set_xlim([0,Thermo.Cl0_1])
    ax_diag.set_ylim([0,Thermo.Cl0_2])
    ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

    ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*scale) % 5 == 0])
    ax_diag.set_xlabel("$C_1$")
    ax_diag.set_ylabel("$C_2$")

    for i in range(1,sol1.Np+1):
        ax_diag.plot(sol1.tCA[i], sol1.tCB[i], linewidth=3, linestyle="-", label="Phase %d"%i)
    
    ax_diag.plot(sol1.tCA["tot"], sol1.tCB["tot"], linewidth=1, linestyle="--")

    ax_diag.set_title("$\Xi$=% 5.3f, $\Xi num$=% 5.3f"%(sol1.Xi, sol1.Xi_Omega))
    ax_diag.legend()

    plt.show()
    
    
