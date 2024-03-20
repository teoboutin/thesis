

import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import scipy.optimize
import ternary
import configparser
import random
tab20=plt.get_cmap("tab20")



class ThermoTernaryQuadratic:
	def __init__(self,KAA,KAB,KBB,mA,mB,Q):
		self.KAA=KAA
		self.KAB=KAB
		self.KBB=KBB
		self.mA=mA
		self.mB=mB
		self.Q=Q
		self.K=KAA*KBB-KAB**2
		assert(self.K>=0)
	
	def f(self,cA,cB):
		return 0.5*self.KAA*(cA-self.mA)**2+self.KAB*(cB-self.mB)*(cA-self.mA)+0.5*self.KBB*(cB-self.mB)**2+self.Q
	def muA(self,cA,cB):
		return self.KAA*(cA-self.mA)+self.KAB*(cB-self.mB)
	def muB(self,cA,cB):
		return self.KAB*(cA-self.mA)+self.KBB*(cB-self.mB)
	def cA(self, muA,muB):
		return self.mA+ (self.KBB*muA-self.KAB*muB)/self.K
	def cB(self, muA,muB):
		return self.mB+ (self.KAA*muB-self.KAB*muA)/self.K
	def w(self,muA,muB):
		return self.f(self.cA(muA,muB),self.cB(muA,muB)) - muA*self.cA(muA,muB) - muB*self.cB(muA,muB)
		
	def print_params(self):
		print("mA=%6.3f, Q=%6.3f" % (self.mA,self.Q))
		
		
class ThermoTernaryDilute:
	def __init__(self,eA,eB,eC):

		self.eA=eA
		self.eB=eB
		self.eC=eC

	def f(self,cA, cB):
		return cA*(self.eA-1+np.log(cA)) + cB*(self.eB-1+np.log(cB)) +self.eC
	def muA(self,cA,cB):
		return self.eA+np.log(cA)
	def muB(self,cA,cB):
		return self.eB+np.log(cB)
	def cA(self, muA, muB):
		return np.exp(muA-self.eA)
	def cB(self, muA, muB):
		return np.exp(muB-self.eB)
		
	def w(self,muA, muB):
		return self.f(self.cA(muA, muB), self.cB(muA, muB)) - muA*self.cA(muA, muB) - muB*self.cB(muA, muB)
		
	def print_params(self):
		print("eA=%6.3f, eB=%6.3f, eC=%6.3f" % (self.eA,self.eB,self.eC))





class Solver2Phase1Int:
	
	def __init__(self,ThermoP0,ThermoP1,DA0,DB0,DA1,DB1,C0infA,C0infB,C1infA,C1infB):
		self.ThermoP0=ThermoP0
		self.ThermoP1=ThermoP1
		self.DA0=DA0
		self.DB0=DB0
		
		self.DB1=DB1
		self.DA1=DA1
		
		self.C0infA=C0infA
		self.C1infA=C1infA
		self.C0infB=C0infB
		self.C1infB=C1infB
		
		if (isinstance(self.ThermoP0,ThermoTernaryDilute) and isinstance(self.ThermoP1,ThermoTernaryDilute)):
			
			kA=np.exp(self.ThermoP0.eA-self.ThermoP1.eA)
			kB=np.exp(self.ThermoP0.eB-self.ThermoP1.eB)
			
			deC=self.ThermoP1.eC-self.ThermoP0.eC
			
			self.Cl0A=kA*deC/(kA-1)
			self.Cs0A=deC/(kA-1)
			self.Cl0B=kB*deC/(kB-1)
			self.Cs0B=deC/(kB-1)
		



	def print_params(self):
		self.ThermoP0.print_params()
		self.ThermoP1.print_params()

		
		print("DA0=%6.3f | DB0=%6.3f | DA1=%6.3f | DB1=%6.3f" % (self.DA0,self.DB0,self.DA1,self.DB1))
		# ~ print("C0inf=%6.3f | C2inf=%6.3f" % (self.C0inf,self.C2inf))
		
		if (isinstance(self.ThermoP0,ThermoTernaryDilute) and isinstance(self.ThermoP1,ThermoTernaryDilute)):
			
			
			print("Cl0A=%6.3f | Cs0A=%6.3f | Cl0B=%6.3f | Cs0B=%6.3f" % (self.Cl0A,self.Cs0A,self.Cl0B,self.Cs0B))
			
	def print_AS(self):

		keylist=["A1A","B1A","A0A","B0A","A1B","B1B","A0B","B0B","xi"]
		for key in keylist:
			print(key+"="+str(self.AnalyticalSolution[key]))

		
		
	# analytical composition profiles
	def fC0A(self,l):
		return self.AnalyticalSolution["A0A"]+self.AnalyticalSolution["B0A"]*special.erfc(-l/2/self.DA0**0.5)
	def fC0B(self,l):
		return self.AnalyticalSolution["A0B"]+self.AnalyticalSolution["B0B"]*special.erfc(-l/2/self.DB0**0.5)
	def fC1A(self,l):
		return self.AnalyticalSolution["A1A"]+self.AnalyticalSolution["B1A"]*special.erfc(l/2/self.DA1**0.5)
	def fC1B(self,l):
		return self.AnalyticalSolution["A1B"]+self.AnalyticalSolution["B1B"]*special.erfc(l/2/self.DB1**0.5)

	

	
	
	# intermediary functions
	def u(self,x):
		return (2/np.pi**0.5)*np.exp(-x**2)/special.erfc(x)

	def u0A(self,x):
		return (self.DA0)**0.5 * self.u(-x/2/self.DA0**0.5)
	def u0B(self,x):
		return (self.DB0)**0.5 * self.u(-x/2/self.DB0**0.5)
		
	def u1A(self,x):
		return (self.DA1)**0.5 * self.u(x/2/self.DA1**0.5)
	def u1B(self,x):
		return (self.DB1)**0.5 * self.u(x/2/self.DB1**0.5)


	

	
	# individual equations to solve
	def f1(self,xi, muA,muB):
		lhs=xi*(self.ThermoP1.cA(muA,muB)-self.ThermoP0.cA(muA,muB))
		rhs=self.u1A(xi)*(self.ThermoP1.cA(muA,muB)-self.C1infA) + self.u0A(xi)*(self.ThermoP0.cA(muA,muB)-self.C0infA)
		return rhs-lhs
		
	def f2(self,xi, muA,muB):
		lhs=xi*(self.ThermoP1.cB(muA,muB)-self.ThermoP0.cB(muA,muB))
		rhs=self.u1B(xi)*(self.ThermoP1.cB(muA,muB)-self.C1infB) + self.u0B(xi)*(self.ThermoP0.cB(muA,muB)-self.C0infB)
		return rhs-lhs
		
	def f3(self,muA,muB):
		lhs=self.ThermoP0.w(muA,muB)
		rhs=self.ThermoP1.w(muA,muB)
		return rhs-lhs
		

	
	
	# transcendental function to solve for 0
	def f(self,vect):
		xi,muA,muB=vect
		
		res1=self.f1(xi, muA,muB)
		res2=self.f2(xi, muA,muB)
		res3=self.f3(muA,muB)
	
		res=[res1,res2,res3]
		return res





	def compute_all_params(self,sol):
		xi, muA,muB=sol.x
		
		# ~ print("xi   :",xi)
		# ~ print("muA  :",muA)
		# ~ print("muB  :",muB)
		
		AnalyticalSolution=dict()
		
		AnalyticalSolution["xi"]=xi
		
		AnalyticalSolution["muA"]=muA
		AnalyticalSolution["muB"]=muB
		
		# ~ AnalyticalSolution["XI"]=-dfC*(x12+x34-x23)
		# ~ print('gp vvar:', 	AnalyticalSolution["XI"])
		
		AnalyticalSolution["c0A"]=self.ThermoP0.cA(muA,muB)
		AnalyticalSolution["c0B"]=self.ThermoP0.cB(muA,muB)
		
		AnalyticalSolution["c1A"]=self.ThermoP1.cA(muA,muB)
		AnalyticalSolution["c1B"]=self.ThermoP1.cB(muA,muB)
		
	
		
	
		# first phase
		AnalyticalSolution["A0A"]=self.C0infA
		AnalyticalSolution["A0B"]=self.C0infB
		AnalyticalSolution["B0A"]=(AnalyticalSolution["c0A"]-self.C0infA)/special.erfc(-xi/2/self.DA0**0.5)
		AnalyticalSolution["B0B"]=(AnalyticalSolution["c0B"]-self.C0infB)/special.erfc(-xi/2/self.DB0**0.5)
	
		# second phase
		AnalyticalSolution["A1A"]=self.C1infA
		AnalyticalSolution["A1B"]=self.C1infB
		AnalyticalSolution["B1A"]=(AnalyticalSolution["c1A"]-self.C1infA)/special.erfc(xi/2/self.DA1**0.5)
		AnalyticalSolution["B1B"]=(AnalyticalSolution["c1B"]-self.C1infB)/special.erfc(xi/2/self.DB1**0.5)

	
		self.AnalyticalSolution=AnalyticalSolution
		
	def search_all_xi(self,start, maxfev, ntry):
		
		self.solution_list=[]
		for i in range(ntry):
		

			# ~ print("Searching xi starting from x0 = ",start)
			sol=scipy.optimize.root(self.f,start,method='hybr')
			success=sol.success

			
			
			if success:
				# ~ print(sol.message)
				isnew=True
				for prev_sol in self.solution_list:
					isnew= isnew and not( np.linalg.norm(sol.x-prev_sol.x) < 1e-4)
				
				xi, muA,muB=sol.x
				
				c0A=self.ThermoP0.cA(muA,muB)
				c0B=self.ThermoP0.cB(muA,muB)
				c1A=self.ThermoP1.cA(muA,muB)
				c1B=self.ThermoP1.cB(muA,muB)
				
				isvalid= (c0A+c0B)<1 and (c1A+c1B)<1 and c0A>0 and c0B>0 and c1A>0 and c1B>0
				if isnew and isvalid:
					print("xi = %6.3f | muA = %6.3f | muB = %6.3f" % (xi,muA,muB))
					self.solution_list.append(sol)
					
				
			start=(np.random.normal(0,1), np.random.normal(0,1),np.random.normal(0,1))
				
		self.solution_list.sort(key=lambda sol:sol.x[0])
		return sol
	
	
	def plot_profile_comp(self,t,x01ini,x12ini, minf=0, pinf=0):
		if minf==0:
			minf=xi01-5*self.D0**0.5
		if pinf==0:
			pinf=xi12+5*self.D2**0.5
		xi01=self.AnalyticalSolution["xi01"]
		xi12=self.AnalyticalSolution["xi12"]
		fig, ax = plt.subplots()
		
		# plot first phase solution
		nptxi=1000
		tl0=np.linspace(minf-x01ini/t**0.5, xi01, nptxi)
		tC0=self.fC0(tl0)
		# ~ tmu=ThermoP0.mu(tC0)
		# ~ tdw=ThermoP0.dw(tmuA,tmuB)
		# ~ dw=tdw
		ax.plot(tl0+x01ini/t**0.5,tC0, color="red", linewidth=2)
		# plot second phase sol
		nptxi=1000
		tl1=np.linspace(xi01, xi12, nptxi)
		tC1=self.fC1(tl1)
		# ~ tmuA=thermo.mu1A(tC2A,tC2B)
		# ~ tmuB=thermo.mu1B(tC2A,tC2B)
		# ~ tdw=thermo.dw(tmuA,tmuB)
		# ~ dw=np.append(dw,tdw)
		tl1r= xi01 + x01ini/t**0.5 + (tl1 - tl1[0])/(tl1[-1]-tl1[0]) * ((x12ini/t**0.5+xi12) - (x01ini/t**0.5+xi01))
		ax.plot(tl1r,tC1, color="red", linewidth=2)
	
	
		# plot fourth phase sol
		nptxi=1000
		tl2=np.linspace(xi12, pinf-x12ini/t**0.5, nptxi)
		tC2=self.fC2(tl2)
		# ~ tmuA=thermo.mu1A(tC4A,tC4B)
		# ~ tmuB=thermo.mu1B(tC4A,tC4B)
		# ~ tdw=thermo.dw(tmuA,tmuB)
		# ~ dw=np.append(dw,tdw)
		ax.plot(tl2+x12ini/t**0.5,tC2, color="red", linewidth=2)
		ax.axvline(x=xi01+x01ini/t**0.5, color="black")
		ax.axvline(x=xi12+x12ini/t**0.5, color="black")
		
		ax.set_xlim([minf,pinf])
		
		return fig,ax
		
	def plot_tax(self,t, minf=0, pinf=0, W=0):
		xi=self.AnalyticalSolution["xi"]
		muA=self.AnalyticalSolution["muA"]
		muB=self.AnalyticalSolution["muB"]
				
				
		if minf==0:
			minf=xi-5*self.DA0**0.5
		if pinf==0:
			pinf=xi+5*self.DA1**0.5

		
		
		scale=1
		# create tern axes, and add all beautiful things
		figure, tax = ternary.figure(scale=scale)
		tax.boundary(linewidth=2.0)
		tax.gridlines(multiple=scale/10, color="k")
		tax.set_axis_limits({'b': [0, scale], 'l': [0, scale], 'r': [0, scale]})
		# get and set the custom ticks:
		tax.get_ticks_from_axis_limits(multiple=scale/10)
		tick_formats = "%.1f"
		tax.set_custom_ticks(fontsize=10, offset=0.02,multiple=scale/10, tick_formats=tick_formats)
		
		tax.ax.axis("off")
		tax.ax.set_aspect('equal', adjustable='box')
		tax._redraw_labels()
		fontsize=14
		tax.left_corner_label("C", fontsize=fontsize, offset=0.17, color='k')
		tax.right_corner_label("A", fontsize=fontsize, offset=0.07, color='k')
		tax.top_corner_label("B", fontsize=fontsize, offset=0.2, color='k')

		figure.suptitle("xi = %6.3f | muA = %6.3f | muB = %6.3f" % (xi,muA,muB))
		# plot first phase solution
		nptxi=1000
		tl0=np.linspace(minf, xi, nptxi)
		tC0A=self.fC0A(tl0)
		tC0B=self.fC0B(tl0)
		tC0ter = [ (tC0A[i], tC0B[i], 1-tC0A[i]- tC0B[i] ) for  i in range(nptxi)]

		tax.plot(tC0ter)
		
		if (isinstance(self.ThermoP0,ThermoTernaryDilute) and isinstance(self.ThermoP1,ThermoTernaryDilute)):
			Frontier=[(self.Cs0A, 0, 1-self.Cs0A), (0, self.Cs0B, 1-self.Cs0B)]
			tax.plot(Frontier, color='green')
			
		
		Bound=[(self.C0infA, self.C0infB)]
		tax.scatter(Bound, color='green')
		
		# plot second phase sol
		nptxi=1000
		tl1=np.linspace(xi, pinf, nptxi)
		tC1A=self.fC1A(tl1)
		tC1B=self.fC1B(tl1)
		tC1ter = [ (tC1A[i], tC1B[i], 1-tC1A[i]- tC1B[i] ) for  i in range(nptxi)]
		tax.plot(tC1ter)

		if (isinstance(self.ThermoP0,ThermoTernaryDilute) and isinstance(self.ThermoP1,ThermoTernaryDilute)):
			Frontier=[(self.Cl0A, 0, 1-self.Cl0A), (0, self.Cl0B, 1-self.Cl0B)]
			tax.plot(Frontier, color='blue')
		
		Bound=[(self.C1infA, self.C1infB)]
		tax.scatter(Bound, color='blue')
	
		
		
		if W>0:
			# plot interpolated profile
			
			tl=np.append(tl0,tl1)
			
			tmu0A=self.ThermoP0.muA(tC0A, tC0B)
			tmu0B=self.ThermoP0.muB(tC0A, tC0B)
			tmu1A=self.ThermoP1.muA(tC1A, tC1B)
			tmu1B=self.ThermoP1.muB(tC1A, tC1B)
			
			tmuA=np.append(tmu0A,tmu1A)		
			tmuB=np.append(tmu0B,tmu1B)
			
			
			
			
			tphi= self.phi0(tl,t,xi,W)
			
			tcA= tphi * self.ThermoP1.cA(tmuA, tmuB) + (1-tphi) * self.ThermoP0.cA(tmuA, tmuB)
			tcB= tphi * self.ThermoP1.cB(tmuA, tmuB) + (1-tphi) * self.ThermoP0.cB(tmuA, tmuB)
			
			tcter= [(tcA[i], tcB[i]) for i in range(len(tl))]
			# ~ tax.plot(tcter, color="purple", marker='')
			
			
			tc1A=self.fC1A(tl)
			tc1B=self.fC1B(tl)
			tc0A=self.fC0A(tl)
			tc0B=self.fC0B(tl)
			tcA= tphi * tc1A + (1-tphi) * tc0A
			tcB= tphi * tc1B + (1-tphi) * tc0B
			
			tcter= [(tcA[i], tcB[i]) for i in range(len(tl))]
			tax.plot(tcter, color="red", marker='')
			
			
			# ~ fig=plt.figure()
			# ~ plt.plot(tl,tphi)
			# ~ plt.plot(tl,self.ThermoP0.cA(tmuA, tmuB))
			# ~ plt.plot(tl,self.ThermoP1.cA(tmuA, tmuB))
			
			
			
			
		
		return figure,tax
	
	def phi0(self,l,t,xi,W):
		return 0.5*(1 + np.tanh(2* (l-xi)*(t**0.5) / W))

if __name__ == "__main__":
	elA=0.0
	elB=0.0
	esA=elA+np.log(5/4)
	esB=elB+np.log(4/3)
	elC=0.1
	esC=0.0
	


	D0=20
	DA1=36/D0
	DA0=16/D0
	DB1=4/D0
	DB0=1/D0




	ThermoP1=ThermoTernaryDilute(elA,elB,elC)
	ThermoP0=ThermoTernaryDilute(esA,esB,esC)

	# ~ DA1=1.8
	# ~ DB1=0.2
	# ~ DA0=0.8
	# ~ DB0=0.06
	
		
	C1infA=23/50
	C1infB=3/50
	C0infA=3/50
	C0infB=12/50
	
	# ~ C0infA=0.06
	# ~ C0infB=0.24
	
	# ~ C1infA=0.46
	# ~ C1infB=0.06
	
	Solver=Solver2Phase1Int(ThermoP0,ThermoP1,DA0,DB0,DA1,DB1,C0infA,C0infB,C1infA,C1infB)
	
	ntry=10

	x0=[-1,0,0]
	
	sol=Solver.search_all_xi(x0, 1000,ntry)
	
	for sol in Solver.solution_list:
		Solver.compute_all_params(sol)
			
			
		
		Solver.plot_tax(1,W=0.1)
	
	Solver.print_params()
	
	plt.show()

	
