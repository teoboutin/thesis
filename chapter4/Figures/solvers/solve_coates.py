
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize
import scienceplots

from ThermoDilute import ThermoDilute
from solver_stefan_2p2c import Stefan2p2cSolution,SolverStefan2p2c,DiffusionSpec,InitSpec
from solver_stefan_4p2c import Stefan4p2cSolution,SolverStefan4p2c



plt.style.use('science')

fsize=plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = (fsize[0]*2,fsize[1]*2)


def plotdw(sol):
    fig_dw, ax_dw= plt.subplots()
    
    dw=sol.dw
     
    minl=min(sol.tl)
    maxl=max(sol.tl)
    ax_dw.set_xlim([minl,maxl])
    ax_dw.set_ylim([min(dw), max(dw)])
    
    ax_dw.plot(sol.tl,dw)
    
    ax_dw.plot([sol.xi, sol.xi], [min(dw), max(dw)], color="k")
    ax_dw.axvline(x=sol.xi, color="k", linewidth=1.4)
    ax_dw.axhline(y=0, color="k", linewidth=1.4)
    #  ax_dw.plot([min(sol.tl), max(sol.tl)], [0, 0], color="k")
    #  ax_dw.plot([min(sol.tl), max(sol.tl)], [sol.dwsinf,sol.dwlinf], marker="o", linewidth=0)
    alp=0.1
    ax_dw.fill_between([minl, sol.xi],[max(dw), max(dw)], color="green", alpha=alp)
    ax_dw.fill_between([sol.xi, maxl],[max(dw), max(dw)], color="red", alpha=alp)
    ax_dw.fill_between([minl, sol.xi],[min(dw), min(dw)], color="red", alpha=alp)
    ax_dw.fill_between([sol.xi, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw.set_xlabel("$\\lambda$")
    ax_dw.set_ylabel("$\\Delta^\\Phi\\omega$")
    #  ax_dw.title(r"$\Delta\omega$")
    plt.savefig("coates_test_dw.pdf")

def plotdw4p(sol):
    fig_dw, ax_dw= plt.subplots()
    
    dw=sol.dw
     
    minl=min(sol.tl)
    maxl=max(sol.tl)
    ax_dw.set_xlim([minl,maxl])
    ax_dw.set_ylim([min(dw), max(dw)])
    
    ax_dw.plot(sol.tl,dw)
    
    ax_dw.axvline(x=sol.xi12, color="k", linewidth=1.4)
    ax_dw.axvline(x=sol.xi23, color="k", linewidth=1.4)
    ax_dw.axvline(x=sol.xi34, color="k", linewidth=1.4)
    ax_dw.axhline(y=0, color="k", linewidth=1.4)
    alp=0.1
    ax_dw.fill_between([minl, sol.xi12],[max(dw), max(dw)], color="green", alpha=alp)
    ax_dw.fill_between([minl, sol.xi12],[min(dw), min(dw)], color="red", alpha=alp)
    
    ax_dw.fill_between([sol.xi12, sol.xi23],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw.fill_between([sol.xi12, sol.xi23],[max(dw), max(dw)], color="red", alpha=alp)
    
    ax_dw.fill_between([sol.xi23, sol.xi34],[max(dw), max(dw)], color="green", alpha=alp)
    ax_dw.fill_between([sol.xi23, sol.xi34],[min(dw), min(dw)], color="red", alpha=alp)

    ax_dw.fill_between([sol.xi34, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw.fill_between([sol.xi34, maxl],[max(dw), max(dw)], color="red", alpha=alp)
    
    ax_dw.set_xlabel("$\\lambda$")
    ax_dw.set_ylabel("$\\Delta^\\Phi\\omega$")
    #  ax_dw.title(r"$\Delta\omega$")
    plt.savefig("coates_test_dw_4p.pdf")



if __name__=="__main__":
    



    C0Znbeta=39.5
    C0Znalpha=34.5
    es1=0.0
    es2=0.0
    el1=np.log(1.03*(100-C0Znalpha)/(100-C0Znbeta))
    el2=np.log((100-C0Znalpha)/(100-C0Znbeta))
    el3=(C0Znalpha-C0Znbeta)
    
    
    #  # D units are resc-1 cm^2/s
    #  resc=3600 # to get results of xi in microns per h**0.5
    resc=1 # to get results of xi in microns per s**0.5
    
    DalphaZn=4*resc
    DbetaZn=25.5*resc
    
    DalphaNi=1.35e-2*resc
    DbetaNi=6.4*resc
    
    Ds2=0.025
    Ds1=DalphaNi
    
    Dl2=10
    Dl1=DbetaNi
    
    #alpha 1 beta 4
    #  Clinf1=(0.075)*100
    #  Clinf2=(1-0.459-0.075)*100
    #  Csinf1=(0.0001)*100
    #  Csinf2=(1-0.341-0.0001)*100
    #alpha 1 beta ?
    Clinf1=(0.15)*100
    Clinf2=(1-0.48-0.15)*100
    Csinf1=(0.000)*100
    Csinf2=(1-0.341-0.000)*100

    Thermo=ThermoDilute(es1,es2,0,el1,el2,el3)
    Diff=DiffusionSpec(Dl1=Dl1,Ds1=Ds1,Dl2=Dl2,Ds2=Ds2)
    Init=InitSpec(Clinf1=Clinf1,Clinf2=Clinf2,Csinf1=Csinf1,Csinf2=Csinf2)
# print the params

    #  print('DlA={:2.3}, DlB={:2.3}, DsA={:2.3}, DsB={:2.3}'.format(float(DlA),float(DlB),float(DsA),float(DsB)))
    #  print('ClinfA={:2.3}, ClinfB={:2.3}, CsinfA={:2.3}, CsinfB={:2.3}'.format(float(ClinfA),float(ClinfB),float(CsinfA),float(CsinfB)))

    starting_points_2p=[[-3.1,2,4],[-1,2,2],[0.5,1,1]]

    Solver = SolverStefan2p2c(Thermo, Diff, Init)
    #  sol1 = Solver.solve([-3.1,2,4], 1000)
    #  sol2 = Solver.solve([-1,2,2], 1000)
    #  sol3 = Solver.solve([0.5,1,1], 1000)
    
    Solver4p = SolverStefan4p2c(Thermo, Diff, Init)
    x0=[-3,-1,0.5,  0.58731941,  0.52767406,  0.12273025] # only valid one
    sol4p = Solver4p.solve(x0, 1000)

    fig_diag, ax_diag = plt.subplots()

    ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    ax_diag.set_xlim([0,max(Thermo.Cl0_1,Thermo.Cs0_1)])
    ax_diag.set_ylim([0,max(Thermo.Cl0_2,Thermo.Cs0_2)])
    ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

    ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if tick % 1 == 0])
    ax_diag.set_xlabel("$Ni$")
    ax_diag.set_ylabel("$Cu$")
    
    for x0 in starting_points_2p:
        sol1 = Solver.solve(x0, 1000)
        p1color = ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol1.xi,sol1.Xi)))[0].get_color()
    #  ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol2.xi,sol2.Xi)))
    #  ax_diag.plot(sol3.tC1, sol3.tC2, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol3.xi,sol3.Xi)))
    ax_diag.legend()
    fig_diag.savefig("coates_test_terplot.pdf")
    
    ax_diag.plot(sol4p.tC1, sol4p.tC2, linewidth=2, linestyle="--", label="4p, $\\Xi=%f$" % (sol4p.Xi))
    ax_diag.legend()
    fig_diag.savefig("coates_test_terplot_with_4p.pdf")



    fig_diagNiZn, ax_diagNiZn = plt.subplots()
    ax_diagNiZn.set_xlim([30,50])
    ax_diagNiZn.set_ylim([0,20])
    #  ax_diagNiZn.plot([100-Thermo.Cl0_2, 100-Thermo.Cl0_1], [Thermo.Cl0_1,0], color="k")
    ax_diagNiZn.axline(xy1=(100-Thermo.Cl0_2,0), xy2=(100-Thermo.Cl0_1,Thermo.Cl0_1), color="blue")
    ax_diagNiZn.axline(xy1=(100-Thermo.Cs0_2,0), xy2=(100-Thermo.Cs0_1,Thermo.Cs0_1), color="k")
    ax_diagNiZn.set_xlabel("$Zn$")
    ax_diagNiZn.set_ylabel("$Ni$")
    for x0 in starting_points_2p:
        sol1 = Solver.solve(x0, 1000)
        ax_diagNiZn.plot(100-sol1.tC1-sol1.tC2, sol1.tC1, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol1.xi,sol1.Xi)))
   
    ax_diagNiZn.plot(100-sol4p.tC1-sol4p.tC2, sol4p.tC1, linewidth=2, linestyle="--", label="4p, $\\Xi=%f$" % (sol4p.Xi))
    #  ax_diagNiZn.plot(100-sol2.tC1-sol2.tC2, sol2.tC1, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol2.xi,sol2.Xi)))
    #  ax_diagNiZn.plot(100-sol3.tC1-sol3.tC2, sol3.tC1, linewidth=2, linestyle="--", label=("$\\xi$=%f, $\\Xi=%f$" % (sol3.xi,sol3.Xi)))
    ax_diagNiZn.legend()
    ax_diagNiZn.set_yticks([tick for tick in ax_diagNiZn.get_yticks() if tick % 1 == 0])
    #  plotdw(sol1)
    #  plotdw4p(sol4p)
    plt.show()
