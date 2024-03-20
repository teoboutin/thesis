




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
cmap_RG=matplotlib.colors.ListedColormap(["red","green"], name='RG')
from scipy import special
import scipy.optimize
import scienceplots

from ThermoDilute import ThermoDilute
from solver_stefan_2p2c import Stefan2p2cSolution,SolverStefan2p2c,DiffusionSpec,InitSpec
from solver_stefan_4p2c import Stefan4p2cSolution,SolverStefan4p2c



plt.style.use(['science','grid'])

fsize=plt.rcParams["figure.figsize"]
rf=1.5
plt.rcParams["figure.figsize"] = (fsize[0]*rf,fsize[1]*rf)

#  plt.rcParams["figure.autolayout"]=True
#  plt.rcParams["savefig.bbox"]="tight"


def get_label(sol):
    return "$\\xi$=% 5.3f, $\\Xi_\Omega$=% 5.3f, $\\Xi_F$=% 5.3f" % (sol.xi,sol.Xi, sol.Xi_F)

def get_label4p(sol):
    return "$\\xi_{12}$=% 5.3f\n$\\Xi_\Omega$=% 5.3f, $\\Xi_F$=% 5.3f" % (sol.xi12,sol.Xi, sol.Xi_F)

def plotdw(sol, sol4p):
    fig_dw, ax_dw= plt.subplots()
    
    dw=sol.dw
     
    minl=min(sol.tl)
    maxl=max(sol.tl)
    ax_dw.set_xlim([minl,maxl])
    ax_dw.set_ylim([min(dw), max(dw)])
    
    
    
    ax_dw.plot(sol.tl,dw)
    
    ax_dw.plot([sol.xi, sol.xi], [min(dw), max(dw)], color="k")
    ax_dw.axvline(x=sol.xi, color="k", linewidth=1.4)
    ax_dw.set_xticks([sol.xi], minor=True)
    ax_dw.set_xticklabels(["$\\xi$=%4.2f" % sol.xi], minor=True, color='k')
    
    #  ax_dw.axvline(x=sol4p.xi12, color="grey", linewidth=1.4)
    #  ax_dw.axvline(x=sol4p.xi23, color="grey", linewidth=1.4)
    #  ax_dw.axvline(x=sol4p.xi34, color="grey", linewidth=1.4)
    
    ax_dw.axhline(y=0, color="k", linewidth=1.4)
    alp=0.1
    ax_dw.fill_between([minl, sol.xi],[max(dw), max(dw)], color="green", alpha=alp, label="Stable")
    ax_dw.fill_between([sol.xi, maxl],[max(dw), max(dw)], color="red", alpha=alp, label="Instable")
    ax_dw.fill_between([minl, sol.xi],[min(dw), min(dw)], color="red", alpha=alp)
    ax_dw.fill_between([sol.xi, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw.set_xlabel("$\\lambda$")
    ax_dw.set_ylabel("$\\Delta^\\Phi\\omega$")
    

    ax_dw.set_title("Profil du terme source $\\Delta^\\Phi\\omega$")
    ax_dw.legend()
    
    #  fig_dw.savefig("maugis_test_dw.pdf")
    fig_df, ax_df= plt.subplots()
    ax_df.plot(sol.tl,sol.df)
    ax_df.plot([sol.xi, sol.xi], [min(sol.df), max(sol.df)], color="k")
    ax_df.axvline(x=sol.xi, color="k", linewidth=1.4)
    ax_df.axvline(x=sol4p.xi12, color="grey", linewidth=1.4)
    ax_df.axvline(x=sol4p.xi23, color="grey", linewidth=1.4)
    ax_df.axvline(x=sol4p.xi34, color="grey", linewidth=1.4)
    
    ax_df.axhline(y=0, color="k", linewidth=1.4)
    ax_df.fill_between([minl, sol.xi],[max(sol.df), max(sol.df)], color="green", alpha=alp, label="Stable")
    ax_df.fill_between([sol.xi, maxl],[max(sol.df), max(sol.df)], color="red", alpha=alp, label="Unstable")
    ax_df.fill_between([minl, sol.xi],[min(sol.df), min(sol.df)], color="red", alpha=alp)
    ax_df.fill_between([sol.xi, maxl],[min(sol.df), min(sol.df)], color="green", alpha=alp)

    fig_dw4p, ax_dw4p= plt.subplots()
    
    dw=sol4p.dw
     
    minl=min(sol4p.tl)
    maxl=max(sol4p.tl)
    ax_dw4p.set_xlim([minl,maxl])
    ax_dw4p.set_ylim([min(dw), max(dw)])
    
    ax_dw4p.plot(sol4p.tl,dw)
    
    ax_dw4p.axvline(x=sol4p.xi12, color="k", linewidth=1.4)
    ax_dw4p.axvline(x=sol4p.xi23, color="k", linewidth=1.4)
    ax_dw4p.axvline(x=sol4p.xi34, color="k", linewidth=1.4)
    ax_dw4p.axhline(y=0, color="k", linewidth=1.4)
    alp=0.1
    ax_dw4p.fill_between([minl, sol4p.xi12],[max(dw), max(dw)], color="green", alpha=alp)
    ax_dw4p.fill_between([minl, sol4p.xi12],[min(dw), min(dw)], color="red", alpha=alp)
    
    ax_dw4p.fill_between([sol4p.xi12, sol4p.xi23],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw4p.fill_between([sol4p.xi12, sol4p.xi23],[max(dw), max(dw)], color="red", alpha=alp)
    
    ax_dw4p.fill_between([sol4p.xi23, sol4p.xi34],[max(dw), max(dw)], color="green", alpha=alp)
    ax_dw4p.fill_between([sol4p.xi23, sol4p.xi34],[min(dw), min(dw)], color="red", alpha=alp)

    ax_dw4p.fill_between([sol4p.xi34, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    ax_dw4p.fill_between([sol4p.xi34, maxl],[max(dw), max(dw)], color="red", alpha=alp)
    
    ax_dw4p.set_xlabel("$\\lambda$")
    ax_dw4p.set_ylabel("$\\Delta^\\Phi\\omega$")
    
    #  fig_dw4p.savefig("maugis_test_dw_4p.pdf")
    #  ax_dw.legend()
    #  ax_dw.title(r"$\Delta\omega$")
    fig_dwf, ax_df= plt.subplots()
    
    dw=sol4p.df
     
    minl=min(sol4p.tl)
    maxl=max(sol4p.tl)
    ax_df.set_xlim([minl,maxl])
    ax_df.set_ylim([min(dw), max(dw)])
    
    ax_df.plot(sol4p.tl,dw)
    
    ax_df.axvline(x=sol4p.xi12, color="k", linewidth=1.4)
    ax_df.axvline(x=sol4p.xi23, color="k", linewidth=1.4)
    ax_df.axvline(x=sol4p.xi34, color="k", linewidth=1.4)
    ax_df.axhline(y=0, color="k", linewidth=1.4)
    alp=0.1
    ax_df.fill_between([minl, sol4p.xi12],[max(dw), max(dw)], color="green", alpha=alp)
    ax_df.fill_between([minl, sol4p.xi12],[min(dw), min(dw)], color="red", alpha=alp)
    
    ax_df.fill_between([sol4p.xi12, sol4p.xi23],[min(dw), min(dw)], color="green", alpha=alp)
    ax_df.fill_between([sol4p.xi12, sol4p.xi23],[max(dw), max(dw)], color="red", alpha=alp)
    
    ax_df.fill_between([sol4p.xi23, sol4p.xi34],[max(dw), max(dw)], color="green", alpha=alp)
    ax_df.fill_between([sol4p.xi23, sol4p.xi34],[min(dw), min(dw)], color="red", alpha=alp)

    ax_df.fill_between([sol4p.xi34, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    ax_df.fill_between([sol4p.xi34, maxl],[max(dw), max(dw)], color="red", alpha=alp)
    return fig_dw,ax_dw


def find_center_between_curves(x, c1, c2):
    n=100
    dy=max(abs(c1-c2))/n
    A=0
    Qx=0
    Qy=0
    for i in range(len(x)-1):
        dx=x[i+1]-x[i]
        ry=np.linspace(c1[i], c2[i], int(abs(c1[i]-c2[i])/dy))
        dA=dx*dy
        for y in ry:
            xel=x[i]
            yel=y
            A+=dA
            Qx+=xel*dA
            Qy+=yel*dA
    return [Qx/A, Qy/A]
def phi0(l,t,xi,W):
    return 0.5*(1 + np.tanh(2* (l-xi)*(t**0.5) / W))
def h(phi):
    return phi
def p(phi):
    return phi**2 * (3.0 - 2.0 * phi)
def pprime(phi):
    return 6*phi*(1-phi)
def plotNucleationCriteria(sol, sol4p):

    #  fig_tot, ax_tot= plt.subplots()
    W=0.005
    Mphi=120
    factor=4/3*W
    
    dt=1e-8
    step=2000
    tmax=1
    base=10
    time=np.logspace(np.log(dt)/np.log(base),np.log(tmax)/np.log(base),num=step, base=base)
    
    

    
    unstable_dwltot=sol.dw * (sol.dw > 0) * (sol.tl > sol.xi)
    #  unstable_dwltot=(sol.dw - sol.dwsinf) * (sol.dw > sol.dwsinf) * (sol.tl > sol.xi)
    #  unstable_dwltot=(sol.dw  - sol.dwsinf) * (sol.dw > 0) * (sol.tl > sol.xi)
    #  unstable_dwltot=(sol.dw  - sol.dwsinf) * (sol.tl > sol4p.xi23)* (sol.tl < sol4p.xi34)
    #  unstable_dwltot=(sol.dw ) * (sol.tl > sol4p.xi23)* (sol.tl < sol4p.xi34)
    #  unstable_dwltot=(sol.dw ) * (sol.dw > 0) * (sol.tl > sol.xi) - sol.dwsinf
    #  plt.figure()
    #  plt.plot(sol.tl,unstable_dwltot)
    
    unstable_dwstot=sol.dw * (sol.dw < 0) * (sol.tl < sol.xi)
    #  unstable_dwstot=sol.dw * (sol.dw < 0) * (sol.tl < sol.xi) - sol.dwlinf
    
    unstable_dfltot=sol.df * (sol.df > 0) * (sol.tl > sol.xi)
    
    unstable_dfstot=sol.df * (sol.df < 0) * (sol.tl < sol.xi)

    max_dwl=max(unstable_dwltot)
    max_dfl=max(unstable_dfltot)
    print("max dwl", max_dwl)
    print("max dfl", max_dfl)
    excess_dwl_tot=(np.trapz(unstable_dwltot, x=sol.tl))
    excess_dws_tot=-(np.trapz(unstable_dwstot, x=sol.tl))
    excess_dfl_tot=(np.trapz(unstable_dfltot, x=sol.tl))
    excess_dfs_tot=-(np.trapz(unstable_dfstot, x=sol.tl))
    
    t_ex_dwl_tot=excess_dwl_tot* time**0.5
    t_ex_dws_tot=excess_dws_tot* time**0.5
    t_ex_dfl_tot=excess_dfl_tot* time**0.5
    t_ex_dfs_tot=excess_dfs_tot* time**0.5
    
    t_lstar_1=factor/(-(sol4p.Xi-sol.Xi)*time**0.5)
    t_lstar_2=factor/t_ex_dwl_tot
    t_lstar_3=np.array([])
    
    for t in time:

        pp=pprime(phi0(sol.tl,t,sol.xi,W))
        
        unstable_dwl_interf=unstable_dwltot * pp
        unstable_dws_interf=unstable_dwstot * pp
        
        unstable_dfl_interf=unstable_dfltot * pp
        
        excess_dwl=abs(np.trapz(unstable_dwl_interf, x=sol.tl))*t**0.5
        excess_dws=abs(np.trapz(unstable_dws_interf, x=sol.tl))*t**0.5
        
        if excess_dwl>0:
            t_lstar_3=np.append(t_lstar_3,[factor/excess_dwl])
        #  t_ex_dws=np.append(t_ex_dws,[excess_dws])
        

        
        

    
    fig, ax= plt.subplots()
    ax.plot(time, t_lstar_1, label="Condition 1")
    ax.plot(time, t_lstar_2, label="Condition 2")
    ax.plot(time[0:len(t_lstar_3)], t_lstar_3, label="Condition 3")
    #  ax.plot(time, t_lstar_2, label="Condition based on total energy gain")
    #  ax.plot(time, (4/3 * W ) / (t_ex_dwl_tot), label="Condition based on total energy gain")
    #  ax.plot(time, factor / t_ex_dwl, label="Condition based on total energy gain")
    #  ax.plot(time, factor / t_ex_dfl_tot, label="Condition based on total energy gain")
    #  data=np.array([[time[i], t_ex_dwl[i]] for i in range(len(time))])
 
    #  ax.plot(time[:len(t_lambda_correct)], t_lambda_correct, label="Condition corrected for energy inside diffuse interface")
    #  ax.plot(time, factor / t_ex_dws4p, label="Condition with comp of total energy")
    #  ax.plot(t_cond_liq_nucl_time, factor / t_ex_dwl_tot, label="test")
    #  ax.plot(t_cond_liq_nucl_time_correct, factor / t_ex_dwl, label="test")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\lambda^*(t)$")
    
    ax.set_ylim([1e-2, 1e7])
    
    liquid_amin=np.argmin( t_lstar_3)
    #  liquid_amin=np.argmin( -(4/3 * W + t_W4dW2) / t_ex_dwl)
    liquid_lmin=(t_lstar_3)[liquid_amin]
    liquid_tmin=time[liquid_amin]
    print(liquid_amin)
    print(liquid_lmin)
    print(liquid_tmin)
    
    ax.set_yticks([liquid_lmin], minor=True)
    ax.set_yticklabels([int(liquid_lmin)], minor=True, color='r')
    ax.axhline(y=liquid_lmin, color='r', linestyle="--")
    ax.scatter(liquid_tmin,liquid_lmin, color='r', marker="o", zorder=10)
    ax.set_title("Minimum value of $\lambda$ needed for the nucleation of solid in the liquid phase, \nas a function of time.")
    
    ax.scatter([1.25e-3,1.8e-5],[400,500], marker="+", label="Points of nucleation")
    ax.scatter([5e-3],[300], marker="+", label="Points with no nucluation")
    ax.legend()
    
    
    #  fig.savefig("maugis_test_nucleation_condition_pondered.pdf")
    #  fig, ax= plt.subplots()
    #  ax.plot(time, factor / t_ex_dws_tot, label="Condition based on total energy gain")
    #  ax.plot(time, factor / t_ex_dws, label="Condition corrected for energy inside diffuse interface")
    #  ax.set_xscale("log")
    #  ax.set_yscale("log")
    
    #  ax.set_xlabel("$t$")
    #  ax.set_ylabel("$\lambda^*(t)$")
    #  ax.set_yticks([min(factor / t_ex_dws)], minor=True)
    #  ax.set_yticklabels([int(min(factor / t_ex_dws))], minor=True, color='r')
    #  ax.axhline(y=min(factor / t_ex_dws), color='r', linestyle="--")
    #  ax.scatter(time[np.argmin(factor / t_ex_dws)],min(factor / t_ex_dws), color='r', marker="o", zorder=10)
    #  ax.set_title("Minimum value of $\lambda$ needed for the nucleation of liquid in the solid phase, \nas a function of time.")
    
    #  ax.scatter([0.000012],[500], marker="+", label="Points of nucleation")
    #  ax.scatter([5e-3],[300], marker="+", label="Points with no nucluation")
    #  ax.legend()
    
    
    #  print("lambda condition = ", factor/max( t_ex_dwl))
    
    #  def get_nucleation_time(plambda):
        #  print(plambda, next(time[idx] for idx, value in enumerate((factor / t_ex_dws)) if value < plambda))
        #  next(time[idx] for idx, value in enumerate(factor / t_ex_dws) if value < plambda)
    #  get_nucleation_time(300)
    #  get_nucleation_time(400)
    #  get_nucleation_time(500)
    return fig, ax
        




if __name__=="__main__":
    

    
# input parameters
    def Si_solubility(TK): # from Fournier and Rowe 1977
        return 10**(-810/TK+4.82)
    TC=90
    
    TK=273.15+TC
    
    Cl0_1=Si_solubility(TK) * 1e-6 

    #  fig,ax=plt.subplots()
    #  ax.plot(1000/np.linspace(300,700,100), Si_solubility(np.linspace(300,700,100)))
    #  ax.set_yscale('log')
    #  ax.invert_xaxis()
    scale=100
    kB=1.38*10**-23
    R=8.314
    def Si_autodiff(TK): # from Fournier and Rowe 1977
        return (328*np.exp(-579000 /(R*TK)))

    print(Si_autodiff(1110+273))
    print(Si_autodiff(TK))
    #  Cl0_1=400e-6 # at 100°C, 
    
    Cl0_1=400 * 1e-6 
    Cl0_2=0.086
    Cs0_1=0.9
    
    Cs0_2=Cl0_2-(Cl0_1-Cs0_1)
    
    el1=-np.log(Cl0_1/Cs0_1)
    el2=-np.log(Cl0_2/Cs0_2)
    el3=Cl0_1-Cs0_1
    
    es1=0
    es2=0
    es3=0

    Thermo=ThermoDilute(es1,es2,es3,el1,el2,el3)
    
    el1=9.7848
    el2=6.1948
    es1=0
    es2=0
    el3=-0.993
    es3=0.0
    Thermo3=ThermoDilute(es1,es2,es3,el1,el2,el3)
    
    Dm=1
    Dl=1e8
    Ds1=1e-2
    Ds2=1.0
    Diff=DiffusionSpec(Dl1=Dl,Ds1=Ds1,Dl2=Dl,Ds2=Ds2)
    Dl=1e2
    Diff2=DiffusionSpec(Dl1=Dl,Ds1=Ds1,Dl2=Dl,Ds2=Ds2)
    
    Diff3=DiffusionSpec(Dl1=100,Ds1=0.01,Dl2=100,Ds2=0.1)
    
    L_undersat=0.4
    L_mix=0.9
    
    #  L_undersat=0.0
    #  L_mix=0.9
    Cl0init1=(1-L_undersat) * Thermo.Cl0_1
    Cl0init2=(1-L_undersat) * Thermo.Cl0_2
    
    Clinf1=(  L_mix) * Cl0init1
    Clinf2=(1-L_mix) * Cl0init2
    
    #  Clinf1=150e-6
    #  Clinf2=100e-6
    Init=InitSpec(Clinf1=Clinf1,Clinf2=Clinf2,Csinf1=0.7,Csinf2=0.3)
    
    
    Init3=InitSpec(Clinf1=0.00002,Clinf2=0.0003,Csinf1=0.7,Csinf2=0.3)
        
    # print the params
    print(Diff)
    print(Init)

    Solver = SolverStefan2p2c(Thermo, Diff, Init)
    Solver2 = SolverStefan2p2c(Thermo, Diff2, Init)
    Solver3 = SolverStefan2p2c(Thermo3, Diff3, Init3)
    Solver34 = SolverStefan4p2c(Thermo3, Diff3, Init3)
    
    
    sol1 = Solver.solve([-4.7*Dm**0.5,-0.5,-1], 1000)
    sol2 = Solver2.solve([-4.7*Dm**0.5,-0.5,-1], 1000)
    sol3 = Solver.solve([311*Dm**0.5,-0.5,-1], 1000)
    
    sol23 = Solver3.solve([100*Dm**0.5,-0.5,-1], 1000)
    sol34 = Solver34.solve([-1,0,1,0.6,0.5,0.4], 1000)
    
    print("sol1")
    sol1.show()
    print("sol2")
    sol2.show()
    print("sol3")
    sol3.show()
    print("sol23")
    sol23.show()
    sol34.show()
    
    print("4p")
    
    Solver4p = SolverStefan4p2c(Thermo, Diff, Init)
    
    x0=[-2.91461196,-0.35893441,7.21155476, 0.58294038, 0.52270732, 0.12272463] # only valid one
    sol4p = Solver4p.solve(x0, 1000)


    fig_diag, ax_diag = plt.subplots(figsize=(6,6))

    
    ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    ax_diag.plot([1, 0], [0,1], color="k")
    ax_diag.plot([0, 0], [0,1], color="k")
    ax_diag.plot([0, 1], [0,0], color="k")
    ax_diag.set_aspect('equal')
    ax_diag.set_xlim([-0.1,1.1])
    #  ax_diag.set_ylim([-0.1,1.1])
    ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

    ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*scale) % 5 == 0])
    ax_diag.set_xlabel("Composition $C_1$ (Silicium)")
    ax_diag.set_ylabel("Composition $C_2$ (Bore)")
    
    ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="-", label="$D^l$=%3.1e"%sol1.Diff.Dl1 , marker="")
    ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="-", label="$D^l$=%3.1e"%sol2.Diff.Dl1 , marker="")
    ax_diag.legend()


    #  fig_diag.savefig("maugis_test_terplot_with_4p.pdf")
    
    # inset Axes....
    x1, x2, y1, y2 = 0,Thermo.Cl0_1,0,Thermo.Cl0_2  # subregion of the original image
    axins = ax_diag.inset_axes(
        [0.1, 0.1, 0.3, 0.3], # xcorner_sw, ycorner_sw, xsize, ysize
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[], xticks=[], yticks=[])
    axins.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        right=False,         # ticks along the top edge are off
        left=False,         # ticks along the top edge are off
        labelleft=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
        
    axins.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    axins.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    axins.plot([1, 0], [0,1], color="k")
    axins.plot([0, 0], [0,1], color="k")
    axins.plot([0, 1], [0,0], color="k")
    axins.set_xlim([0,max(Thermo.Cl0_1,Thermo.Cl0_1)])
    axins.set_ylim([0,max(Thermo.Cl0_2,Thermo.Cl0_2)])
    axins.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

    axins.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="-")
    axins.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="-")
    #  axins.plot(sol4p.tC1, sol4p.tC2, linewidth=2, linestyle="--", label="Solution 4p" )

    #  axins.set_aspect('equal')
    ax_diag.indicate_inset_zoom(axins)
    
    
    fig_diag.savefig("glass_diag_instability.pdf")
    #  fig_dw,ax_dw=plotdw(sol1, sol4p)
    #  fig_crit,ax_crit=plotNucleationCriteria(sol1, sol4p)

    
    #  fig_dw,ax_dw=plotdw4p(sol4p)
    


    
    ## plot all dw with same ref
    #  fig_dw, ax_dw= plt.subplots()
    
    #  dw=np.append(sol1.dw, [sol2.dw, sol3.dw])
     
    #  all_tl=np.append(sol1.tl-sol1.xi, [sol2.tl-sol2.xi, sol3.tl-sol3.xi])
    #  minl=min(all_tl)
    #  maxl=max(all_tl)
    #  ax_dw.set_xlim([minl,maxl])
    #  ax_dw.set_ylim([min(dw), max(dw)])
    
    
    #  r=1
    #  ax_dw.plot(sol1.tl-r*sol1.xi,sol1.dw, label="Solution 1")

    
    #  ax_dw.axvline(x=0, color="k", linewidth=1.4)
    #  ax_dw.set_xticks([0], minor=True)
    #  ax_dw.set_xticklabels(["$\\xi$"], minor=True, color='k')
    

    
    #  ax_dw.axhline(y=0, color="k", linewidth=1.4)
    #  alp=0.1
    #  ax_dw.fill_between([minl, 0],[max(dw), max(dw)], color="green", alpha=alp, label="Domaine stable")
    #  ax_dw.fill_between([0, maxl],[max(dw), max(dw)], color="red", alpha=alp, label="Domaine instable")
    #  ax_dw.fill_between([minl, 0],[min(dw), min(dw)], color="red", alpha=alp)
    #  ax_dw.fill_between([0, maxl],[min(dw), min(dw)], color="green", alpha=alp)
    #  ax_dw.set_xlabel("$\\lambda-\\xi$")
    #  ax_dw.set_ylabel("$\\Delta^\\Phi\\omega$")
    

    #  ax_dw.set_title("Profil du terme source $\\Delta^\\Phi\\omega$")
    #  ax_dw.legend()
    #  fig_dw.savefig("maugis_test_dw_2p_all.pdf")
    plt.show()
