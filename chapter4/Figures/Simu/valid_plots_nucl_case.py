import math
import fnmatch
import os
import datetime
import subprocess
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.scale
import numpy as np
import configparser
import AnalyticalSolverTernary
import time
tab20=plt.get_cmap("tab20")
tab10=plt.get_cmap("tab10")
from solvers.ThermoDilute import ThermoDilute
from solvers.solver_stefan_2p2c import Stefan2p2cSolution,SolverStefan2p2c,DiffusionSpec,InitSpec
from solvers.solver_stefan_4p2c import Stefan4p2cSolution,SolverStefan4p2c
from solvers.envelope import envelope3 as envelope
from solvers.envelope import envelope_monotonous

import scienceplots
plt.style.use('science')

fsize=plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = (fsize[0]*2,fsize[1]*2)


colors=plt.rcParams['axes.prop_cycle'].by_key()['color']

script_start_time=time.time()

#  dir="test"
dir="case_insta_l300_metastable"
dir="case_insta_l400_slow_nucleation"
dir="case_insta_l500_nucleation"


os.chdir(dir)
    
# ~ print(datetime.datetime.fromtimestamp(os.path.getmtime("integrals_0.csv")))
# ~ print(os.path.getmtime("integrals_1.csv"))

if (os.path.exists("integrals_0.csv") and os.path.getmtime("integrals_0.csv")<os.path.getmtime("Simu_insta_maugis.xmf")):
    subprocess.call("rm integrals*", shell=True)
    subprocess.call("rm profile*", shell=True)
    
if (not(os.path.exists("integrals_0.csv"))):
    print("rebuilding output files from paraview")
    subprocess.call("$HOME/Apps/ParaView-5.10.1-MPI-Linux-Python3.9-x86_64/bin/pvpython ../integrate.pvpy", shell=True)



# No need to touch after this

if (not(os.path.exists("Figs"))):
    os.mkdir("Figs")

#  subprocess.call("rm Figs/*", shell=True)


cfg=configparser.ConfigParser()
filename="settings.ini"
cfg.read(filename)


elA=cfg.getfloat("params","elA")
elB=cfg.getfloat("params","elB")
elC=cfg.getfloat("params","elC")
esA=cfg.getfloat("params","esA")
esB=cfg.getfloat("params","esB")
esC=cfg.getfloat("params","esC")


ThermoP0=AnalyticalSolverTernary.ThermoTernaryDilute(esA,esB,esC)
ThermoP1=AnalyticalSolverTernary.ThermoTernaryDilute(elA,elB,elC)


DA0=cfg.getfloat("params","DA0")
DA1=cfg.getfloat("params","DA1")
DB0=cfg.getfloat("params","DB0")
DB1=cfg.getfloat("params","DB1")




plambda=cfg.getfloat("params","lambda")
print("lambda = %f" % plambda)

W=cfg.getfloat("params","W")

sigma=2/3*W/plambda
C0infA=cfg.getfloat("init","initCsA")
C1infA=cfg.getfloat("init","initClA")
C0infB=cfg.getfloat("init","initCsB")
C1infB=cfg.getfloat("init","initClB")

Thermo=ThermoDilute(esA,esB,0,elA,elB,elC)
Diff=DiffusionSpec(Dl1=DA1,Ds1=DA0,Dl2=DB1,Ds2=DB0)

Init=InitSpec(Clinf1=C1infA,Clinf2=C1infB,Csinf1=C0infA,Csinf2=C0infB)
x0=[-2.91461196,-0.35893441,7.21155476, 0.58294038, 0.52270732, 0.12272463] # only valid one


#  Init=InitSpec(Clinf1=23/scale,Clinf2=10/scale,Csinf1=3/scale,Csinf2=5/scale)
#  x0=[-0.8,-0.9,-1,  0.5,  0.5,  0.5] # only valid one
# print the params
print(Diff)
print(Init)
#  print('DlA={:2.3}, DlB={:2.3}, DsA={:2.3}, DsB={:2.3}'.format(float(DlA),float(DlB),float(DsA),float(DsB)))
#  print('ClinfA={:2.3}, ClinfB={:2.3}, CsinfA={:2.3}, CsinfB={:2.3}'.format(float(ClinfA),float(ClinfB),float(CsinfA),float(CsinfB)))

Solver = SolverStefan2p2c(Thermo, Diff, Init)
sol1 = Solver.solve([-4.7,2,2], 1000)
sol2 = Solver.solve([-1,2,2], 1000)
sol3 = Solver.solve([3.1,2,2], 1000)

Solver4p = SolverStefan4p2c(Thermo, Diff, Init)


sol4p = Solver4p.solve(x0, 1000)

#  Solver=AnalyticalSolverTernary.Solver2Phase1Int(ThermoP0,ThermoP1,DA0,DB0,DA1,DB1,C0infA,C0infB,C1infA,C1infB)
#  sol=Solver.search_all_xi((0,0,0), 1000,100)
#  sol=Solver.solution_list[0]
#  Solver.compute_all_params(sol)
#  Solver.print_params()
xi=sol1.xi
muA=sol1.muco1
muB=sol1.muco2


#  # write parameters in latex tabular
#  fileparams= open("Figs/params.tex","w")
#  d1={"m0":m0, "m1":m1,"m2":m2, "lambda01":lambda01,"lambda12":lambda12}
#  d2={"q0":q0,"q1":q1,"q2":q2,}
#  d3={"D0":D0, "D1":D1,"D2":D2, "xi01":xi01, "xi12":xi12}
#  d4={"C0inf":C0inf, "C2inf":C2inf,"mu01":mu01, "mu12":mu12}
#  txt="""\\renewcommand{\\arraystretch}{1.3}
#  \\begin{tabular}{|cc||cc||cc|c|cc||cc|}
#  \\cline{1-6} \\cline{2-6} \\cline{3-6} \\cline{4-6} \\cline{5-6} \\cline{6-6} \\cline{8-11} \\cline{9-11} \\cline{10-11} \\cline{11-11} 
#  \\multicolumn{2}{|c||}{Phase 0} & \\multicolumn{2}{c||}{Phase 1} & \\multicolumn{2}{c|}{Phase 2} &  & \\multicolumn{2}{c||}{Interface 0/1} & \\multicolumn{2}{c|}{Interface 1/2}\\tabularnewline
#  \\cline{1-6} \\cline{2-6} \\cline{3-6} \\cline{4-6} \\cline{5-6} \\cline{6-6} \\cline{8-11} \\cline{9-11} \\cline{10-11} \\cline{11-11}"""

#  txt+="$m_{{0}}$ & {m0} & $m_{{1}}$ & {m1} & $m_{{2}}$ & {m2} &  & $\\lambda_{{0/1}}$ & {lambda01} & $\\lambda_{{1/2}}$ & {lambda12}\\tabularnewline\n".format(**d1)
#  txt+="$Q_{{0}}$ & {q0} & $Q_{{1}}$ & {q1} & $Q_{{2}}$ & {q2} &  &  &  &  & \\tabularnewline\n".format(**d2)
#  txt+="$D_{{0}}$ & {D0} & $D_{{1}}$ & {D1} & $D_{{2}}$ & {D2} &  & ${{\\color{{red}}\\xi_{{0/1}}}}$ & {{\\color{{red}} {xi01:6.3f} }}& ${{\\color{{red}}\\xi_{{1/2}}}}$ &{{\\color{{red}} {xi12:6.3f}}}\\tabularnewline\n".format(**d3)
#  txt+="$C^{{0,\\infty}}$ & {C0inf} &  &  & $C^{{2,\\infty}}$ & {C2inf} &  & {{\\color{{red}}$\\mu_{{0/1}}$}} & {{\\color{{red}}{mu01:6.3f} }}& {{\\color{{red}}$\\mu_{{1/2}}$ }}& {{\\color{{red}}{mu12:6.3f}}}\\tabularnewline\n".format(**d4)
#  txt+="""\\cline{1-6} \\cline{2-6} \\cline{3-6} \\cline{4-6} \\cline{5-6} \\cline{6-6} \\cline{8-11} \\cline{9-11} \\cline{10-11} \\cline{11-11} 
#  \\end{tabular}"""
#  fileparams.write(txt)
#  fileparams.close()


nfiles=len(fnmatch.filter(os.listdir("."),"integrals*"))
dt=cfg.getfloat("run","dt")
t0=cfg.getfloat("init","init_time", fallback=0.0)
noutput=cfg.getint("run","nOutput")






xmin=cfg.getfloat("mesh","xmin")
xmax=cfg.getfloat("mesh","xmax")
Lx=xmax-xmin
Ly=cfg.getfloat("mesh","ymax")-cfg.getfloat("mesh","ymin")
nx=cfg.getfloat("mesh","nx")
ny=cfg.getfloat("mesh","ny")
dx=Lx/nx
dy=Ly/ny
if abs(dx-dy)>1e-10:
    print("ERROR: dx != dy")
Vcell=Lx*Ly/nx/ny
# init topology

x01_init=cfg.getfloat("init","x01",fallback=0.0)
x12_init=cfg.getfloat("init","x12",fallback=0.0)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read integral data

print("[%d s]" % (time.time()-script_start_time),"Reading integral data files")
data=[]

for i in range(nfiles):
    name="integrals_"+str(i)+".csv"
    with open(name) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:
            rowt=row
            rowt["step"]=i
            rowt["time"]=i*noutput*dt+t0
            data.append(rowt)


tsteps=[]
tx=[]
tw=[]
tie=[]
t_n_interf=[]
t_F=[]
times=[]
for row in data:
    area=float(row["Area"])
    tsteps.append(int(row["step"]))
    times.append(float(row["time"]))
    tx.append((float(row["Area"])-float(row["phi"])/Ly-(xmin)))
    gp=float(row["grand_potential"])
    ie=float(row["interfacial_energy"])
    F=float(row["free_energy"])
    
    nb_interf=math.floor(ie/(sigma)/2)*2+1
    t_n_interf.append(nb_interf)
    tw.append(((gp)/max(plambda,1e-8))/Ly)
    tie.append(((ie)/max(plambda,1e-8))/Ly)
    t_F.append(((F)/max(plambda,1e-8))/Ly)


tsteps=np.array(tsteps)
t_n_interf=np.array(t_n_interf)

times=np.array(times)
tx=np.array(tx)
tw=np.array(tw)
tie=np.array(tie)
t_F=np.array(t_F)
txold=tx-tx[0]+xi*t0**0.5



interf_color_choices=["blue","orange","green","red","purple","k","k","k"]
interf_color_choices=([mpl.colormaps["tab10"](i) for i in range(10)])
interf_colors=[]
for val in t_n_interf:
    interf_colors.append(interf_color_choices[min(val-1,6)])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read profile data

def get_profile(filenumber):
    profile_data={}
    profile_data["time"]=times[filenumber]
    filename="profile_"+str(filenumber)+".csv"
    with open(filename) as csvfile:
        reader=csv.DictReader(csvfile)
        for field in (reader.fieldnames):
            profile_data[field]=np.array([])
        for row in reader:
            if row["vtkValidPointMask"]=="1":
                for field in (reader.fieldnames):
                    profile_data[field]=np.append(profile_data[field],float(row[field]))
    return profile_data
    



tv=(tx[1:]-tx[:-1])/(dt*noutput)




markevery=1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# figure with total energy


tdw=(tw-tw[0])

tdwspeed=(tdw[1:]-tdw[0:-1])/(dt*noutput)



cutoff=1000
gpValues=[]
gpValues.append(times[1:cutoff]**0.5*sol1.Xi)
gpValues.append(times[1:cutoff]**0.5*sol2.Xi)
gpValues.append(times[1:cutoff]**0.5*sol3.Xi)
gpValues.append(times[1:cutoff]**0.5*sol4p.Xi)









cutoff=1000
gpVarValues=[]
gpVarValues.append(1/2/times[1:cutoff]**0.5*sol1.Xi)
gpVarValues.append(1/2/times[1:cutoff]**0.5*sol2.Xi)
gpVarValues.append(1/2/times[1:cutoff]**0.5*sol3.Xi)
gpVarValues.append(1/2/times[1:cutoff]**0.5*sol4p.Xi)


linthresh=min(abs(gpValues[0][1:cutoff]))
for gpvar in gpVarValues:
    linthresh=min(linthresh,min(abs(gpvar[1:cutoff])))
    
    
fig,ax_energy=plt.subplots(2,1,figsize=(8,8))


ax=ax_energy[0]
ax.plot(times[1:cutoff],gpVarValues[0],label="Solution 1",)
ax.plot(times[1:cutoff],gpVarValues[1],label="Solution 2",)
ax.plot(times[1:cutoff],gpVarValues[2],label="Solution 3",)
ax.plot(times[1:cutoff],gpVarValues[3],label="Solution 4 phases",)



#  ax.plot(times[1:cutoff],tdwspeed[1:cutoff], marker='',color="k")
ax.scatter(times[1:cutoff],tdwspeed[1:cutoff],label="Simulation",c=interf_colors[1:cutoff], marker='o')

ax.set_title("(a) Vitesse de variation du grand-potentiel")
ax.legend()

ax.set_xscale("log")
linthresh=1
ax.set_yscale("symlog", linthresh=linthresh)

ax.set_ylabel("Vitesse (échelle log symétrique, \nlinéaire entre -1 et 1)")

ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

#######################################
#######################################
#######################################
## plot free energy var
t_dF=t_F-t_F[0]

t_dF_var=(t_dF[1:]-t_dF[0:-1])/(dt*noutput)



cutoff=1000
FVarValues=[]
FVarValues.append(-1/2/times[1:cutoff]**0.5*sol1.Xi_F)
FVarValues.append(-1/2/times[1:cutoff]**0.5*sol2.Xi_F)
FVarValues.append(-1/2/times[1:cutoff]**0.5*sol3.Xi_F)
FVarValues.append(-1/2/times[1:cutoff]**0.5*sol4p.Xi_F)


linthresh=min(abs(gpValues[0][1:cutoff]))
for gpvar in FVarValues:
    linthresh=min(linthresh,min(abs(gpvar[1:cutoff])))
    
    
#  fig,ax=plt.subplots()
ax=ax_energy[1]
ax.plot(times[1:cutoff],FVarValues[0],label="Solution 1",)
ax.plot(times[1:cutoff],FVarValues[1],label="Solution 2",)
ax.plot(times[1:cutoff],FVarValues[2],label="Solution 3",)
ax.plot(times[1:cutoff],FVarValues[3],label="Solution 4 phases",)


ax.scatter(times[1:cutoff],-t_dF_var[1:cutoff],label="Simulation",c=interf_colors[1:cutoff], marker='o')

ax.set_title("(b) Vitesse de variation de l'énergie libre (valeur absolue)")
ax.legend()

ax.set_yscale("log")
ax.set_xscale("log")

ax.set_ylabel("Vitesse (échelle log)")
ax.set_xlabel("Temps (échelle log)")
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
fig.savefig("Figs/TemporalPlot_energy_variation.pdf")

#######################################
#######################################
#######################################
## plot excess dw with time
def phi0(l,t,xi,W):
    return 0.5*(1 + np.tanh(2* (l-xi)*(t**0.5) / W))
def pprime(phi):
    return 6*phi*(1-phi)
def plot_excess(sol,times):
    fig_excess_dw, ax_excess_dw=plt.subplots()
    fig_excess_dw_tot, ax_excess_dw_tot=plt.subplots()
    t_ex_dws_tot=np.array([])
    t_ex_dwl_tot=np.array([])
    t_ex_dws=np.array([])
    t_ex_dwl=np.array([])
    for t in times:
        unstable_dwltot=sol.dw * (sol.dw > 0) * (sol.tl > sol.xi)
        unstable_dwstot=sol.dw * (sol.dw < 0) * (sol.tl < sol.xi)
        
        r=1
        #  unstable_dwl=unstable_dwltot * (sol.tl < (sol.xi+r*W/t**0.5))
        #  unstable_dws=unstable_dwstot * (sol.tl > (sol.xi-r*W/t**0.5))
        unstable_dwl=unstable_dwltot * pprime(phi0(sol.tl,t,sol.xi,W))
        unstable_dws=unstable_dwstot * pprime(phi0(sol.tl,t,sol.xi,W))
        
        #  ax_excess_dw.plot(unstable_dws)
        excess_dwl_tot=abs(np.trapz(unstable_dwltot, x=sol.tl))*t**0.5
        excess_dws_tot=abs(np.trapz(unstable_dwstot, x=sol.tl))*t**0.5
        excess_dwl=abs(np.trapz(unstable_dwl, x=sol.tl))*t**0.5
        excess_dws=abs(np.trapz(unstable_dws, x=sol.tl))*t**0.5
        
        t_ex_dwl_tot=np.append(t_ex_dwl_tot,[excess_dwl_tot])
        t_ex_dws_tot=np.append(t_ex_dws_tot,[excess_dws_tot])
        t_ex_dwl=np.append(t_ex_dwl,[excess_dwl])
        t_ex_dws=np.append(t_ex_dws,[excess_dws])
    



    ax_excess_dw.plot(times,times**0, label="threshold")
    ax_excess_dw.plot(times, 3/4*plambda/W * t_ex_dwl, label="dwl")
    ax_excess_dw.plot(times, 3/4*plambda/W * t_ex_dws, label="dws")
    ax_excess_dw.legend()
    
    print("lambda condition = ", 1/max(3/4/W * t_ex_dws))
    
    ax_excess_dw_tot.plot(times,times**0, label="threshold")
    ax_excess_dw_tot.plot(times, 3/4*plambda/W * t_ex_dwl_tot, label="dwl tot")
    ax_excess_dw_tot.plot(times, 3/4*plambda/W * t_ex_dws_tot, label="dws tot")
    ax_excess_dw_tot.legend()
    #  ax_excess_dw.plot(time[1:], plambda*time[1:]**0.5 * abs(sol1.excess_dws))
    #  ax_excess_dw.plot(time[1:], plambda*time[1:]**0.5 * abs(sol1.excess_dwl))


#  alltimes=np.array([dt*i for i in range(20*noutput)])
#  plot_excess(sol1,alltimes[1:])


fig_ie,ax_ie=plt.subplots()
ax_ie.plot(tie/sigma)
ax_ie.plot(t_n_interf )
ax_ie.set_title("ie")


markevery=1


fig_diag, axs = plt.subplots(nrows=3, ncols=2,sharex=True, sharey=True,layout="constrained",figsize=(10,12))

for row in axs:
    for ax_diag in row:
        ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
        ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
        ax_diag.set_xlim([0,Thermo.Cl0_1])
        ax_diag.set_ylim([0,Thermo.Cl0_2])
        ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)
        
        ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*100) % 5 == 0])

        
        plot1 = ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="--", label="Solution 1")
        plot2 = ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="--", label="Solution 2")
        plot3 = ax_diag.plot(sol3.tC1, sol3.tC2, linewidth=2, linestyle="--", label="Solution 3")
        ax_diag.set_aspect('equal')
    
    
        plot4 = ax_diag.plot(sol4p.tC1, sol4p.tC2, linewidth=2, linestyle="--", label="Solution 4 phases", color=colors[3])
        ax_diag.legend()
    
axs[-1,0].set_xlabel("Composition $C_1$")
axs[-1,1].set_xlabel("Composition $C_1$")
axs[0,0].set_ylabel("Composition $C_2$")
axs[1,0].set_ylabel("Composition $C_2$")
axs[2,0].set_ylabel("Composition $C_2$")


def plot_profile_at_step(timestep, ax, title_prefix=""):
    tC1=get_profile(timestep)["composition"]
    tC2=get_profile(timestep)["cB"]
    ax.plot(tC1,tC2,linewidth=0,marker="o",markevery=markevery,zorder=-1, label="Simulation", color=colors[4])
    ax.set_title("%s Simulation à t=%5.4f" % (title_prefix, times[timestep]))
    ax.legend()
    ax_ie.axvline(x=timestep, color="k")




timestep=7 # just before nucleation number 1
plot_profile_at_step(timestep, axs[0,0], title_prefix="(a)")


timestep=30 # just before first transition to 5 interf
plot_profile_at_step(timestep, axs[0,1], title_prefix="(b)")


timestep=75 # first 5 interfs
plot_profile_at_step(timestep, axs[1,0], title_prefix="(c)")

timestep=150 # just after disparition of second interface
plot_profile_at_step(timestep, axs[1,1], title_prefix="(d)")

timestep=250 # 
plot_profile_at_step(timestep, axs[2,0], title_prefix="(e)")

timestep=400 # 
plot_profile_at_step(timestep, axs[2,1], title_prefix="(f)")




ax_diag.legend()
fig_diag.savefig("Figs/TernaryPlot.pdf")

plt.show()
