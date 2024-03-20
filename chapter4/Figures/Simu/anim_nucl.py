import math
import fnmatch
import os
import datetime
import subprocess
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.scale

from matplotlib.animation import FuncAnimation

import json
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



print("read interfaces")

with open("interfaces_by_count.json", "r") as jsfile:
    #  print(jsfile.read())
    data=json.loads(jsfile.read())
    interfaces_by_count=data["interfaces_by_count"]
    nucleation_events=data["nucleation_events"]
    disparition_events=data["disparition_events"]
    jsfile.close()

max_interfaces=len(interfaces_by_count)
print("max interfaces", max_interfaces)

min_val=np.inf
for i in range(max_interfaces):
    for j in range(len(interfaces_by_count[i]["interfs"])):
        for val in interfaces_by_count[i]["interfs"][j]:
            min_val=min(min_val, abs(val))

print(min_val)
min_val=5e-4

labels=["1st interface", "2nd interface", "3rd interface", "4th interface", "5th interface"]



markevery=1



fig_anim = plt.figure(figsize=(16, 9), layout="tight")
fig_anim.tight_layout()
spec = fig_anim.add_gridspec(2, 2,height_ratios=[1,3])


axC= fig_anim.add_subplot(spec[0,1])
ax_phi = fig_anim.add_subplot(spec[0,0])
ax_diag = fig_anim.add_subplot(spec[1,1])
ax_multi = fig_anim.add_subplot(spec[1,0])

#  axC = ax_phi.twinx()  # instantiate a second Axes that shares the same x-axis

color = tab10(0)
#  axC.set_ylabel('Composition')  # we already handled the x-label with ax1
#  axC.tick_params(axis='y', labelcolor=color)

plot_ana, = ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="--", label="Solution analytique")
    
plot_profile,=ax_diag.plot([],[],linewidth=0,marker="o",markevery=markevery,zorder=-1, label="Simulation", color=colors[4])
plot_phi,=ax_phi.plot([],[],label=r'$\phi$')
plot_C1,=axC.plot([],[],label=r"$C_1$", color=tab10(2))
plot_C2,=axC.plot([],[],label=r"$C_2$", color=tab10(3))


n_int=min(10,max_interfaces)
for max_int in range(n_int):
    time_loc=np.array(interfaces_by_count[max_int]["times"])/dt
    print(len(time_loc))
    if len(time_loc)>0:
        for j in range(max_int+1):
            print("%d interfaces, %de int" % (max_int,max_int-j+1))
            #  axins.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            if max_int < n_int-1:
                ax_multi.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
                #  ax.errorbar(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,xerr=10*W, yerr=0,elinewidth=1,ecolor="k",linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            else:
                ax_multi.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,label=labels[min(j, len(labels)-1)],linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            
            exclude_times=[e["time"] for e in nucleation_events]
            where=np.append([True], [not(t in exclude_times) for t in time_loc[1:]])
            kwargs=dict(where=where, alpha=1, edgecolor="none")
            csol="tab:olive"
            cliq="tab:cyan"
            if j<max_int:
                if j%2==0:
                    color=cliq
                else:
                    color=csol
                


artists=[
plot_ana,
plot_profile,
plot_phi,
plot_C1,
plot_C2
]

ax_phi.legend(loc="upper left")
axC.legend(loc="center right")

ax_phi.set_title("Indicatrice de phase")
ax_multi.set_title("Position des interfaces (axe vertical : temps)")
axC.set_title("Profils de composition")
ax_diag.set_title("Profil de composition dans le diagramme de phase")

fig_anim.tight_layout()


print("[%d s]" % (time.time()-script_start_time),"Reading profiles data files")

stride=1
anim_end=500
profile_list=[]
for i in range(0,anim_end,stride):
    profile_list.append(get_profile(i))

print(nfiles)
print(profile_list[0].keys())
tX=profile_list[0]['Points:1']
print(min(tX), max(tX))

def init():
    
    ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
    ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
    ax_diag.set_xlim([0,Thermo.Cl0_1])
    ax_diag.set_ylim([0,Thermo.Cl0_2])
    ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)
    
    ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*100) % 5 == 0])
    ax_diag.set_aspect('equal')
    
    
    #  plot2 = ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="--", label="Solution 2")
    #  plot3 = ax_diag.plot(sol3.tC1, sol3.tC2, linewidth=2, linestyle="--", label="Solution 3")
    
    
    #  plot4 = ax_diag.plot(sol4p.tC1, sol4p.tC2, linewidth=2, linestyle="--", label="Solution 4 phases", color=colors[3])
    ax_diag.legend()
    ax_multi.legend()
    
    xr=[-0.25,0.25]
    
    ax_phi.set_ylim([-0.05,1.05])
    ax_phi.set_xlim(xr)
    axC.set_xlim(xr)
    ax_multi.set_xlim(xr)
    
    
    axC.set_ylim([0, 0.25])
    
    return artists



def update(frame):
    
    plot_ana.set_data(sol1.tC1, sol1.tC2)
    
    profile=profile_list[frame]
    tC1=profile["composition"]
    tC2=profile["cB"]
    
    plot_profile.set_data(tC1,tC2)
    
    tX=profile['Points:1']-1
    plot_phi.set_data(tX,profile['phi'])
    plot_C1.set_data(tX,tC1)
    plot_C2.set_data(tX,tC2)
    
    ax_multi.set_ylim([-1, noutput*stride*frame])
    
    
    #  fig_anim.suptitle(noutput*stride*frame)
    
    return artists

ani = FuncAnimation(fig_anim, update, frames=range(0,len(profile_list)),
                    init_func=init, blit=False, repeat=False, interval=40)

# saving to m4 using ffmpeg writer 
print("saving video")
writervideo = mpl.animation.FFMpegWriter() 
ani.save('Figs/anim.mp4', writer=writervideo) 
#  ani.save(filename="Figs/anim.gif", writer="pillow")
print("saved video")

plt.show()
