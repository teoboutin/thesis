
import fnmatch
import os
import datetime
import subprocess
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.scale
from matplotlib.patches import ConnectionPatch
import numpy as np
import configparser
import AnalyticalSolverTernary
import time
import json
tab20=plt.get_cmap("tab20")
from solvers.ThermoDilute import ThermoDilute
from solvers.solver_stefan_2p2c import Stefan2p2cSolution,SolverStefan2p2c,DiffusionSpec,InitSpec
from solvers.solver_stefan_4p2c import Stefan4p2cSolution,SolverStefan4p2c
from solvers.envelope import envelope3 as envelope
from solvers.envelope import envelope_monotonous

import scienceplots
#  plt.style.use(['science','light'])
plt.style.use(['science','vibrant'])
#  plt.style.use(['science','std-colors'])

fsize=plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = (fsize[0]*2,fsize[1]*2)


colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
colors=([mpl.colormaps["Set1"](i) for i in range(10)])
colors=([mpl.colormaps["tab10"](i) for i in range(10)])
#  colors=(["red"])

script_start_time=time.time()

#  dir="test"
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

#  nfiles=100

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
time=[]
for row in data:
    area=float(row["Area"])
    tsteps.append(int(row["step"]))
    time.append(float(row["time"]))
    tx.append((float(row["Area"])-float(row["phi"])/Ly-(xmin)))
    gp=float(row["grand_potential"])
    ie=float(row["interfacial_energy"])
    tw.append((gp/max(plambda,1e-8))/Ly)
    # ~ tx01.append((float(row["phi3"])-float(data[-1]["phi3"]))*2)
    # ~ tx12.append(-(float(row["phi2"])-float(data[-1]["phi2"]))*2)

#  print(area,area/(Lx*Ly))
tsteps=np.array(tsteps)

time=np.array(time)
tx=np.array(tx)
tw=np.array(tw)
txold=tx-tx[0]+xi*t0**0.5



with open("interfaces_by_count.json", "r") as jsfile:
    #  print(jsfile.read())
    data=json.loads(jsfile.read())
    interfaces_by_count=data["interfaces_by_count"]
    nucleation_events=data["nucleation_events"]
    disparition_events=data["disparition_events"]
    jsfile.close()

markevery=1

max_interfaces=len(interfaces_by_count)
print("max interfaces", max_interfaces)

min_val=np.inf
for i in range(max_interfaces):
    for j in range(len(interfaces_by_count[i]["interfs"])):
        for val in interfaces_by_count[i]["interfs"][j]:
            min_val=min(min_val, abs(val))

print(min_val)
min_val=5e-4

max_time=0
for i in range(max_interfaces):
    for val in interfaces_by_count[i]["times"]:
        max_time=max(max_time, abs(val))


time=np.array([val for val in time if val<=max_time])

#  fig, ax = plt.subplots()
fig, ax = plt.subplots(layout="constrained")
#  fig.tight_layout()
# inset Axes....
x1, x2, y1, y2 = 1e-4, 2.5e-4, -1.1e-1, -3e-2  # subregion of the original image
axins = ax.inset_axes(
    [0.4, 0.3, 0.3, 0.3],
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

    
ax.set_title("Interfaces positions" )


labels=["1st interface", "2nd interface", "3rd interface", "4th interface", "5th interface"]

## plot 1 int phase

n_int=min(10,max_interfaces)
for max_int in range(n_int):
    time_loc=np.array(interfaces_by_count[max_int]["times"])
    print(len(time_loc))
    if len(time_loc)>0:
        for j in range(max_int+1):
            print("%d interfaces, %de int" % (max_int,max_int-j+1))
            axins.plot(time_loc,interfaces_by_count[max_int]["interfs"][max_int-j],linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            if max_int < n_int-1:
                ax.plot(time_loc,interfaces_by_count[max_int]["interfs"][max_int-j],linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            else:
                ax.plot(time_loc,interfaces_by_count[max_int]["interfs"][max_int-j],label=labels[j],linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])


ax.annotate("Nucleation", (interfaces_by_count[2]["times"][0], (interfaces_by_count[2]["interfs"][1][0]+interfaces_by_count[2]["interfs"][2][0])/2), xytext=(4e-6, 0), arrowprops=dict(arrowstyle="->"),alpha=0)
ax.annotate("Nucleation", (interfaces_by_count[4]["times"][0], (interfaces_by_count[4]["interfs"][1][0]+interfaces_by_count[4]["interfs"][0][0])/2), xytext=(4e-6, 0),xycoords=(axins.transData), textcoords=(ax.transData), arrowprops=dict(arrowstyle="->"),zorder=100)

    
ax.plot(time[1:], (sol4p.xi12)*time[1:]**0.5, zorder=100, color="k", label="Analytical")
ax.plot(time[1:], (sol4p.xi23)*time[1:]**0.5, zorder=100, color="k")
ax.plot(time[1:], (sol4p.xi34)*time[1:]**0.5, zorder=100, color="k")
ax.plot(time[1:], (xi)*time[1:]**0.5, zorder=100, color="k")

axins.plot(time[1:], (sol4p.xi12)*time[1:]**0.5, zorder=100, color="k")
axins.plot(time[1:], (sol4p.xi23)*time[1:]**0.5, zorder=100, color="k")
axins.plot(time[1:], (sol4p.xi34)*time[1:]**0.5, zorder=100, color="k")
axins.plot(time[1:], (xi)*time[1:]**0.5, zorder=100, color="k")
        
        
t=1.2e-5
ax.annotate("$\\xi$", (t, (xi)*t**0.5), xytext=(t, -1.1e-1), arrowprops=dict(arrowstyle="->"),ha="right")
t=1.2e-5
ax.annotate("$\\xi_{12}$", (t, (sol4p.xi12)*t**0.5), xytext=(t*1.1, -0.5e-2), arrowprops=dict(arrowstyle="->"), ha="left")
t=3e-3
ax.annotate("$\\xi_{23}$", (t, (sol4p.xi23)*t**0.5), xytext=(t, 2e-4), arrowprops=dict(arrowstyle="->"))
t=6e-3
ax.annotate("$\\xi_{34}$", (t, (sol4p.xi34)*t**0.5), xytext=(1.2e-2, 1e-3), arrowprops=dict(arrowstyle="->", relpos=(0.8,1)), ha="right")

ax.legend(frameon=True, loc='upper left',fancybox=True, facecolor='white', framealpha=1).set_zorder(102)
ax.set_xscale("log")
ax.set_yscale('symlog', linthresh=min_val/2)
axins.set_xscale("log")
axins.set_yscale('symlog', linthresh=min_val/2)

ax.set_xlabel("Time (log scale)")

ax.indicate_inset_zoom(axins, edgecolor="black")
plt.savefig("Figs/TemporalPlot_multi_interfaces_positions.pdf")

##################################################
## plot the other way
fig, ax = plt.subplots(layout="constrained")
#  fig, ax = plt.subplots()
#  fig.set_tight_layout(True)

ax.set_title("Interfaces positions" )


labels=["1st interface", "2nd interface", "3rd interface", "4th interface", "5th interface", "6th interface", "7th interface", "8th interface", "9th interface", "10th interface"]

## plot 1 int phase

n_int=min(10,max_interfaces)
for max_int in range(n_int):
    time_loc=np.array(interfaces_by_count[max_int]["times"])
    print(len(time_loc))
    if len(time_loc)>0:
        for j in range(max_int+1):
            print("%d interfaces, %de int" % (max_int,max_int-j+1))
            #  axins.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            if max_int < n_int-1:
                ax.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
                #  ax.errorbar(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,xerr=10*W, yerr=0,elinewidth=1,ecolor="k",linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            else:
                ax.plot(interfaces_by_count[max_int]["interfs"][max_int-j],time_loc,label=labels[min(j, len(labels)-1)],linestyle="", marker='o',markevery=markevery, color=colors[min(j, len(colors)-1)])
            
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
                
                #  ax.fill_betweenx(time_loc, interfaces_by_count[max_int]["interfs"][j], x2=interfaces_by_count[max_int]["interfs"][j+1], color=color,**kwargs)
            #  ax.fill_betweenx(time_loc, interfaces_by_count[max_int]["interfs"][0], x2=-0.24, color=csol,**kwargs)
            #  ax.fill_betweenx(time_loc, interfaces_by_count[max_int]["interfs"][-1], x2=0.24, color=cliq,**kwargs)
                
#  time_loc=np.array(interfaces_by_count[2]["times"])
#  ax.fill_betweenx(time_loc, interfaces_by_count[2]["interfs"][0], x2=interfaces_by_count[2]["interfs"][1], color="b", alpha=0.5, step="post")

ax.plot((sol4p.xi12)*time**0.5,time, zorder=100, color="k", label="Analytical")
ax.plot((sol4p.xi23)*time**0.5,time, zorder=100, color="k")
ax.plot((sol4p.xi34)*time**0.5,time, zorder=100, color="k")
ax.plot((xi)*time**0.5,time, zorder=100, color="k")


#  t=0.0015
#  ax.annotate("$\\xi$", ((xi)*t**0.5, t), xytext=(-0.2, 0.001), arrowprops=dict(arrowstyle="->"),ha="right")
#  t=0.001
#  ax.annotate("$\\xi_{12}$", ( (sol4p.xi12)*t**0.5, t), xytext=(-0.1, 0.0015), arrowprops=dict(arrowstyle="->"), ha="left")
#  t=0.002
#  ax.annotate("$\\xi_{23}$", ( (sol4p.xi23)*t**0.5, t), xytext=(-0.05, 0.002), arrowprops=dict(arrowstyle="->"))
#  t=0.001
#  ax.annotate("$\\xi_{34}$", ( (sol4p.xi34)*t**0.5, t), xytext=(1.2e-2, 1e-3), arrowprops=dict(arrowstyle="->", relpos=(0.8,1)), ha="right")

t=time[-1]
ax.text((xi)*t**0.5, t+0.00003, "$\\xi\\sqrt{t}$", va="bottom", ha="center")
ax.text((sol4p.xi12)*t**0.5, t+0.00003, "$\\xi_{12}\\sqrt{t}$", va="bottom", ha="center")
ax.text((sol4p.xi23)*t**0.5, t+0.00003, "$\\xi_{23}\\sqrt{t}$", va="bottom", ha="center")
ax.text((sol4p.xi34)*t**0.5, t+0.00003, "$\\xi_{34}\\sqrt{t}$", va="bottom", ha="center")


#  nucl_points=[((interfaces_by_count[2]["interfs"][1][0]+interfaces_by_count[2]["interfs"][2][0])/2, interfaces_by_count[2]["times"][0]), ((interfaces_by_count[4]["interfs"][1][0]+interfaces_by_count[4]["interfs"][0][0])/2, interfaces_by_count[4]["times"][0])]

#  alpha=1
#  for nucl_p in nucl_points:
    #  ax.annotate("Nucleation",nucl_p , xytext=(-0.05, -0.00025), arrowprops=dict(arrowstyle="->"),alpha=alpha, va="top")
    #  alpha=0
#  ax.annotate("Nucleation", , xytext=(-0.2, 0.0), arrowprops=dict(arrowstyle="->"),zorder=100)

s=1*plt.rcParams['lines.markersize'] ** 2
ax.scatter([e["pos"] for e in nucleation_events], [e["time"] for e in nucleation_events], marker="d", s=s, label="Nucleation event", color="k", zorder=1000)
ax.scatter([e["pos"] for e in disparition_events], [e["time"] for e in disparition_events], marker="d", s=s, label="Disparition event", color="yellow", zorder=1000)
#  do_label=True
#  for nucl_event in nucleation_events:
    #  t=nucl_event["time"]
    #  pos=nucl_event["pos"]
    #  if do_label:
        #  ax.axhline(y=t, color="grey", linestyle="--", label="Nucleation event")
        #  do_label=False
    #  else:
        #  ax.axhline(y=t, color="grey", linestyle="--")
        #  ax.axvline(x=pos, color="grey", linestyle="--")

ax.legend(frameon=True, loc='lower right',fancybox=True, facecolor='white', framealpha=1).set_zorder(102)
#  ax.set_yscale("log")
#  ax.set_xscale('symlog', linthresh=min_val/2)
ax.set_xlabel("Position")
ax.set_ylabel("Time")
#  ax.set_ylim([-0.0005, 0.003])
plt.savefig("Figs/TemporalPlot_multi_interfaces_positions.pdf")

plt.show()
