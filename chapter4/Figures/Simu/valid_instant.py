
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
#  dir="case_insta_l300_metastable"
dir="case_instant_diff_valid"
#  dir="case_insta_l400_slow_nucleation"



os.chdir(dir)
    
# ~ print(datetime.datetime.fromtimestamp(os.path.getmtime("integrals_0.csv")))
# ~ print(os.path.getmtime("integrals_1.csv"))

if (os.path.exists("integrals_0.csv") and os.path.getmtime("integrals_0.csv")<os.path.getmtime("Simu_insta_maugis.xmf")):
    subprocess.call("rm integrals*", shell=True)
    subprocess.call("rm profile*", shell=True)
    subprocess.call("rm interf*", shell=True)
    
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

Dl=1e13
DA1=Dl
DB1=Dl

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
#  sol1 = Solver.solve([10,-2,-2], 1000)
sol1 = Solver.solve([-5,-0.2,-2], 1000)

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
    print("dx",dx,"dy",dy)
    print('correct ymax is ', dx*Ly/dy)
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
times=[]
for row in data:
    area=float(row["Area"])
    tsteps.append(int(row["step"]))
    times.append(float(row["time"]))
    #  tx.append((float(row["Area"])-float(row["phi"])/Ly-(xmin)))
    #  tx.append(xmin+Lx*((area-float(row["phi"]))/area))
    tx.append(1.0)
    gp=float(row["grand_potential"])
    #  ie=float(row["interfacial_energy"])
    tw.append((gp/max(plambda,1e-8))/Ly)
    # ~ tx01.append((float(row["phi3"])-float(data[-1]["phi3"]))*2)
    # ~ tx12.append(-(float(row["phi2"])-float(data[-1]["phi2"]))*2)

#  print(area,area/(Lx*Ly))
tsteps=np.array(tsteps)

times=np.array(times)
tx=np.array(tx)
tw=np.array(tw)
#  txold=tx-tx[0]+xi*t0**0.5

for i in range(nfiles):
    name="./interf_"+str(i)+".csv"
    #  print(name)
    with open(name) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:
            tx[i]=-(float(row["Points:1"])+xmin)
            #  print("txval",tx[i])
            #  print("txth",float(row["Points:1"])+xmin)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read profile data

def get_profile(filenumber):
    filenumber=(filenumber % nfiles)
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
    
def phi0(l,t,xi,W):
    return 0.5*(1 + np.tanh(2* (l-xi)*(t**0.5) / W))

def find_W_from_profile(profile):
    tphi=profile["phi"]
    tpx=profile["Points:1"]
    lb=phi0(-1,1,0,1)
    ub=phi0(1,1,0,1)
    lp=0
    up=1
    for i in range(len(tphi)-1):
        if (tphi[i]-lb)*(tphi[i+1]-lb)<=0:
            lp=tpx[i]
        if (tphi[i]-ub)*(tphi[i+1]-ub)<=0:
            up=tpx[i]
    print(lp, up)
    Wreal=abs(up-lp)/2
    print("real W",Wreal, "rate:", W/Wreal)
    return Wreal
#  print("[%d s]" % (time.time()-script_start_time),"Reading profile data files")
#  profile_list=[]
#  for filenumber in range(nfiles):                    
    #  profile_list.append(get_profile(filenumber))

# determine interface positions in profiles:
#  print("Determine interface positions in profiles")
#  for filenumber in range(nfiles):
    #  profile_data=profile_list[filenumber]
    #  tphi=profile_data["phi"]
    #  interflist=[]
    #  for i in range(len(tphi)-1):
        #  phi0=tphi[i]
        #  phi1=tphi[i+1]
        #  if ((phi0<0.5 and phi1>0.5) or (phi0>0.5 and phi1<0.5)):
            #  x0=profile_data["Points:1"][i]
            #  x1=profile_data["Points:1"][i+1]
            #  posint=x0+(0.5-phi0)*(x1-x0)/(phi1-phi0)+xmin
            #  # ~ posint=
            #  interflist.append(posint)
            
    #  profile_list[filenumber]["interf_pos"]=interflist



#  tx=np.array([ profile_list[filenumber]["interf_pos"][0] for filenumber in range(nfiles) ])


tv=(tx[1:]-tx[:-1])/(dt*noutput)




markevery=1
#plot interface positions

err_x=sum(abs(tx-(xi)*times**0.5))
errp_x=sum(abs(tx-(xi)*times**0.5))/sum(abs(tx))*100

print("Error on x :",err_x," (",errp_x,"%)")
err_v=sum(abs(tv-(xi)/2/times[1:]**0.5))
errp_v=sum(abs(tv-(xi)/2/times[1:]**0.5))/sum(abs(tv))*100


plt.figure()
plt.plot(times[1:],abs(tx[1:]), label="Simulation",linestyle="", marker='x', color=colors[0],markevery=markevery)
# ~ plt.plot(time,txold, label="x int",linestyle="", marker='x', color=colors[0],markevery=markevery)

plt.plot(times[1:], abs((xi)*times[1:]**0.5), color=colors[0], label="Analytical solution")
#  plt.plot(times[1:], (sol4p.xi_mean)*times[1:]**0.5, label="sol 4p")
#  plt.title("Interface position, L1 relative error is "+("%5.2f" %errp_x)+"$\%$" )
plt.title("Interface position" )
plt.legend()
plt.xscale("log")
plt.yscale("log")
#  plt.yscale('symlog', linthresh=min(abs(tx[1:]))/2)
#  plt.xlim([min((times[1:])),max((times[1:]))])
#  plt.ylim([min((tx[1:])),max((tx[1:]))])
plt.savefig("Figs/TemporalPlot_interface_position.pdf")


#plot speeds
plt.figure()
plt.plot(times[1:],abs(tv), label="Simulation",linestyle="", marker='x', color=colors[0],markevery=markevery)

plt.plot(times[1:], abs((xi)/2/times[1:]**0.5), color=colors[0],label="Analytique ($\\xi$=%g)" % xi)

plt.legend()
plt.title("Vitesse de l'interface en fonction du temps" )
plt.xscale("log")
plt.yscale("log")
#  plt.yscale('symlog', linthresh=min(abs(tv))/2)
#  plt.xlim([min((times[1:])),max((times[1:]))])
#  plt.ylim([min((tv[1:])),max((tv[1:]))])
plt.savefig("Figs/TemporalPlot_interface_speed.pdf")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# figure with total energy
XI1int=-xi*(elC-esC)

xi1,xi2,xi3=[-0.6517270481111545,-0.08026017339409147,1.6125526660195686]
XI3int=-(xi1+xi3-xi2)*(elC-esC)
tdwth=times**0.5*XI1int
tdwth3int=times**0.5*XI3int
figw=plt.figure()
tdw=(tw-tw[0]+tdwth[0])

tdwspeed=(tdw[1:]-tdw[0:-1])/(dt*noutput)
rgp=tdwth[-1]/tdw[-1]


gpvarValues=[]
gpvarValues.append(1/2/times[1:]**0.5*sol1.Xi)
plt.plot(times[1:],1/2/times[1:]**0.5*sol1.Xi,label="gp th",)
plt.plot(times[1:],tdwspeed,label="gp",linestyle="", marker='x',markevery=1)
#  plt.plot(tdw-tdw[-1]+time[-1]**0.5*sol4p.Xi,label="gp th",)
# ~ plt.plot(time, tdwth3int,label="gp th", color=colors[0])

err_gp=sum(abs(tdw-tdwth))
errp_gp=sum(abs(tdw-tdwth))/sum(abs(tdw))*100
plt.title("error on gp: %6.3f, %6.3f percent" % (err_gp,errp_gp))

plt.xscale("log")
linthresh=min(abs(gpvarValues[0]))
for gpvar in gpvarValues:
    linthresh=min(linthresh,min(abs(gpvar)))
#  linthresh=min(min(min(abs(1/2/time[1:]**0.5*sol1.Xi)),min(1/2/time[1:]**0.5*sol2.Xi)),min(min(1/2/time[1:]**0.5*sol3.Xi),min(1/2/time[1:]**0.5*sol4p.Xi)))
plt.yscale('symlog', linthresh=linthresh)
plt.legend()
plt.savefig("Figs/TemporalPlot_grand_potential.pdf")







markevery=1


fig_diag, ax_diag = plt.subplots()

ax_diag.plot([Thermo.Cl0_1, 0], [0,Thermo.Cl0_2], color="k")
ax_diag.plot([Thermo.Cs0_1, 0], [0,Thermo.Cs0_2], color="k")
ax_diag.plot([0,1,0,0], [0,0,1,0], color="k")
ax_diag.set_xlim([0,1])
ax_diag.set_ylim([0,1])
ax_diag.scatter([Init.Csinf1, Init.Clinf1], [Init.Csinf2,Init.Clinf2], color="k",zorder=5)

ax_diag.set_yticks([tick for tick in ax_diag.get_yticks() if int(tick*100) % 5 == 0])
ax_diag.set_xlabel("$C_1$")
ax_diag.set_ylabel("$C_2$")


ax_diag.set_aspect('equal')



plot1 = ax_diag.plot(sol1.tC1, sol1.tC2, linewidth=2, linestyle="--", label="Analytical solution")
#  plot2 = ax_diag.plot(sol2.tC1, sol2.tC2, linewidth=2, linestyle="--", label=("$\\xi$=% 5.3f" % (sol2.xi)))
#  plot3 = ax_diag.plot(sol3.tC1, sol3.tC2, linewidth=2, linestyle="--", label=("$\\xi$=% 5.3f" % (sol3.xi)))





#  plot4 = ax_diag.plot(sol4p.tC1, sol4p.tC2, linewidth=2, linestyle="--", label="4p, $\\Xi=%f$" % (sol4p.Xi))

    

timestep=-1
profile=get_profile(timestep)
tC1=profile["composition"]
tC2=profile["cB"]
tphi=profile["phi"]

Wreal=find_W_from_profile(profile)

c1,c2=sol1.get_diffuse_profile(W, times[timestep])
ax_diag.plot(c1,c2, label="Analytical solution interpolated")

#  Init_corrected=InitSpec(Clinf1=tC1[-1], Clinf2=tC2[-1], Csinf1=Init.Csinf1, Csinf2=Init.Csinf2)
#  Solver_corrected= SolverStefan2p2c(Thermo, Diff, Init_corrected)
#  sol_corrected = Solver_corrected.solve([-4.7,2,2], 1000)
#  ax_diag.plot(sol_corrected.tC1, sol_corrected.tC2, linewidth=2, linestyle="--", label=("$\\xi$=% 5.3f" % (sol_corrected.xi)))

ax_diag.plot(tC1,tC2,linewidth=0,marker="o",markevery=markevery,zorder=-1, label="Simulation at t=%3.1e" % (times[timestep]))
#  ax_diag.plot(tC1,tC2,linewidth=1,marker="",markevery=markevery,zorder=-1, label="Simulation at %5.3e" % (times[timestep]))
#  ax_diag.scatter(tC1,tC2,marker="o",zorder=-1, c=(1-tphi), cmap="coolwarm")

print(tC1)
#  timestep=70
#  tC1=get_profile(timestep)["composition"]
#  tC2=get_profile(timestep)["cB"]
#  ax_diag.plot(tC1,tC2,linewidth=0,marker="o",markevery=markevery,zorder=-1, label="Simulation at %d" % (timestep))






ax_diag.legend()
fig_diag.savefig("Figs/TernaryPlot.pdf")

plt.show()
