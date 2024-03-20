
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
import json
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
times=[]
for row in data:
    area=float(row["Area"])
    tsteps.append(int(row["step"]))
    times.append(float(row["time"]))
    tx.append((float(row["Area"])-float(row["phi"])/Ly-(xmin)))
    gp=float(row["grand_potential"])
    ie=float(row["interfacial_energy"])
    tw.append((gp/max(plambda,1e-8))/Ly)
    # ~ tx01.append((float(row["phi3"])-float(data[-1]["phi3"]))*2)
    # ~ tx12.append(-(float(row["phi2"])-float(data[-1]["phi2"]))*2)

#  print(area,area/(Lx*Ly))
tsteps=np.array(tsteps)

times=np.array(times)
tx=np.array(tx)
tw=np.array(tw)
txold=tx-tx[0]+xi*t0**0.5


nfiles=1001


###############
# Define a function named 'test' that takes a list of numbers as input
def find_two_closest(nums):
    # Sort the unique elements in the list in ascending order
    s = sorted(set(nums))
    
    # Use list comprehension to find pairs of adjacent elements
    # Then, find the pair with the smallest difference
    return min([[a, b] for a, b in zip(s, s[1:])], key=lambda x: x[1] - x[0])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read profile data

print("[%d s]" % (time.time()-script_start_time),"Reading profile data files")

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
    
#  print("[%d s]" % (time.time()-script_start_time),"Reading profile data files")
profile_list=[]
nucleation_events=[]
disparition_events=[]


interflist_prev=[]
for filenumber in range(nfiles):
    print("Processing profile #%d" % (filenumber))
    profile_data=get_profile(filenumber)
    
    tphi=profile_data["phi"]
    interflist=[]
    for i in range(len(tphi)-1):
        phi0=tphi[i]
        phi1=tphi[i+1]
        if ((phi0<0.5 and phi1>0.5) or (phi0>0.5 and phi1<0.5)):
            x0=profile_data["Points:1"][i]
            x1=profile_data["Points:1"][i+1]
            posint=x0+(0.5-phi0)*(x1-x0)/(phi1-phi0)+xmin
            # ~ posint=
            interflist.append(posint)

    if len(interflist)>len(interflist_prev):
        if len(interflist)>1:
            a,b=find_two_closest(interflist)
            nucleation_events.append(dict(time=profile_data["time"], pos=(a+b)/2))
    elif len(interflist)<len(interflist_prev):
        if len(interflist_prev)>1:
            a,b=find_two_closest(interflist_prev)
            disparition_events.append(dict(time=profile_data["time"], pos=(a+b)/2))

    interflist_prev=interflist
    profile_data["interf_pos"]=interflist
    print("Found %d interfaces" % (len(interflist)))
    profile_list.append(profile_data)

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



tx=np.array([ profile_list[filenumber]["interf_pos"][0] for filenumber in range(nfiles) ])


tv=(tx[1:]-tx[:-1])/(dt*noutput)


interfs_list=[{"time":profile_list[filenumber]["time"], "interfs":profile_list[filenumber]["interf_pos"]} for filenumber in range(nfiles)]

#  print(interfs_list)

max_interfaces=len(max(interfs_list, key=lambda interf: len(interf["interfs"]))["interfs"])
#  print(max_interfaces)
interfaces_by_count=[{"times":[], "interfs":[[] for j in range(i+1)]} for i in range(max_interfaces)]
for interfs in interfs_list:
    count=len(interfs["interfs"])
    #  print(count)
    interfaces_by_count[count-1]["times"].append(interfs["time"])
    for j in range(count):
        interfaces_by_count[count-1]["interfs"][j].append(interfs["interfs"][j])
    
#  print(interfaces_by_count)


data=dict(interfaces_by_count=interfaces_by_count, nucleation_events=nucleation_events, disparition_events=disparition_events)

with open("interfaces_by_count.json", "w") as jsfile:
    jsfile.write(json.dumps(data))
    jsfile.close()
    
