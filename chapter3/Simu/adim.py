import os
import shutil
import configparser
import forallpeople as si
def phys_to_scientific(phys):
    val,unit=str(phys).split(" ")
    # ~ print(val)
    # ~ print(unit)
    sv=f"{float(val):.3e}"
    return (sv+" "+unit)
    

import numpy as np

print(si.__file__)
path="/home/catA/tb266682/.local/lib/python3.8/site-packages/forallpeople/environments"

shutil.copyfile("./adim_units.json", path+"/adim_units.json")

si.environment('adim_units', top_level=False)

cfg=configparser.ConfigParser()


filename="simu1.ini"
filename="simu2.ini"

cfg.read(filename)

ly=(cfg.getfloat("mesh","ymax")-cfg.getfloat("mesh","ymin"))
lx=(cfg.getfloat("mesh","xmax")-cfg.getfloat("mesh","xmin"))
vv=cfg.getfloat("ccparams","virtual_volume")



Angstrom=10**-10 * si.m
nm=10**-9 * si.m
Dunit=si.m*si.m/si.s

Values=["Length", "Diff", "Time", "Time_per_output", "Time_per_step","E_volu","sigma","TimePF", "S_V", "total_volume"]
types=["sim","ref","reel"]

# ~ dd=dict([(t,0) for t in types])
data=dict([(name,dict([(t,0) for t in types])) for name in Values])


data["Length"]["reel"]=100 * nm
data["Diff"]["reel"]=3*10**-22*Dunit




data["Length"]["sim"]=(ly)*si.Dimensionless
data["Length"]["ref"]=data["Length"]["reel"]/data["Length"]["sim"]

data["Diff"]["sim"]=cfg.getfloat("params","DB0")*si.Dimensionless
data["Diff"]["ref"]=data["Diff"]["reel"]/data["Diff"]["sim"]

data["Time"]["sim"]=cfg.getfloat("run","dt")*cfg.getfloat("run","nStepmax")*si.Dimensionless
data["Time"]["ref"]=(data["Length"]["ref"]**2/data["Diff"]["ref"]).to(unit_name="day")
data["Time"]["reel"]=data["Time"]["sim"]*data["Time"]["ref"].to(unit_name="day")

data["Time_per_output"]["sim"]=cfg.getfloat("run","dt")*cfg.getfloat("run","nOutput")*si.Dimensionless
data["Time_per_output"]["ref"]=data["Time"]["ref"]
data["Time_per_output"]["reel"]=data["Time_per_output"]["sim"]*data["Time_per_output"]["ref"].to(unit_name="minute")

data["Time_per_step"]["sim"]=cfg.getfloat("run","dt")*si.Dimensionless
data["Time_per_step"]["ref"]=data["Time"]["ref"]
data["Time_per_step"]["reel"]=data["Time_per_step"]["sim"]*data["Time_per_step"]["ref"].to(unit_name="seconds")

data["TimePF"]["sim"]=cfg.getfloat("params","W")**2/cfg.getfloat("params","Mphi")*si.Dimensionless
data["TimePF"]["ref"]=data["Time"]["ref"]
data["TimePF"]["reel"]=data["TimePF"]["sim"]*data["TimePF"]["ref"]
kB=1.38*10**-23 * si.J / si.K
Na=6.022*10**23 * 1/si.mol
Rgp=kB*Na
T=(273.15+100)*si.K
rm=1.22*Angstrom
Vm=4/3*np.pi*rm**3
V_molaire_eau = 18*1e-6 * si.m**3 /si.mol
V_molaire_silicium = 12.06*1e-6*si.m**3/si.mol
V_molaire_bore= 61.8*1e-3 * si.kg / si.mol / (1.435*1e-3*si.kg/(1e-6*si.m**3))
V_moleculaire_eau=V_molaire_eau/Na
V_moleculaire_bore=V_molaire_bore/Na
V_moleculaire_silicium=V_molaire_silicium/Na

V_molec_franck= 20*1e-6 * si.m**3 /si.mol / Na
T_franck = (273.15+20)*si.K
print("V_molaire_eau",V_molaire_eau)
print("V_molaire_bore",V_molaire_bore)
print("V_molaire_silicium",V_molaire_silicium)
#  data["E_volu"]["ref"]=kB*T/V_moleculaire_eau
data["E_volu"]["sim"]=1.0*si.Dimensionless
#  data["E_volu"]["reel"]=data["E_volu"]["sim"]*data["E_volu"]["ref"]
#  data["E_volu"]["reel"]=kB*T/V_moleculaire_eau
data["E_volu"]["reel"]=kB*T/V_moleculaire_eau
data["E_volu"]["reel"]=kB*T_franck/V_molec_franck
data["E_volu"]["ref"]=data["E_volu"]["reel"]/data["E_volu"]["sim"]

data["S_V"]["ref"]=1/data["Length"]["ref"]


print("real volume", (ly*(lx-ly)))
lv= (ly*(lx-ly))+vv
data["S_V"]["sim"]=ly/lv*si.Dimensionless
data["S_V"]["reel"]=(data["S_V"]["sim"]*data["S_V"]["ref"])


data["total_volume"]["ref"]=((data["Length"]["ref"])**3).to(unit_name="m")
data["total_volume"]["sim"]=((ly*lx)+vv)*si.Dimensionless
data["total_volume"]["reel"]=(data["total_volume"]["sim"]*data["total_volume"]["ref"])


data["sigma"]["sim"]=2/3*cfg.getfloat("params","W")/cfg.getfloat("params","lambda")*si.Dimensionless
data["sigma"]["ref"]=(data["E_volu"]["ref"]*data["Length"]["ref"])

data["sigma"]["reel"]=(data["sigma"]["sim"]*data["sigma"]["ref"]).to(unit_name="sigma")

width=18
s=''.rjust(width,' ')+"|"
for t in types:
    s+=t.rjust(width,' ')+"|"
print(s)
for name in Values:
    s=name+" :"
    s=s.ljust(width,' ')+"|"
    for t in types:
        s+=  str(data[name][t]).rjust(width,' ')
        # ~ s+=  phys_to_scientific(data[name][t]).rjust(width,' ')
        s+= "|"
    print(s)
    
print((2*si.m).to(unit_name="cm"))

M=62/1000*si.kg/si.mol
Mv=1430 * si.kg/(si.m**3)
Vm=M/Mv
Va=Vm/Na
r=(3/4*Va/np.pi).sqrt(3)
print(Vm*si.mol)

print((Va)/Angstrom**3, "A**3")
print(r/Angstrom, "A")
