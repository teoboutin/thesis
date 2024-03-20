
import numpy as np
import forallpeople as si
import matplotlib.pyplot as plt
gpmol=1e-3 #* si.kg/si.mol
Mm=dict()
Mm["H3BO3"]=61.833 * gpmol
Mm["H2O"]=19.015 * gpmol


T=(np.array([0,15,20,25,50,75])+273.15)
H3BO3_solubility_wtp=np.array([2.7,4.17,4.65,5.44,10.24,17.41])/100




def convert_wtp_to_molp(wtpB):
    nH3BO3=wtpB / Mm["H3BO3"]
    nH2O=(1-wtpB) / Mm["H2O"]
    return nH3BO3/ (nH3BO3+nH2O)
    
for i in range(len(T)):
    print(str(T[i]).ljust(15), "->", convert_wtp_to_molp(H3BO3_solubility_wtp[i]))

H3BO3_solubility_molp=convert_wtp_to_molp(H3BO3_solubility_wtp)
plt.plot(T,H3BO3_solubility_molp)
plt.xscale("log")
plt.yscale("log")



a,b=np.polyfit(np.log(T),np.log(H3BO3_solubility_molp),1)

def solubility_th(t):
    return np.exp(a*np.log(t)+b)
plt.plot(T, solubility_th(T))

tk=273.15+90
print(str(tk).ljust(15), "->", solubility_th(tk))
print()


plt.show()


