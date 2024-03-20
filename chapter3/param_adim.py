
import numpy as np
import forallpeople as si
import matplotlib.pyplot as plt


vals=["L", "W", "dt", "t", ]

vals=dict(sim={}, ref={}, real={})

vals["sim"]["L"]=1
vals["real"]["L"]=100 * 10**-9 * si.m
vals["ref"]["L"]=vals["real"]["L"]/vals["sim"]["L"]

vals["sim"]["D"]=1
vals["real"]["D"]=3 * 10**-22 * si.m**2 / si.s
vals["ref"]["D"]=vals["real"]["D"]/vals["sim"]["D"]

vals["sim"]["T"]=1e-8

vals["ref"]["T"]=vals["ref"]["L"]**2/vals["ref"]["D"]
vals["real"]["T"]=vals["sim"]["T"]*vals["ref"]["T"]


print(vals["real"]["T"])

vals["sim"]["K"]=1e-8

vals["ref"]["K"]=vals["ref"]["L"]**2/vals["ref"]["D"]
vals["real"]["K"]=vals["sim"]["T"]*vals["ref"]["T"]












