#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize


from ThermoK import ThermoK, mat


def plot_diagram(ax, TK, title=""):
    TK.compute_diagram_from_c()
    
    ax.plot([0,1,0,0], [0,0,1,0], color="k")
    
    ax.plot(TK.front0sA, TK.front0sB, color="b", label="Frontiere phase 0")
    
    ax.plot(TK.front1sA, TK.front1sB, color="r", label="Frontiere phase 1")
    

    # ~ x, rx, ry = fit_ellipse(TK.front0[Na2O], TK.front0[1-Na2O])
    
    # ~ rr = 10000
    # ~ x_coord = np.linspace(rx[0],rx[1],rr)
    # ~ y_coord = np.linspace(ry[0],ry[1],rr)
    # ~ X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
    # ~ Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
    # ~ plt.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('red'), linewidths=2)

    # ~ x, rx, ry = fit_ellipse(TK.front1[Na2O], TK.front1[1-Na2O])
    
    # ~ rr = 10000
    # ~ x_coord = np.linspace(rx[0],rx[1],rr)
    # ~ y_coord = np.linspace(ry[0],ry[1],rr)
    # ~ X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
    # ~ Z_coord = x[0] * X_coord ** 2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
    # ~ plt.contour(X_coord, Y_coord, Z_coord, levels=[1], colors=('blue'), linewidths=2)
    
    ax.plot([ceqB0, ceqB1], [ceqA0, ceqA1], color="green", label="Ceq", marker="o")
    
    for tl in TK.tielines:
        ax.plot(tl[0], tl[1], color = "gray", zorder= -10)
    # ~ ax.set_xlabel("C_A")
    # ~ ax.set_ylabel("C_B")
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_title(title)
    ax.legend()
    
    
if __name__ == "__main__":
    
    fig, axs=plt.subplots(nrows=2, ncols=3, figsize=[18, 10])
    
    
    # diagram with parabolas: deltaKinv has exaclty 1 non zero eigenvalue
    K0 = np.linalg.inv(mat(1,2,1))
    K1 = np.linalg.inv(mat(1,3,1))
    K0 = mat(2,1,0)
    K1 = mat(2,4,0)
    ceqA0 = 0.3
    ceqB0 = 0.3
    ceqA1 = 0.4
    ceqB1 = 0.4
    
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[0,0], TK, title="Parabolas")
    
    # ~ # diagram with straigth lines: deltaKinv = 0, degenerate parabolas
    K0 = mat(1,1,0)
    K1 = mat(1,1,0)
    ceqA0 = 0.1
    ceqB0 = 0.1
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[1,0], TK, title="Degenerate parabola: straight lines")
    
    # diagram with circles deltaKinv != 0, eigenvalues are equals
    # ~ K0 = mat(1,1,0)
    # ~ K1 = 2*mat(1,1,0)
    # ~ ceqA0 = 0.1
    # ~ ceqB0 = 0.1
    # ~ ceqA1 = 0.2
    # ~ ceqB1 = 0.2
    
    # ~ TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    # ~ plot_diagram(axs[0,1], TK)
    
    # diagram with ellipses: det(deltaKinv) > 0
    K0 = mat(1,1,0)
    K1 = mat(2,2,0.2)
    ceqA0 = 0.1
    ceqB0 = 0.1
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[0,1], TK, title="Ellipses (circle if K have same eigenvalues)")
    
    # ~ # diagram with a single point: degenerate ellipses. Numerical resolution is very unstable
    K0 = mat(1,1,0)
    K1 = 2*mat(1,1,0)
    ceqA0 = 0.2
    ceqB0 = 0.2
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[1,1], TK, title="Degenerate ellipses: a single point")
    
    # diagram with hyperbolas: det(deltaKinv) < 0. managed to get both section in the gibbs triangle, yay !
    K0 = mat(2,1,0)
    K1 = mat(1,2,1)
    ceqA0 = 0.1
    ceqB0 = 0.1
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[0,2], TK, title="Hyperbolas")
    
    # diagram with intersecting lines (degenerate hyperbola)
    K0 = mat(2,1,0)
    K1 = mat(1,2,0)
    ceqA0 = 0.1
    ceqB0 = 0.1
    ceqA1 = 0.2
    ceqB1 = 0.2
    
    TK = ThermoK(K0, ceqA0, ceqB0, K1, ceqA1, ceqB1)
    
    plot_diagram(axs[1,2], TK, title="Degenerate hyperbolas: intersecting lines")
    

    
   
    plt.savefig("catalog.pdf")
    plt.show()