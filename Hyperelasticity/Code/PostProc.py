import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

font = {'family': 'serif', 'weight': 'normal', 'size': 18}
FontLegend = {'family': 'serif', 'weight': 'normal', 'size': 15}
plt.rc('font', **font)
colors= [ 'black', 'red', "blue","royalblue", "lightblue", ]
style = ['-', '-', '--', '-.']
# Read in text file
DirStart = '../Results'
DirList = ["/EnergyComp/S_20_Nx_15_N_15_nu_4.9e-01",
           "/EnergyPenalty/S_20_Nx_15_N_15_nu_4.9e-01",
           "/WeakFormU/S_20_Nx_15_N_15_nu_4.9e-01",
           "/WeakFormUp/S_20_Nx_10_N_10_nu_4.9e-01"]
# Title = [r"$N_x = 10, \, N_z = 10$", r"$N_x=15, \, N_z = 20$", r"$N_x=20, \, N_z = 30$"]

DirLen = len(DirList)
# Title = [r"$\nu = 0.4$", r"$\nu = 0.45$", r"$\nu = 0.49$"]
# LegTitle = [r"$0.4$", r"$0.45$", r"$0.49$"]
Title = ["Compressible", "Penalty", r"Penalty Weak Form: $\mathbf{u}$", r"Weak Form: $\mathbf{u}$, $p$"]
LegTitle = ["Compressible", "Penalty", r"Penalty Weak Form: $\mathbf{u}$", r"Weak Form: $\mathbf{u}$, $p$"]

ThetaArray = np.linspace(0, (np.pi)/3, 20) # For expressions
RampArray = np.linspace(0, 0.5, 20)   # For expressions

DirNo = 0
if DirLen == 2:
    fig, ax = plt.subplots(1,DirLen, figsize=(10,5), sharey=True)
    fig2, ax2 = plt.subplots(1,1, figsize=(5,5), sharey=True)
elif DirLen == 3:
    fig, ax = plt.subplots(1,DirLen, figsize=(15,5), sharey=True)
    fig2, ax2 = plt.subplots(1,1, figsize=(6,5), sharey=True)
else:
    fig, ax = plt.subplots(1,DirLen, figsize=(15,5), sharey=True)
    fig2, ax2 = plt.subplots(1,1, figsize=(6,5), sharey=True)

ax2.set_xlim((0,20))
plt.axhline(1.0,color='black', linestyle='--', alpha=0.5)
for Dir in DirList:
    # Step, Scale of expression, Displacement x, y, and z
    df1 = pd.read_table(DirStart + Dir + '/PostProc.txt', delimiter=' ', header=None)

    df1_convert = df1.to_numpy()
    Step = df1_convert[:,0]
    Scale = df1_convert[:,1]
    Disp2 = df1_convert[:,2]
    Disp3 = df1_convert[:,3]
    DetF = df1_convert[:,4]

    Exp1, Exp2 = [], []
    for Ind, Val in enumerate(ThetaArray):
        ExpNum1 = RampArray[Ind]*(0.5 + 0.5*np.cos(Val) - 0.5*np.sin(Val) - 1)
        ExpNum2 = RampArray[Ind]*(0.5 + 0.5*np.sin(Val) + 0.5*np.cos(Val) - 1)
        Exp1.append(ExpNum1)
        Exp2.append(ExpNum2)

    ax[DirNo].plot(Step, Disp2, linestyle='-', color=colors[1], linewidth=2, label='Disp. ($u_y$)')
    ax[DirNo].plot(Step, Disp3, linestyle='-', color=colors[2], linewidth=2, label='Disp. ($u_z$)')

    ax[DirNo].plot(Step, Exp1, linestyle='--', color=colors[1], linewidth=2)
    ax[DirNo].plot(Step, Exp2, linestyle='--', color=colors[2], linewidth=2)
    ax2.plot(Step, DetF, linestyle=style[DirNo], color=colors[DirNo], linewidth=2, label=r'{0}'.format(LegTitle[DirNo]) )
    DirNo += 1

for i in [0,1,2,3]:
    ax[i].set_xlim((0, 20))
    ax[i].set_ylim((-0.4, 0.1))

    ax[i].set_xlabel("Simulation Steps", fontdict=font)
    ax[i].set_ylabel("Displacement", fontdict=font)
    ax[i].set_title(Title[i], fontdict=font)
    ax[i].legend(loc="best", frameon=True, prop=FontLegend)


plt.xlabel("Simulation Steps", fontdict=font)
plt.ylabel(r"$J=det \mathbf{F}$", fontdict=font)
plt.legend(loc="best", frameon=True, prop=FontLegend)

plt.show()
# # Save in order of last changed image
# plt.savefig("../Images/FormComp_J_N_15_Nu_049.pdf", transparent=True, bbox_inches='tight')
# plt.close()
# plt.savefig("../Images/FormComp_Disp_N_15_Nu_049.pdf", transparent=True, bbox_inches='tight')
# plt.close()
