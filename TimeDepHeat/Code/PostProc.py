import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

font = {'family': 'serif', 'weight': 'normal', 'size': 18}
FontLegend = {'family': 'serif', 'weight': 'normal', 'size': 15}
plt.rc('font', **font)
colors= [ 'black', 'red', "blue","violet", "lightblue", ]
style = ['-', '--', '-.', '-', ':']
# Read in text file
DirStart = '../Results/Interpolate'
DirList = ["/S_50_N_10/", "/S_50_N_20/","/S_50_N_40/", "/S_50_N_50/", "/S_50_N_100/"]
Title = [r"$N = 10$", r"$N=20$", "S_50_N_40", r"$N=50$", r"$N=100"]

DirLen = len(DirList)
# Title = [r"$\nu = 0.4$", r"$\nu = 0.45$", r"$\nu = 0.49$"]
LegTitle = [r"$10$", r"$20$",r"$40$",r"$50$", r"$100$"]
# LegTitle = ["Compressible", "Penalty", r"Penalty Weak Form: $\mathbf{u}$", r"Weak Form: $\mathbf{u}$, $p$"]

fig, ax = plt.subplots(1,1, figsize=(6,5), sharey=True)
fig2, ax2 = plt.subplots(1,1, figsize=(6,5), sharey=True)

# ax2.set_xlim((0,20))
# plt.axhline(1.0,color='black', linestyle='--', alpha=0.5)
DirNo = 0
for Dir in DirList:
    # Step, Scale of expression, Displacement x, y, and z
    df1 = pd.read_table(DirStart + Dir + '/PostProc.txt', delimiter=' ', header=None)

    df1_convert = df1.to_numpy()
    Step = df1_convert[:,0]
    Time = df1_convert[:,1]
    Error = df1_convert[:,2]

    ax.plot(Step, Error, linestyle=style[DirNo],color=colors[DirNo], linewidth=2, label=LegTitle[DirNo])
    # ax[DirNo].plot(Step, Disp3, linestyle='-', color=colors[2], linewidth=2, label='Disp. ($u_z$)')

    # ax[DirNo].plot(Step, Exp1, linestyle='--', color=colors[1], linewidth=2)
    # ax[DirNo].plot(Step, Exp2, linestyle='--', color=colors[2], linewidth=2)
    ax2.plot(Time, Error, linestyle=style[DirNo], color=colors[DirNo], linewidth=2, label=LegTitle[DirNo])
    DirNo += 1

ax2.set_xlabel("Simulation Time", fontdict=font)
ax2.set_ylabel("Error", fontdict=font)
ax2.legend(loc="best", title="N", frameon=True, prop=FontLegend)
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.savefig("../Images/CompIntSimTime.pdf", transparent=True, bbox_inches='tight')
plt.close()

ax.set_xlabel("Simulation Steps", fontdict=font)
ax.set_ylabel("Error", fontdict=font)
ax.legend(loc="best", title="N", frameon=True, prop=FontLegend)
# ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig("../Images/CompIntSimSteps.pdf", transparent=True, bbox_inches='tight')
plt.close()


plt.close()
