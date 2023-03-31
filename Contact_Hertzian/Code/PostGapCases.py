import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

font = {'family': 'serif', 'weight': 'normal', 'size': 18}
FontLegend = {'family': 'serif', 'weight': 'normal', 'size': 15}
plt.rc('font', **font)
colors= [ 'black', 'red', "royalblue", "lightblue", ]

# Read in text file
DirStart = '../Result/2DPenalty/'
DirList = ['Pen_1e+01', 'Pen_1e+02', 'Pen_1e+03', 'Pen_1e+04' ]
Title = ["$k_{pen} = 1e+01$","$k_{pen} = 1e+02$","$k_{pen} = 1e+03$","$k_{pen} = 1e+04$" ]
StepNumber = [10]
TotSteps = 10
Points = np.linspace(0, 1.0, 121)         # Points along the profile

fig, ax = plt.subplots(1,4, figsize=(18,5), sharey=True)
# fig2, ax2 = plt.subplots(figsize=(5,5))

DirNo = 0
for Dir in DirList:
    df1 = pd.read_table(DirStart + Dir + '/PostProc.txt', delimiter=' ', header=None)
    df2 = pd.read_table(DirStart + Dir + '/DispProf.txt', delimiter=' ', header=None)
    df3 = pd.read_table(DirStart + Dir + '/GapProf.txt', delimiter=' ', header=None)

    df1_convert = df1.to_numpy()
    df2_convert = df2.to_numpy()
    df3_convert = df3.to_numpy()

    Step = df1_convert[:,0]
    Depth = df1_convert[:,1]

    ProfY = -Depth[8] + (Points-0.5)**2
    for LineNo in StepNumber:
        ax[DirNo].plot(Points, df2_convert[LineNo,:], 'k-', linewidth=2, label='Displacement ($u_y$)')
        ax[DirNo].plot(Points, df3_convert[LineNo,:], 'b-', linewidth=2, label='Gap')
        ax[DirNo].plot(Points, ProfY, 'r--', label='Indenter')
    # for LineNo in [10]:
    #     ax2.plot(Points, df4_convert[LineNo,:], 'k-', linewidth=2, color=colors[DirNo])
    DirNo += 1

for i in [0,1,2,3]:
    ax[i].set_xlim((0, 1))
    ax[i].set_xlabel("Profile", fontdict=font)
    ax[i].set_ylabel("Displacement", fontdict=font)
    ax[i].set_title(Title[i], fontdict=font)
ax[3].legend(loc="best", frameon=True, prop=FontLegend)
plt.savefig("../Images/PenaltyVarying.pdf", transparent=True, bbox_inches='tight')
# plt.close()
plt.show()
