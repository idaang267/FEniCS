import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

font = {'family': 'serif', 'weight': 'normal', 'size': 18}
FontLegend = {'family': 'serif', 'weight': 'normal', 'size': 15}
plt.rc('font', **font)
colors= [ 'black', 'red', "royalblue", "lightblue", ]

# Read in text file
DirStart = '../Result/'
DirList = ["plane_stress/S_50/", "plane_strain/S_50/"]
Title = ["Plane Stress", "Plane Strain"]
ErrMeanArr = []

DirNo = 0
fig, ax = plt.subplots(1,2, figsize=(12,5))
for Dir in DirList:

    # Step, Rho, Displacement, Analytical
    df1 = pd.read_table(DirStart + Dir + '/PostProc.txt', delimiter=' ', header=None)

    df1_convert = df1.to_numpy()
    Step = df1_convert[:,0]
    Rho = df1_convert[:,1]
    Disp = df1_convert[:,2]
    DispAny = df1_convert[:,3]
    ErrArr = []
    for (Count, Val) in enumerate(DispAny):
        if Disp[Count] != 0:
            PerErr = 100*abs(Disp[Count]-Val)/Disp[Count]
            ErrArr.append(PerErr)

    ErrMean = np.mean(ErrArr)
    ErrMeanArr.append(ErrMean)
    ax[DirNo].plot(Rho, Disp, 'r--', linewidth=2, label='Simulation')
    ax[DirNo].plot(Rho, DispAny, 'k-', linewidth=2, label='Analytical')

    DirNo += 1

for i in [0,1]:
    ax[i].set_xlim((0, 0.2))

    ax[i].set_xlabel("Body Force", fontdict=font)
    ax[i].set_ylabel("Displacement ($u_y$)", fontdict=font)
    ax[i].set_title(Title[i] + ": Error {0:.1f}%".format(ErrMeanArr[i]), fontdict=font)
    ax[i].legend(loc="best", frameon=True, prop=FontLegend)


plt.savefig("../Images/PlaneStrainVsStress.pdf", transparent=True, bbox_inches='tight')

plt.show()
