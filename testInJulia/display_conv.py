import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 14})
plt.rc('legend',fontsize=12.5)

file1 = open('data_conv_lap_div_vor.txt','r')

nCol = 2
nRow = 7

i = 0
data = np.zeros([nRow,nCol])
for line in file1.readlines():
    temp = line.split('\t')
    temp_array = []
    for j in range(nCol):
        temp_array.append(temp[j])
    data[i,:] = temp_array
    i = i+1
file1.close

print('=== Display data ===')
print(data)

ls1 = np.zeros(nRow)
ls2 = np.zeros(nRow)
d1 = np.zeros(nRow)
d2 = np.zeros(nRow)

dx = np.zeros(nRow)
for i in range(nRow):
    dx[i] = data[i,0]
    d1[i] = data[i,1]

# --- Linear scaling
ls1[0] = d1[0]*0.95
ls2[0] = d1[0]*1.05
for i in range(nRow-1):
    ls1[i+1] = ls1[0]  / ((dx[0]/dx[i+1])**2)
    ls2[i+1] = ls2[0]  / ((dx[0]/dx[i+1])   )

fig,axs = plt.subplots(1,1,figsize=(5.5,5.0),constrained_layout=True)
ax = axs

ax.plot(dx,d1,marker='o',color='r',ls='-',label='div(grad(scalar))',lw=2.5,markersize=7)
ax.plot(dx,ls1,color='grey',label='2nd-order reference',lw=1.5,ls='--')
ax.plot(dx,ls2,color='grey',label='1st-order reference',lw=1.5,ls='dashdot')

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(([480,240,120,60,30,15,7.5]))
ax.set_xticklabels(([480,240,120,60,30,15,7.5]))
ax.invert_xaxis()
#ax.set_yticks(([0.1,0.2,0.3,0.4,0.5,0.8,1.0]))
#ax.set_yticklabels(([0.1,0.2,0.3,0.4,0.5,0.8,1.0]))
ax.tick_params(axis='x',which='minor',bottom=False)
ax.tick_params(axis='x',which='major',bottom=True)
#ax.set_ylim([0.1,1.1])

ax.grid(True,axis='y',which='both',ls='--',lw=0.5)
ax.grid(True,axis='x',which='major',ls='--',lw=0.5)

#ax.legend(loc='upper center',labelspacing=0.25,handlelength=3,ncol=2)
ax.legend(loc='upper right')#,labelspacing=0.25,handlelength=3)
#
ax.set_xlabel("dx (km)")
ax.set_ylabel("L2 norm error")

plt.savefig('fig_conv.png',bbox_inches='tight',dpi=200)
