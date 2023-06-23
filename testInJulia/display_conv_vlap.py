import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 14})
plt.rc('legend',fontsize=11.0)

file1 = open('data_conv_lap_div_vor.txt','r')

nCol = 6
nRow = 6

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
d3 = np.zeros(nRow)
d4 = np.zeros(nRow)
d5 = np.zeros(nRow)

dx = np.zeros(nRow)
for i in range(nRow):
    dx[i] = data[i,0]
    d1[i] = data[i,1]
    d2[i] = data[i,2]
    d3[i] = data[i,3]
    d4[i] = data[i,4]
    d5[i] = data[i,5]

# --- Linear scaling
ls1[0] = d5[0]*1.30
ls2[0] = d4[0]*1.10
for i in range(nRow-1):
    ls1[i+1] = ls1[0]  / ((dx[0]/dx[i+1])**2)
    ls2[i+1] = ls2[0]  / ((dx[0]/dx[i+1])   )


fig,axs = plt.subplots(1,1,figsize=(5.5,5.0),constrained_layout=True)
ax = axs

ax.plot(dx,d1,marker='o',color='k',ls='-',label='Del2 (default) ',lw=2.0,markersize=5)
ax.plot(dx,d4,marker='o',color='r',ls='-',label='Del2 (rVorIntp)',lw=2.0,markersize=5)
ax.plot(dx,d2,marker='o',color='g',ls='-',label='Div            ',lw=2.0,markersize=5)
ax.plot(dx,d3,marker='o',color='k',ls='--',label='rVor (default) ',lw=2.0,markersize=5)
ax.plot(dx,d5,marker='o',color='r',ls='--',label='rVor (Intp)    ',lw=2.0,markersize=5)

ax.plot(dx,ls1,color='grey',label='2nd-order Ref. ',lw=1.0,ls='--')
ax.plot(dx,ls2,color='grey',label='1st-order Ref. ',lw=1.2,ls=':')

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(([480,240,120,60,30,15]))
ax.set_xticklabels(([480,240,120,60,30,15]))
ax.invert_xaxis()
#ax.set_yticks(([0.1,0.2,0.3,0.4,0.5,0.8,1.0]))
#ax.set_yticklabels(([0.1,0.2,0.3,0.4,0.5,0.8,1.0]))
ax.tick_params(axis='x',which='minor',bottom=False)
ax.tick_params(axis='x',which='major',bottom=True)
#ax.set_ylim([0.1,1.1])

ax.grid(True,axis='y',which='both' ,ls='--',lw=0.5)
ax.grid(True,axis='x',which='major',ls='--',lw=0.5)

#ax.legend(loc='upper center',labelspacing=0.25,handlelength=3,ncol=2)
ax.legend(loc='lower left',labelspacing=0.30,handlelength=3.5)
#
ax.set_xlabel("dx (km)")
ax.set_ylabel("L2 norm error")

plt.savefig('fig_conv.png',bbox_inches='tight',dpi=200)
