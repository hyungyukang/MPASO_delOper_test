import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr



plt.rcParams.update({'font.size': 14})
plt.rc('legend',fontsize=12.5)

#====================================================#

runDir = './'
#initFileName = '../../meshes/oQU240.nc'
#initFileName = 'x1.2562.grid.nc'
#initFileName = 'x1.10242.grid.nc'
initFileName = 'x1.40962.grid.nc'
outputFileName = 'output.nc'
#outputFileName = 'output_time.nc'

initDS = xr.open_dataset(runDir+initFileName)
outDS = xr.open_dataset(runDir+outputFileName)

pi180 = np.pi / 180.0

nCells = initDS.dims['nCells']
nEdges = initDS.dims['nEdges']

N = int(np.sqrt(nCells))
lonCell = initDS.variables['lonCell']/pi180
latCell = initDS.variables['latCell']/pi180
lonEdge = initDS.variables['lonEdge']/pi180
latEdge = initDS.variables['latEdge']/pi180
lonVertex = initDS.variables['lonVertex']/pi180
latVertex = initDS.variables['latVertex']/pi180

dcEdge = initDS.variables['dcEdge']
dvEdge = initDS.variables['dvEdge']
#output_init = outDS.variables['normalVelocity']
#output_init = outDS.variables['init']
#output_var = outDS.variables['divergence']
#output_init = outDS.variables['uEdge_anl']
#output_var = outDS.variables['uEdge']
#output_init = outDS.variables['u_anl']
#output_init = outDS.variables['divergence']
#output_init1 = outDS.variables['vLapAnl']

#output_init1 = outDS.variables['relativeVorticityAnlEdge']
#output_init1 = outDS.variables['uLapAnl']
output_init1 = outDS.variables['normalVelocityAnl']

#output_var = outDS.variables['divergenceVert']
#output_var = outDS.variables['divergenceEdge']
#output_var = outDS.variables['divergenceSL']

#output_var1 = outDS.variables['relativeVorticityEdge']
#output_var1 = outDS.variables['u']
output_var1 = outDS.variables['del2u']

#output_var2 = outDS.variables['uLapAnl']
#====================================================#

#k=0
#iTime = 1
#figdpi = 200
#varNames = ['init']
#nVars = len(varNames)
#loc = ['cell']

#xmin = 0.0
#xmax = 60.0
#ymin =   30
#ymax =   60

xmin = 0.0
xmax = 360.0
ymin =   -90
ymax =    90


Umin = -2000.0
Umax =  2000.0
Vmin = -2000.0
Vmax =  2000.0
Emin = -200.0
Emax =  200.0
#######################

fig = plt.figure(figsize=(9,12))
axs = fig.subplots(3,1)


#------
#ax = plt.subplot(4,1,1)
#ax = axs[0,0]
ax = axs[0]
im = ax.scatter(lonEdge,latEdge,c=output_var1,s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Umin,vmax=Umax)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
#------
#ax = plt.subplot(4,2,1)
#ax = axs[0,1]
ax = axs[1]
im = ax.scatter(lonEdge,latEdge,c=output_init1,s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Vmin,vmax=Vmax)
#im = ax.scatter(lonEdge,latEdge,c=dvEdge,s=10.0,marker='o')#,cmap=plt.cm.bwr)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
#------
#------
#------
#ax = axs[1,0]
#im = ax.scatter(lonEdge,latEdge,c=output_init1,s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Umin,vmax=Umax)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='2%', pad=0.10)
#plt.colorbar(im,cax=cax)
#ax.set_xlim([xmin,xmax])
#ax.set_ylim([ymin,ymax])
##------
#ax = axs[1,1]
#im = ax.scatter(lonEdge,latEdge,c=output_init2,s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Vmin,vmax=Vmax)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='2%', pad=0.10)
#ax.set_xlim([xmin,xmax])
#ax.set_ylim([ymin,ymax])
#plt.colorbar(im,cax=cax)
#------
#------
#------
#ax = axs[2,0]
ax = axs[2]
im = ax.scatter(lonEdge,latEdge,c=(output_var1-output_init1),s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Emin,vmax=Emax)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
#------
#ax = axs[2,1]
#im = ax.scatter(lonEdge,latEdge,c=(output_var2-output_init2),s=10.0,marker='o',cmap=plt.cm.bwr,vmin=Umin,vmax=Umax)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='2%', pad=0.10)
#plt.colorbar(im,cax=cax)
#ax.set_xlim([xmin,xmax])
#ax.set_ylim([ymin,ymax])


#print(output_var2/output_init2)

###############
#figdpi = 200
plt.savefig('fig_scatter.png',bbox_inches='tight')
plt.close()
