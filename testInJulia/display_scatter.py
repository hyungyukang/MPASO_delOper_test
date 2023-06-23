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
#initFileName = 'x1.2562.grid.nc'
initFileName = 'x1.2562.grid.nc'
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

#output_init = outDS.variables['normalVelocity']
#output_init = outDS.variables['init']
#output_var = outDS.variables['divergence']
#output_init = outDS.variables['uEdge_anl']
#output_var = outDS.variables['uEdge']
#output_init = outDS.variables['u_anl']
#output_init = outDS.variables['divergence']
output_init = outDS.variables['uLapAnl']
#output_init = outDS.variables['vAnl']
#output_init = outDS.variables['relativeVorticityEdge']
#output_initAnl = outDS.variables['relativeVorticityAnlEdge']
#output_var = outDS.variables['divergenceVert']
#output_var = outDS.variables['divergenceEdge']
#output_var = outDS.variables['divergenceSL']
#output_var = outDS.variables['divergence']
#output_var = outDS.variables['divergence']
output_var = outDS.variables['u']
#====================================================#

#k=0
#iTime = 1
#figdpi = 200
#varNames = ['init']
#nVars = len(varNames)
#loc = ['cell']

#######################

figdpi = 200
fig = plt.figure(figsize=(11,11))
ax = plt.subplot(2,1,1)
#im = plt.scatter(lonCell,latCell,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
#im = plt.scatter(lonCell,latCell,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.1,vmax=0.1)
#im = plt.scatter(lonCell,latCell,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.005,vmax=0.005)
#im = plt.scatter(lonCell,latCell,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr)
#im = plt.scatter(lonVertex,latVertex,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.005,vmax=0.005)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.100,vmax=0.100)
im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.005,vmax=0.005)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr)

#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.jet)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-1.0,vmax=1.0)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.1,vmax=0.1)#,vmin=-10,vmax=10)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-10,vmax=10)
#im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
#im = plt.scatter(lonEdge,latEdge,c=output,s=8.0,marker='o',cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)

#xtick = np.arange(-180,210,60)
#ytick = np.arange(-90,120,30)
#ax.set_xlim([-180,180])
#ax.set_xticks(xtick)
#ax.set_ylim([-75,75])
ax.set_ylim([-90,90])
#ax.set_yticks(ytick)

ax = plt.subplot(2,1,2)
#im = plt.scatter(lonCell,latCell,c=output_init,s=13.0,marker='o',cmap=plt.cm.bwr)
#im = plt.scatter(lonCell,latCell,c=output_init,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
#im = plt.scatter(lonCell,latCell,c=output_init,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.005,vmax=0.005)
#im = plt.scatter(lonEdge,latEdge,c=output_init,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-10,vmax=10)
#im = plt.scatter(lonEdge,latEdge,c=output_init-output_initAnl,s=13.0,marker='o',cmap=plt.cm.bwr)
im = plt.scatter(lonEdge,latEdge,c=output_init,s=13.0,marker='o',cmap=plt.cm.bwr)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
#ax.set_ylim([-75,75])
ax.set_ylim([-90,90])
# Add the color bar
plt.colorbar(im,cax=cax)
plt.savefig('fig_scatter.png',bbox_inches='tight')
plt.close()
