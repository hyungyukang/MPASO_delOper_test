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
from glob import glob



plt.rcParams.update({'font.size': 14})
plt.rc('legend',fontsize=12.5)

#====================================================#

runDir = './'
#initFileName = 'x1.2562.grid.nc'
#initFileName = 'x1.2562.grid.nc'
outputFileName = 'y1.*.nc' 
#outputFileName = 'output_time.nc'

#initDS = xr.open_mfdataset(runDir+initFileName)
#outDS = xr.open_mfdataset(runDir+glob(outputFileName))
outDS = xr.open_mfdataset('y1.*.nc',concat_dim='Time',combine='nested')
print(outDS)

#pi180 = np.pi / 180.0
#
#nCells = initDS.dims['nCells']
#nEdges = initDS.dims['nEdges']
#
#N = int(np.sqrt(nCells))
#lonCell = initDS.variables['lonCell']/pi180
#latCell = initDS.variables['latCell']/pi180
#lonEdge = initDS.variables['lonEdge']/pi180
#latEdge = initDS.variables['latEdge']/pi180
#
##output = outDS.variables['normalVelocity']
##output_init = outDS.variables['init']
##output_var = outDS.variables['divergence']
##output_init = outDS.variables['uEdge_anl']
##output_var = outDS.variables['uEdge']
#output_init = outDS.variables['alp']
#output_var = outDS.variables['alp']
##====================================================#
#
##k=0
##iTime = 1
##figdpi = 200
##varNames = ['init']
##nVars = len(varNames)
##loc = ['cell']
#
########################
#
#figdpi = 200
#fig = plt.figure(figsize=(11,12))
#ax = plt.subplot(2,1,1)
#im = plt.scatter(lonCell,latCell,c=output_var[1,:],s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-1.0,vmax=1.0)
##im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
##im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-1.0,vmax=1.0)
##im = plt.scatter(lonEdge,latEdge,c=output_var[2,:],s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-20,vmax=20)
##im = plt.scatter(lonEdge,latEdge,c=output_var,s=13.0,marker='o',cmap=plt.cm.bwr)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='2%', pad=0.10)
#plt.colorbar(im,cax=cax)
##im = plt.scatter(lonEdge,latEdge,c=output,s=8.0,marker='o',cmap=plt.cm.bwr,vmin=-0.6,vmax=0.6)
#
#xtick = np.arange(-180,210,60)
#ytick = np.arange(-90,120,30)
#ax.set_xlim([-180,180])
#ax.set_xticks(xtick)
#ax.set_ylim([-90,90])
#ax.set_yticks(ytick)
#
#ax = plt.subplot(2,1,2)
#im = plt.scatter(lonCell,latCell,c=output_init[2,:],s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-1.0,vmax=1.0)
##im = plt.scatter(lonEdge,latEdge,c=output_init[0,:],s=13.0,marker='o',cmap=plt.cm.bwr,vmin=-20,vmax=20)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes('right', size='2%', pad=0.10)
## Add the color bar
#xtick = np.arange(-180,210,60)
#ytick = np.arange(-90,120,30)
#ax.set_xlim([-180,180])
#ax.set_xticks(xtick)
#ax.set_ylim([-90,90])
#ax.set_yticks(ytick)
#
#plt.colorbar(im,cax=cax)
#plt.savefig('fig_scatter.png',bbox_inches='tight')
#plt.close()
