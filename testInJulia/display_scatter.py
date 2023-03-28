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
initFileName = 'x1.2562.grid.nc'
outputFileName = 'output.nc'

initDS = xr.open_dataset(runDir+initFileName)
outDS = xr.open_dataset(runDir+outputFileName)

pi180 = np.pi / 180.0

nCells = initDS.dims['nCells']
N = int(np.sqrt(nCells))
lonCell = initDS.variables['lonCell']/pi180
latCell = initDS.variables['latCell']/pi180

output = outDS.variables['init']
#====================================================#

k=0
iTime = 1
figdpi = 200
varNames = ['init']
nVars = len(varNames)
loc = ['cell']

#######################


figdpi = 200
fig = plt.figure(figsize=(11,6))
ax = plt.subplot(1,1,1)
im = plt.scatter(lonCell,latCell,c=output,s=13.0,marker='o',cmap=plt.cm.bwr)

xtick = np.arange(-180,210,60)
ytick = np.arange(-90,120,30)
ax.set_xlim([-180,180])
ax.set_xticks(xtick)
ax.set_ylim([-90,90])
ax.set_yticks(ytick)

# Add the color bar
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
plt.savefig('fig_scatter.png',bbox_inches='tight')
plt.close()
