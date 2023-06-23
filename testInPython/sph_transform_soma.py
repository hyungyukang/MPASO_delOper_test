import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr
import pyshtools as pysh
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyshtools.expand import SHExpandLSQ

plt.rcParams.update({'font.size': 16})
plt.rc('legend',fontsize=14)

#====================================================#
# Read files

runDir = './'

globFileName = 'x1.40962.grid.nc'
outputFileName = 'output.0001-01-01_00.00.00_ori_nu2.0e10.nc'

initDS = xr.open_dataset(runDir+globFileName) # For mesh info
outDS = xr.open_dataset(runDir+outputFileName) # For output file read

pi180 = np.pi / 180.0

nCellsGlob = initDS.dims['nCells']
lonCellGlob = initDS.variables['lonCell']/pi180 # RAD -> DEG
latCellGlob = initDS.variables['latCell']/pi180 # RAD -> DEG

nCellsOut =  outDS.dims['nCells']
lonCellOut = outDS.variables['lonCell']/pi180 # RAD -> DEG
latCellOut = outDS.variables['latCell']/pi180 # RAD -> DEG

nTime  = outDS.dims['Time']

#####

tt = nTime - 1
zz = 0

#output = outDS.variables['divergence'][tt,:,zz]
output = outDS.variables['temperature'][tt,:,zz]
#output = outDS.variables['kineticEnergyCell'][tt,:,zz]

#####

# To rotate longitude...
for iCell in range(nCellsOut):
    if ( lonCellOut[iCell] > 180.0 ):
         lonCellOut[iCell] = lonCellOut[iCell] - 360.0
#    lonCellOut[iCell] = lonCellOut[iCell]+180


# Lon & Lat extension for SOMA simulation domain

minlat = np.min(latCellOut)
maxlat = np.max(latCellOut)
minlon = np.min(lonCellOut)
maxlon = np.max(lonCellOut)

meanLat = 0.5 * (minlat + maxlat)

latCellOut = latCellOut - meanLat
lonCellOut = lonCellOut + (-minlon)

minlat = np.min(latCellOut)
maxlat = np.max(latCellOut)
minlon = np.min(lonCellOut)
maxlon = np.max(lonCellOut)

#ratioLat =  89.9  / maxlat
#ratioLon =  359.9 / maxlon
ratioLat =  69.9  / maxlat
#ratioLon =  329.9 / maxlon

dx = 16.0 * ratioLat
print(dx)

#hWavelength = 2.0 * np.pi * 6371.0 / 2.0
#WN = int(hWavelength/dx)-1
#WN = int(WN/2)
WN = 120
#
lonCellOut = lonCellOut * ratioLat + 90.0
latCellOut = latCellOut * ratioLat

# Sort
#sortArgOut = lonCellOut.argsort()
#lonCellOut = lonCellOut[sortArgOut]
#latCellOut = latCellOut[sortArgOut]
#
#sortArgGlob = lonCellGlob.argsort()
#lonCellGlob = lonCellGlob[sortArgGlob]
#latCellGlob = latCellGlob[sortArgGlob]

minLonOut = np.min(lonCellOut)#-1.0
maxLonOut = np.max(lonCellOut)#+1.0
minLatOut = np.min(latCellOut)#-1.0
maxLatOut = np.max(latCellOut)#+1.0

print(minLonOut,maxLonOut)
print(minLatOut,maxLatOut)

#############################################

nGlob = nCellsGlob

outGlob = np.zeros(nCellsGlob)

for iCell in range(nCellsGlob):
    if ( lonCellGlob[iCell] >= minLonOut and lonCellGlob[iCell] <= maxLonOut and
         latCellGlob[iCell] >= minLatOut and latCellGlob[iCell] <= maxLatOut ):
         nGlob = nGlob - 1 

nCellsMerge = nGlob + nCellsOut
lonCell = np.zeros(nCellsMerge)
latCell = np.zeros(nCellsMerge)
outCell = np.zeros(nCellsMerge)
print(nCellsMerge)
outCell[:] = -999.99

i = 0
for iCell in range(nCellsGlob):
    if ( lonCellGlob[iCell] < minLonOut or lonCellGlob[iCell] > maxLonOut or 
         latCellGlob[iCell] < minLatOut or latCellGlob[iCell] > maxLatOut ):
         lonCell[i] = lonCellGlob[iCell] 
         latCell[i] = latCellGlob[iCell] 
         outCell[i] = 0.0
         i = i + 1

lonCell[nGlob:] = lonCellOut[:]
latCell[nGlob:] = latCellOut[:]
outCell[nGlob:] = output[:]

outCellList = outCell.tolist()
lonCellList = lonCell.tolist()
latCellList = latCell.tolist()

print(len(outCellList))

############

#for iCell in range(nCellsMerge):
#    if ( outCell[iCell] == -999.99):
#         outCell[iCell] = 0.0
#         lonCell[iCell] = 

# 
#for iCell 

#
##output = outDS.variables['init']
#
#
##output = outDS.variables['kineticEnergyCell']
#output = outDS.variables['divergence']
##WN = outDS.dims['WN']
##====================================================#
#
## Lon & Lat extension for SOMA simulation domain
#
#minlat = np.min(latCell)
#maxlat = np.max(latCell)
#minlon = np.min(lonCell)
#maxlon = np.max(lonCell)
#
#meanLat = 0.5 * (minlat + maxlat)
#
#latCell = latCell - meanLat
#lonCell = lonCell + (-minlon)
#
#minlat = np.min(latCell)
#maxlat = np.max(latCell)
#minlon = np.min(lonCell)
#maxlon = np.max(lonCell)
#
##ratioLat =  89.9  / maxlat
##ratioLon =  359.9 / maxlon
#ratioLat =  69.9  / maxlat
#ratioLon =  329.9 / maxlon
#
#dx = 16.0 * ratioLat
#print(dx)
#
#hWavelength = 2.0 * np.pi * 6371.0 / 2.0
#WN = int(hWavelength/dx)-1
#WN = int(WN/2)
#WN = 40
#
##lonCell = lonCell * ratioLon
##latCell = latCell * ratioLat
#
#minlat = np.min(latCell)
#maxlat = np.max(latCell)
#minlon = np.min(lonCell)
#maxlon = np.max(lonCell)
#
#print(minlat,maxlat)
#print(minlon,maxlon)
#
#print(WN)
#print(ratioLat,ratioLat)
##print(ratioLonM,ratioLonP)
#print('val output : ',output.values[tt,0:20,zz])
#
##====================================================#
#
##Add grid patch to expand current SOMA domain into a whole sphere
#
## --- Sort first
#sortArg = lonCell.argsort()
#lonCellSort = lonCell[sortArg]
#latCellSort = latCell[sortArg]
#outputSort = output[tt,sortArg,zz]
#
#

##====================================================#

fig = plt.figure(figsize=(11,5))
ax = plt.subplot(1,1,1)
im = ax.scatter(lonCell,latCell,c=outCell,s=1.0,marker='o',cmap=plt.cm.jet)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.10)
plt.colorbar(im,cax=cax)
#
#

##====================================================#

# Get SPH coefficient on the MPAS grid
cilm, chi2 = SHExpandLSQ(outCell,latCell,lonCell,WN,4,1)
 
# Get power spectrum
shcoeffs = pysh.SHCoeffs.from_array(cilm)
psectrum = shcoeffs.spectrum()
print(psectrum[:])

# Plot power spectrum
#pysh.utils.figstyle(rel_width=0.25)
#pysh.utils.figstyle(aspect_ratio=2/1.1,screen_dpi=10,rel_width=0.3)
#pysh.utils.figstyle(figsize=[2,1],screen_dpi=50,)
#pysh.utils.figstyle(figsize=[2,1])
#fig, ax = shcoeffs.plot_spectrum(xscale='lin',yscale='log',show=False,marker='o',color='r')

########################
#fig = plt.figure(figsize=(11,12))
#ax = plt.subplot(1,1,1)
#im = plt.plot(psectrum,marker='o',color='r',label='Power per degree')
#ax.set_yscale('log')
#ax.set_xlabel('Spherical harmonic degree')
#ax.set_ylabel('Power')
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.xaxis.grid(True,which='major')
#ax.yaxis.grid(True,which='both')
#ax.legend(loc='upper left')
#ax.set_xlim([0,41])
#ax.set_ylim([1e-35,10e+1])

#
#################
## Get SPH coefficient on the MPAS grid
#cilm, chi2 = SHExpandLSQ(output[2,:],latCell,lonCell,40,4,1)
## Get power spectrum
#shcoeffs = pysh.SHCoeffs.from_array(cilm)
#psectrum = shcoeffs.spectrum()
#print(psectrum[:])
#
#ax = plt.subplot(2,1,2)
#im = plt.plot(psectrum,marker='o',color='r',label='Power per degree')
#ax.set_yscale('log')
#ax.set_xlabel('Spherical harmonic degree')
#ax.set_ylabel('Power')
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.xaxis.grid(True,which='major')
#ax.yaxis.grid(True,which='both')
#ax.legend(loc='upper left')
#ax.set_xlim([0,41])
#ax.set_ylim([1e-35,10e+1])
#
plt.savefig('fig_spectrum.png',bbox_inches='tight')
plt.close()
