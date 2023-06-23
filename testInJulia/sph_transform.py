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
from pyshtools.expand import SHExpandLSQ
from matplotlib.ticker import AutoMinorLocator

plt.rcParams.update({'font.size': 16})
plt.rc('legend',fontsize=14)

#====================================================#
# Read files

runDir = './'
initFileName = 'x1.2562.grid.nc'
#outputFileName = 'output.nc'
outputFileName = 'ALP_3modes_QU480.0km_Cell.nc'

initDS = xr.open_dataset(runDir+initFileName) # For mesh info
outDS = xr.open_dataset(runDir+outputFileName) # For output file read

pi180 = np.pi / 180.0

nCells = initDS.dims['nCells']

lonCell = initDS.variables['lonCell']/pi180 # RAD -> DEG
latCell = initDS.variables['latCell']/pi180 # RAD -> DEG

#output = outDS.variables['init']
output = outDS.variables['alp']
#WN = outDS.dims['WN']
#====================================================#

# Get SPH coefficient on the MPAS grid
cilm, chi2 = SHExpandLSQ(output[1,:],latCell,lonCell,40,4,1)
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

fig = plt.figure(figsize=(11,12))

ax = plt.subplot(2,1,1)
im = plt.plot(psectrum,marker='o',color='r',label='Power per degree')
ax.set_yscale('log')
ax.set_xlabel('Spherical harmonic degree')
ax.set_ylabel('Power')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.grid(True,which='major')
ax.yaxis.grid(True,which='both')
ax.legend(loc='upper left')
ax.set_xlim([0,41])
ax.set_ylim([1e-35,10e+1])


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

