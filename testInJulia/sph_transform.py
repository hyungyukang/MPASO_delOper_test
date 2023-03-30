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

plt.rcParams.update({'font.size': 14})
plt.rc('legend',fontsize=12.5)

#====================================================#
# Read files

runDir = './'
initFileName = 'x1.2562.grid.nc'
#outputFileName = 'output.nc'
outputFileName = 'output_time.nc'

initDS = xr.open_dataset(runDir+initFileName) # For mesh info
outDS = xr.open_dataset(runDir+outputFileName) # For output file read

pi180 = np.pi / 180.0

nCells = initDS.dims['nCells']
N = int(np.sqrt(nCells))
lonCell = initDS.variables['lonCell']/pi180 # RAD -> DEG
latCell = initDS.variables['latCell']/pi180 # RAD -> DEG

#output = outDS.variables['init']
output = outDS.variables['psi_n0']
#====================================================#

# Get SPH coefficient on the MPAS grid
cilm, chi2 = SHExpandLSQ(output,latCell,lonCell,30,4,1)

# Get power spectrum
shcoeffs = pysh.SHCoeffs.from_array(cilm)
psectrum = shcoeffs.spectrum()

print(psectrum[:])

# Plot power spectrum
fig, ax = shcoeffs.plot_spectrum(xscale='lin',yscale='log',show=False,marker='o',color='r')
ax.legend(loc='upper left')
ax.set_xlim([0,31])
plt.savefig('fig_spectrum.png',bbox_inches='tight')
plt.close()

