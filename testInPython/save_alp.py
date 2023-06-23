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
from pyshtools.legendre import PlmON,PlmIndex

#====================================================#
# Setup grid and total wavenumber

DX = 240 # km

# Half-wavelength spatial resolution (km)
hWavelength = 2.0 * np.pi * 6371.0 / 2.0

if (   DX == 480 ):
         # dx = 480 km ~ 4deg
         ngrid = "2562"
         WN = int(hWavelength/DX)-1
         KM = "480.0"
elif ( DX == 240 ):
         # dx = 240 km ~ 2deg
         ngrid = "10242"
         WN = int(hWavelength/DX)-1
         KM = "240.0"
elif ( DX == 120 ):
         # dx = 120 km ~ 1deg
         ngrid = "40962"
         WN = int(hWavelength/DX)-1
         KM = "120.0"
elif ( DX == 60 ):
         # dx = 60 km ~ 1deg
         ngrid = "163842"
         WN = int(hWavelength/DX)-1
         KM = "060.0"
elif ( DX == 30 ):
         # dx = 30 km ~ 1deg
         ngrid = "655362"
         WN = int(hWavelength/DX)-1
         KM = "030.0"
elif ( DX == 15 ):
         # dx = 15 km ~ 1deg
         ngrid = "2621442"
         WN = int(hWavelength/DX)-1
         KM = "015.0"
elif ( DX == 7.5 ):
         # dx = 15 km ~ 1deg
         ngrid = "10485762"
         WN = int(hWavelength/DX)-1
         KM = "007.5"

WN1 = WN + 1
WNH = int(WN/2.0)

#====================================================#
# Read files

#runDir = '../../meshes/'
#initFileName = "x1."+ngrid+".grid.nc"
#outName = "./ALP_3modes_QU"+KM+"km_Cell.nc"

runDir = './'
initFileName = "oQU240km.nc"
outName = "./ALP_3modes_oQU240km_Cell.nc"


#outputFileName = 'output.nc'
#outputFileName = 'output_time.nc'

initDS = xr.open_dataset(runDir+initFileName) # For mesh info
#outDS = xr.open_dataset(runDir+outputFileName) # For output file read

pi180 = np.pi / 180.0

nCells = initDS.dims['nCells']
lonCell = initDS.variables['lonCell']/pi180 # RAD -> DEG
latCell = initDS.variables['latCell']/pi180 # RAD -> DEG


#output = outDS.variables['init']
#output = outDS.variables['psi_n0']
#====================================================#

# Get SPH coefficient on the MPAS grid
#cilm, chi2 = SHExpandLSQ(output,latCell,lonCell,30,4,1)
#
## Get power spectrum
#shcoeffs = pysh.SHCoeffs.from_array(cilm)
#psectrum = shcoeffs.spectrum()
#
#print(psectrum[:])
#
## Plot power spectrum
#fig, ax = shcoeffs.plot_spectrum(xscale='lin',yscale='log',show=False,marker='o',color='r')
#ax.legend(loc='upper left')
#ax.set_xlim([0,31])
#plt.savefig('fig_spectrum.png',bbox_inches='tight')
#plt.close()
#

totNwave = int((WN+1)*(WN+2)/2)

alp = np.zeros(totNwave)


nFields = 3
alp_save = np.zeros([nFields,nCells])

for iCell in range(nCells):
    z = np.sin(latCell[iCell] * pi180)
    lon = lonCell[iCell] * pi180

    alp = PlmON(WN,z)

    #----------------
    m = 0 ; n = WN
    index = PlmIndex(n,m)
    alp_save[0,iCell] = np.real(np.exp(m*1j*lon)*alp[index])
   
    #----------------
    m = int(WN/4.0) ; n = int(WN/2.0)
    #m = 0 ; n = WN*2
    index = PlmIndex(n,m)
    alp_save[1,iCell] = np.real(np.exp(m*1j*lon)*alp[index])

    #----------------
    m = WN ; n = WN
    index = PlmIndex(n,m)
    alp_save[2,iCell] = np.real(np.exp(m*1j*lon)*alp[index])

#=====================================================#

# Save init_alps
ds = nc.Dataset(outName, 'w', format='NETCDF4')
nCell = ds.createDimension('nCells',nCells)
nMode = ds.createDimension('nMode',nFields)
nWN = ds.createDimension('WN',WN)
lats = ds.createVariable('latCell','f8',('nCells',))
lats[:] = latCell[:]
lons = ds.createVariable('lonCell','f8',('nCells',))
lons[:] = lonCell[:]
Fields = ds.createVariable('alp','f8',('nMode','nCells',))
Fields[:,:] = alp_save[:,:]

ds.close()

#=====================================================#

cilm, chi2 = SHExpandLSQ(alp_save[2,:],latCell,lonCell,WN,4,1)
shcoeffs = pysh.SHCoeffs.from_array(cilm)
psectrum = shcoeffs.spectrum()
print(psectrum[:])
