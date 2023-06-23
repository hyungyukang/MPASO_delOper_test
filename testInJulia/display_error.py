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

nw = 5


coord_480 = np.zeros(5)
coord_480[0] = 5
coord_480[1] = 10
coord_480[2] = 20
coord_480[3] = 40
coord_480[4] = 80

error_480 = np.zeros(5)
error_240 = np.zeros(5)
ls_2 = np.zeros(5)
ls_1 = np.zeros(5)


error_480[0] = 0.011000579511857
error_480[1] = 0.0389727989634848
error_480[2] = 0.142099077015883
error_480[3] = 0.450264770219896
error_480[4] = np.nan

error_240[0] = 0.003028011
error_240[1] = 0.009930329
error_240[2] = 0.037536507
error_240[3] = 0.139513698
error_240[4] = 0.445805049

ls_2[0] = error_240[0]+0.0005
ls_1[0] = error_240[0]-0.0005
for i in range(nw-1):
    ls_2[i+1] = ls_2[0]*((coord_480[i+1]/coord_480[0])**2)
    ls_1[i+1] = ls_1[0]*((coord_480[i+1]/coord_480[0])   )
    

fig = plt.figure(figsize=(11,10))

ax = plt.subplot(1,1,1)
im = plt.plot(coord_480,error_480,marker='o',markersize=8.0,lw=2.5,color='r',label='Error-QU480')
im = plt.plot(coord_480,error_240,marker='o',markersize=8.0,lw=2.5,color='b',label='Error-QU240')
im = plt.plot(coord_480,ls_2,color='grey',lw=1.5,ls='--',label='2nd-order reference')
im = plt.plot(coord_480,ls_1,color='grey',lw=1.5,ls='-.',label='2nd-order reference')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Spherical harmonic degree')
ax.set_ylabel('Power')

ax.set_xticks([5,10,20,40,80])
ax.set_xticklabels([5,10,20,40,80])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.grid(True,which='both')
ax.yaxis.grid(True,which='both')
ax.legend(loc='upper left')
#ax.set_xlim([0,41])
ax.set_ylim([1e-3,1e-0])


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

