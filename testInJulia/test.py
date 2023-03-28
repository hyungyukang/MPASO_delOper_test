import xarray
import numpy
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from mpas_tools.viz import mesh_to_triangles
#import mpas_tools

dsMesh = xarray.open_dataset('output.nc')
dsTris = mesh_to_triangles(dsMesh, periodicCopy=True)

sst = dsMesh.init.isel(Time=0, nVertLevels=0).values

sstTri = sst[dsTris.triCellIndices]

nTriangles = dsTris.sizes['nTriangles']
tris = numpy.arange(3*nTriangles).reshape(nTriangles, 3)

lonNode = numpy.rad2deg(dsTris.lonNode.values).ravel()
latNode = numpy.rad2deg(dsTris.latNode.values).ravel()
#sstNode = sstNode.values.ravel()

triangles = Triangulation(lonNode, latNode, tris)

#plt.figure(1)
#plt.tripcolor(triangles, sstNode, shading='gouraud')
#plt.xlim([0., 360.])
#plt.ylim([-90., 90.])

plt.figure(1)
plt.tripcolor(triangles, sstTri, shading='flat')
plt.xlim([0., 360.])
plt.ylim([-90., 90.])

plt.show()
