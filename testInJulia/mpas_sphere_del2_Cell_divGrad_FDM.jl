using NCDatasets
#using Plots
using PyPlot

###########################################################

mesh_file_name = "../../meshes/x1.2562.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.10242.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.40962.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.163842.grid.nc"   # Miniweather output file name
ds = Dataset(mesh_file_name,"r")

###########################################################

latCell           = ds["latCell"][:]
lonCell           = ds["lonCell"][:]
meshDensity       = ds["meshDensity"][:]
xCell             = ds["xCell"][:]
yCell             = ds["yCell"][:]
zCell             = ds["zCell"][:]
indexToCellID     = ds["indexToCellID"][:]
latEdge           = ds["latEdge"][:]
lonEdge           = ds["lonEdge"][:]
xEdge             = ds["xEdge"][:]
yEdge             = ds["yEdge"][:]
zEdge             = ds["zEdge"][:]
indexToEdgeID     = ds["indexToEdgeID"][:]
latVertex         = ds["latVertex"][:]
lonVertex         = ds["lonVertex"][:]
xVertex           = ds["xVertex"][:]
yVertex           = ds["yVertex"][:]
zVertex           = ds["zVertex"][:]
indexToVertexID   = ds["indexToVertexID"][:]
cellsOnEdge       = ds["cellsOnEdge"][:]
nEdgesOnCell      = ds["nEdgesOnCell"][:]
nEdgesOnEdge      = ds["nEdgesOnEdge"][:]
edgesOnCell       = ds["edgesOnCell"][:]
edgesOnEdge       = ds["edgesOnEdge"][:]
weightsOnEdge     = ds["weightsOnEdge"][:]
dvEdge            = ds["dvEdge"][:]
#dv1Edge           = ds["dv1Edge"][:]
#dv2Edge           = ds["dv2Edge"][:]
dcEdge            = ds["dcEdge"][:]
angleEdge         = ds["angleEdge"][:]
areaCell          = ds["areaCell"][:]
areaTriangle      = ds["areaTriangle"][:]
cellsOnCell       = ds["cellsOnCell"][:]
verticesOnCell    = ds["verticesOnCell"][:]
verticesOnEdge    = ds["verticesOnEdge"][:]
edgesOnVertex     = ds["edgesOnVertex"][:]
cellsOnVertex     = ds["cellsOnVertex"][:]
kiteAreasOnVertex = ds["kiteAreasOnVertex"][:]

nCells = ds.dim["nCells"]
nEdges = ds.dim["nEdges"]
nVertices = ds.dim["nVertices"]
maxEdges = ds.dim["maxEdges"]
maxEdges2 = ds.dim["maxEdges2"]
vertexDegree = ds.dim["vertexDegree"]


edgeSignOnCell = zeros(Float64, maxEdges,nCells)
 
for iCell in 1:nCells
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        
        if ( iCell == cellsOnEdge[1,iEdge])
           edgeSignOnCell[i,iCell] = -1.0
        else
           edgeSignOnCell[i,iCell] =  1.0
        end
    end
end 


###########################################################

# Initial field = Spherical harmonics (Y_n^m)
# del2 Analytic solution = -n(n+1)*Y_n^m

init = zeros(Float64, nCells)
anl = zeros(Float64, nCells)
num = zeros(Float64, nCells)
diff = zeros(Float64, nCells)

# Spherical harmonics (Y_n=4^m=3)
m = 3
n = 4

    # P_n=4_m=3 (z) = -105z(1-z^2)^(3/2)

for iCell in 1:nCells
    z = sin(latCell[iCell])
    a = exp(m*lonCell[iCell]im) * (-105.0*z*(1.0-z^2)^(1.5))
    init[iCell] = a.re

    # Analytic solution
    anl[iCell] = -n*(n+1)*a.re
end

###########################################################

# Compute div(grad(phi)) : Finite-volume

#for iCell in 1:nCells
#
#    grad = 0
#    for i in 1:nEdgesOnCell[iCell]
#        iEdge = edgesOnCell[i,iCell]
#        cell1 = cellsOnEdge[1,iEdge]
#        cell2 = cellsOnEdge[2,iEdge]
#
#        # grad
#        div = - (init[cell2] - init[cell1]) / dcEdge[iEdge]
# 
#        # div
#        grad = grad + edgeSignOnCell[i,iCell] * div * dvEdge[iEdge]
#    end
#
#    num[iCell] = grad / areaCell[iCell]
#end

# Compute div(grad(phi)) : Finite-difference


#println(size(cellsOnCell))

for iCell in 1:nCells
    sum1 = 0.0
    sum2 = 0.0 
    sumArea = 0.0
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        cell0 = cellsOnCell[i,iCell]

        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]

        alpha = dvEdge[iEdge] / dcEdge[iEdge]
 
#       sum1 = sum1 + alpha * init[cell1]
#       sum2 = sum2 + alpha

#       sumArea = sumArea + dvEdge[iEdge] * dcEdge[iEdge]
         
        num[iCell] = num[iCell] + ( alpha * init[cell0] 
                                -   alpha * init[iCell] ) / areaCell[iCell]


        aa = - (init[cell2] - init[cell1]) * alpha * edgeSignOnCell[i,iCell]
        bb = ( alpha * init[cell0] - alpha * init[iCell] )

        println(aa-bb)

    end
#   num[iCell] = (sum1 - sum2*init[iCell]) / (sumArea/4.0)

    #num[iCell] = (sum1 - sum2*init[iCell]) / areaCell[iCell]
end



###########################################################

###########################################################

# Output
diff[:] = num[:]-anl[:]

out_file_name = "./output.nc"
ds_out = Dataset(out_file_name,"c")

defDim(ds_out,"nCells",nCells)

nc_var = defVar(ds_out,"lonCell",Float64,("nCells",))
nc_var[:] = lonCell[:]
nc_var = defVar(ds_out,"latCell",Float64,("nCells",))
nc_var[:] = latCell[:]
nc_var = defVar(ds_out,"init",Float64,("nCells",))
nc_var[:] = init[:]
nc_var = defVar(ds_out,"anl",Float64,("nCells",))
nc_var[:] = anl[:]
nc_var = defVar(ds_out,"num",Float64,("nCells",))
nc_var[:] = num[:]
nc_var = defVar(ds_out,"diff",Float64,("nCells",))
nc_var[:] = diff[:]

close(ds_out)

#print(latCell)

###########################################################

# Compute L2 norm error

for iCell in 1:nCells
    diff[iCell] = diff[iCell]^2
    anl[iCell] = anl[iCell]^2
end

asum = sum(diff)
bsum = sum(anl)

l2norm = sqrt(asum/bsum)

println("L2 norm error =" , l2norm)
