using NCDatasets
#using Plots
using PyPlot

###########################################################

mesh_file_name = "../../meshes/x1.2562.grid.nc"   # Miniweather output file name
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
areaEdge = zeros(Float64, nEdges)
 
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

for iEdge in 1:nEdges
    areaEdge[iEdge] = dcEdge[iEdge] * dvEdge[iEdge]
end

############
println(sum(areaEdge))
println(4.0*pi)

numCell = zeros(Float64,  nCells)
for iCell in 1:nCells
    sum1 = 0.0
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        sum1 = sum1 + areaEdge[iEdge] * areaEdge[iEdge] / 4.0
#       sum1 = sum1 + initEdge[iEdge] * areaEdge[iEdge] / 4.0
    end
    numCell[iCell] = sum1 / areaCell[iCell]
end
println(sum(numCell))
############


###########################################################

# Initial field = Spherical harmonics (Y_n^m)
# del2 Analytic solution = -n(n+1)*Y_n^m

initCell = zeros(Float64, nCells)
initEdge = zeros(Float64, nEdges)
numEdge = zeros(Float64,  nEdges)
numCell = zeros(Float64,  nCells)

# Spherical harmonics (Y_n=4^m=3)
m = 3
n = 4

    # P_n=4_m=3 (z) = -105z(1-z^2)^(3/2)

for iEdge in 1:nEdges
    z = sin(latEdge[iEdge])
    a = exp(m*lonEdge[iEdge]im) * (-105.0*z*(1.0-z^2)^(1.5))
    initEdge[iEdge] = a.re
end

for iCell in 1:nCells
    z = sin(latCell[iCell])
    a = exp(m*lonCell[iCell]im) * (-105.0*z*(1.0-z^2)^(1.5))
    initCell[iCell] = a.re
end

###########################################################

# Compute grad(div(phi))

#for iCell in 1:nCells
#
#    invAreaCell = 1.0 / areaCell[iCell]
#
#    grad = 0
#    for i in 1:nEdgesOnCell[iCell]
#        iEdge = edgesOnCell[i,iCell]
#        cell1 = cellsOnEdge[1,iEdge]
#        cell2 = cellsOnEdge[2,iEdge]
#
#        # grad
#        div = (init[cell2] - init[cell1]) / dcEdge[iEdge]
# 
#        # div
#        grad = grad - edgeSignOnCell[i,iCell] * div * dvEdge[iEdge] * invAreaCell
#    end
#
#    num[iCell] = grad
#end


###########################################################

# Interpolation - Cell => Edge

for iEdge in 1:nEdges
    cell1 = cellsOnEdge[1,iEdge]
    cell2 = cellsOnEdge[2,iEdge]
    numEdge[iEdge] = (initCell[cell1]+initCell[cell2])/2.0
end

# Compute L2 norm error
global asum = 0.0
global bsum = 0.0
for iEdge in 1:nEdges
    global asum += (numEdge[iEdge] - initEdge[iEdge])^2.0
    global bsum +=                   initEdge[iEdge]^2.0
end
l2norm = sqrt(asum/bsum)
println("L2 norm error (Cell => Edge) =" , l2norm)

###########################################################

# Interpolation - Edge => Cell

for iCell in 1:nCells
    sum1 = 0.0
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
#       sum1 = sum1 + numEdge[iEdge] * areaEdge[iEdge] / 4.0
        sum1 = sum1 + initEdge[iEdge] * areaEdge[iEdge] / 4.0
    end
    numCell[iCell] = sum1 / areaCell[iCell]
end

# Compute L2 norm error
global asum = 0.0
global bsum = 0.0
for iCell in 1:nCells
    global asum += (numCell[iCell] - initCell[iCell])^2.0
    global bsum +=                   initCell[iCell]^2.0
end
l2norm = sqrt(asum/bsum)
println("L2 norm error (Edge => Cell) =" , l2norm)
###########################################################

###########################################################


###########################################################

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

close(ds_out)


#function l2norm(exp::Array{Float,1},sol::Array{Float,2})
#end
