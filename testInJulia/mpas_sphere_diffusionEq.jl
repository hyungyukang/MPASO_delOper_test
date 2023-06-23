using NCDatasets
#using Plots
using PyPlot

###########################################################

mesh_file_name = "x1.40962.grid.nc"   # Miniweather output file name
init_file_name = "ALP_3modes_QU120.0km_Cell.nc"

ds = Dataset(mesh_file_name,"r")
ds_init = Dataset(init_file_name,"r")


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

init = ds_init["alp"][1,:]

###########################################################

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
 
    areaCell[iCell] = areaCell[iCell] * 6371000^2
end 


###########################################################

# Initial field = Spherical harmonics (Y_n^m)
# del2 Analytic solution = -n(n+1)*Y_n^m

#init = zeros(Float64, nCells)
anl = zeros(Float64, nCells)
num = zeros(Float64, nCells)
diff = zeros(Float64, nCells)
psi_n0 = zeros(Float64, nCells)
psi_n1 = zeros(Float64, nCells)
psi_d2 = zeros(Float64, nCells)

###########################################################

# Time stepping
ntime = 16
dt = 7200.0
nu = 180000.0

# Initial condition
psi_n0[:] = init[:]
println(psi_n0[1:10])

for n = 1:ntime
 
    println("Ntime",n)
    # Compute div(grad(phi))
    for iCell in 1:nCells
        grad = 0
        for i in 1:nEdgesOnCell[iCell]
            iEdge = edgesOnCell[i,iCell]
            cell1 = cellsOnEdge[1,iEdge]
            cell2 = cellsOnEdge[2,iEdge]
            # grad
            div = - (psi_n0[cell2] - psi_n0[cell1]) / dcEdge[iEdge]
            # div
            grad = grad + edgeSignOnCell[i,iCell] * div * dvEdge[iEdge]
        end
        psi_d2[iCell] = grad / areaCell[iCell]
    end

    psi_n1[:] = psi_n0[:] + dt*nu*psi_d2[:]

    psi_n1[:] = 0.5*(psi_n0[:]+psi_n1[:])


    # Compute div(grad(phi))
    for iCell in 1:nCells
        grad = 0
        for i in 1:nEdgesOnCell[iCell]
            iEdge = edgesOnCell[i,iCell]
            cell1 = cellsOnEdge[1,iEdge]
            cell2 = cellsOnEdge[2,iEdge]
            # grad
            div = - (psi_n1[cell2] - psi_n1[cell1]) / dcEdge[iEdge]
            # div
            grad = grad + edgeSignOnCell[i,iCell] * div * dvEdge[iEdge]
        end
        psi_d2[iCell] = grad / areaCell[iCell]
    end

    psi_n1[:] = psi_n0[:] + dt*nu*psi_d2[:]

    # Stepping
    psi_n0[:] = psi_n1[:]
end

###########################################################

println(psi_n0[1:10])

###########################################################

# Output
diff[:] = num[:]-anl[:]

out_file_name = "./output_time.nc"
ds_out = Dataset(out_file_name,"c")

defDim(ds_out,"nCells",nCells)

nc_var = defVar(ds_out,"lonCell",Float64,("nCells",))
nc_var[:] = lonCell[:]
nc_var = defVar(ds_out,"latCell",Float64,("nCells",))
nc_var[:] = latCell[:]
nc_var = defVar(ds_out,"psi_n0",Float64,("nCells",))
nc_var[:] = psi_n0[:]

close(ds_out)
#print(latCell)

###########################################################
#
## Compute L2 norm error
#
#for iCell in 1:nCells
#    diff[iCell] = diff[iCell]^2
#    anl[iCell] = anl[iCell]^2
#end
#
#asum = sum(diff)
#bsum = sum(anl)
#
#l2norm = sqrt(asum/bsum)
#
#println("L2 norm error =" , l2norm)
