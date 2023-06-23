using NCDatasets
#using Plots
using PyPlot
import OffsetArrays.OffsetArray
using SHTOOLS

nxgrid = "2562"

    if ( nxgrid == "2562" )
         # dx = 480 km ~ 4deg
         const WN = 30
elseif ( nxgrid == "10242" )
         # dx = 240 km ~ 2deg
         const WN = 60
elseif ( nxgrid == "40962" )
         # dx = 120 km ~ 1deg
         const WN = 120
elseif ( nxgrid == "655362" )
         # dx = 60 km ~ 1deg
         const WN = 240
elseif ( nxgrid == "2621442" )
         # dx = 30 km ~ 1deg
         const WN = 480
elseif ( nxgrid == "10485762" )
         # dx = 15 km ~ 1deg
         const WN = 960
end

###########################################################

#mesh_file_name = "../meshes/x1."*nxgrid*".grid.nc"   # Miniweather output file name
mesh_file_name = "./x1."*nxgrid*".grid.nc"   # Miniweather output file name
ds = Dataset(mesh_file_name,"r")

###########################################################

latCell           = ds["latEdge"][:]
lonCell           = ds["lonEdge"][:]
#latCell           = ds["latCell"][:]
#lonCell           = ds["lonCell"][:]
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
dv1Edge           = ds["dv1Edge"][:]
dv2Edge           = ds["dv2Edge"][:]
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

#nCells = ds.dim["nCells"]
nCells = ds.dim["nEdges"]
nEdges = ds.dim["nEdges"]
nVertices = ds.dim["nVertices"]
maxEdges = ds.dim["maxEdges"]
maxEdges2 = ds.dim["maxEdges2"]
vertexDegree = ds.dim["vertexDegree"]

edgeSignOnCell = zeros(Float64, maxEdges,nCells)
areaEdge = zeros(Float64, nEdges)

###########################################################
 
#for iCell in 1:nCells
#    for i in 1:nEdgesOnCell[iCell]
#        iEdge = edgesOnCell[i,iCell]
#        
#        if ( iCell == cellsOnEdge[1,iEdge])
#           edgeSignOnCell[i,iCell] = -1.0
#        else
#           edgeSignOnCell[i,iCell] =  1.0
#        end
#
#    end
#    #areaCell[iCell] = areaCell[iCell] * 6371000^2
#end 

for iEdge in 1:nEdges
    areaEdge[iEdge] = dcEdge[iEdge] * dvEdge[iEdge]/2.0
#    #areaEdge[iEdge] = dcEdge[iEdge] * (dvEdge[iEdge]+dv2Edge[iEdge]+dv1Edge[iEdge])/4.0
end

###########################################################
###########################################################

# Allocation array
_alp = zeros(Float64,2*WN+1,WN+1,nCells)
alp  = OffsetArray(_alp,-WN:WN,0:WN,1:nCells)
_dalp = zeros(Float64,2*WN+1,WN+1,nCells)
dalp  = OffsetArray(_dalp,-WN:WN,0:WN,1:nCells)
_norm = zeros(Float64,WN+1,WN+1)
norm = OffsetArray(_norm,0:WN,0:WN)

# Factorials to normalization
for n = 0:WN
    for m = 0:n
        nmm = n-m
        npm = n+m
        fact_nmm = factorial(big(nmm))
        fact_npm = factorial(big(npm))
        norm_scale = sqrt( (2.0*n+1.0)*fact_nmm/(4.0*pi*fact_npm) )
        norm[m,n] = norm_scale
    end
end


println("norm 0,1 ",norm[0,1])


# First few Associated Legendre Polynomials (ALP)
for iCell in 1:nCells
    z = sin(latCell[iCell])
    x = 1.0-z^2

    alp[ 0,0,iCell] = 1.0

    alp[ 0,1,iCell] = z
    alp[ 1,1,iCell] = -sqrt(x)
    alp[-1,1,iCell] = -0.5*alp[1,1,iCell]

    alp[ 0,2,iCell] = 0.5*(3.0*z^2-1.0)
    alp[ 1,2,iCell] = -3.0*z*sqrt(x)
    alp[ 2,2,iCell] =  3.0*(x)
    alp[-2,2,iCell] = (1.0/24.0)*alp[2,2,iCell]
    alp[-1,2,iCell] =-(1.0/6.0)*alp[1,2,iCell]

    alp[ 0,3,iCell] = 0.5*(5.0*z^3-3.0*z)
    alp[ 1,3,iCell] = 1.5*(1.0-5.0*z^2)*sqrt(x)
    alp[ 2,3,iCell] = 15.0*z*(x)
    alp[ 3,3,iCell] =-15.0*x^(3.0/2.0)
    alp[-3,3,iCell] =-(1.0/720.0)*alp[3,3,iCell]
    alp[-2,3,iCell] = (1.0/120.0)*alp[2,3,iCell]
    alp[-1,3,iCell] =-(1.0/12.0)*alp[1,3,iCell]
end


# Recurrence
for iCell in 1:nCells

    z = sin(latCell[iCell])
    a = sqrt(1.0-z^2)

    # m = 0,1,2
    for m = 0:2
        for n = 4:WN
            alp[m,n,iCell] = ( z*(2.0*n-1.0)*alp[m,n-1,iCell]
                                    -(n+m-1)*alp[m,n-2,iCell]) / (n-m)
        end
    end
 

    for n = 4:WN-1
        alp[n,n  ,iCell] =  -(2.0*(n-1)+1.0)*a*alp[n-1,n-1,iCell]
        alp[n,n+1,iCell] = z*(2.0* n   +1.0)  *alp[n,n,iCell]
    end
    alp[WN,WN,iCell] =  -(2.0*(WN-1)+1.0)*a*alp[WN-1,WN-1,iCell]

    for m = 3:WN-2
        for n = m+2:WN
            alp[m,n,iCell] = ( z*(2.0*n-1.0)*alp[m,n-1,iCell]
                                    -(n+m-1)*alp[m,n-2,iCell]) / (n-m)
        end
    end

end 

# Derivation of ALP
for iCell in 1:nCells
    z = sin(latCell[iCell])
    a = sqrt(1.0-z^2)

    m = 0
    for n in m+1:WN
        dalp[m,n,iCell] = ( n*z*alp[m,n,iCell]-(n+m)*alp[m,n-1,iCell] ) / (-(a^2))
    end

    for n in 1:WN
        for m in 1:n
            dalp[m,n,iCell] = ( -(n+m)*(n-m+1)*a*alp[m-1,n,iCell]
                                -m*z*alp[m,n,iCell] ) / (-(a^2))
        end
    end
end

# Save the last mode
init = zeros(Float64, nCells)
anl = zeros(Float64, nCells)
diff = zeros(Float64, nCells)
num = zeros(Float64, nCells)
u = zeros(Float64, nCells)
v = zeros(Float64, nCells)
normalVelocity = zeros(Float64, nCells)
tangentialVelocity = zeros(Float64, nCells)
 = zeros(Float64, nCells)


#m = Integer(WN/2)
m = 5
n = 5
for iCell in 1:nCells
    a = exp(im*m*lonCell[iCell]) * alp[m,n,iCell] * norm[m,n]
    b = im*m*exp(im*m*lonCell[iCell]) * alp[m,n,iCell] * norm[m,n] / cos(latCell[iCell])
    c = exp(im*m*lonCell[iCell]) * dalp[m,n,iCell] * norm[m,n]

    init[iCell] = a.re

    u[iCell] = b.re
    v[iCell] = c.re

    # To check probability
    #anl[iCell] = (abs(a)^2.0) * areaCell[iCell]
    anl[iCell] = (abs(a)^2.0) * areaEdge[iCell] 

#   diff[iCell] = u[iCell]
end


# u,v -> normalVelocity, tangentialVelocity
for iCell in 1:nCells
        normalVelocity[iCell] =  u[iCell]*cos(angleEdge[iCell]) + v[iCell]*sin(angleEdge[iCell])
    tangentialVelocity[iCell] = -u[iCell]*sin(angleEdge[iCell]) + v[iCell]*cos(angleEdge[iCell])
end

# normalVelocity, tangentialVelocity -> u,v
for iCell in 1:nCells
    u[iCell] = normalVelocity[iCell] * cos(angleEdge[iCell]) - tangentialVelocity[iCell]*sin(angleEdge[iCell])
    v[iCell] = normalVelocity[iCell] * sin(angleEdge[iCell]) + tangentialVelocity[iCell]*cos(angleEdge[iCell])
   
#   diff[iCell] = u[iCell] - diff[iCell]
end

# Divergence
for iCell in 1:nCells

    div = 0
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]

        # grad
        div = - (init[cell2] - init[cell1]) / dcEdge[iEdge]
 
        # div
        div = div + edgeSignOnCell[i,iCell] * normalVelocity[iCell] * dvEdge[iEdge]
    end

    num[iCell] = grad / areaCell[iCell]
end

#          normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
#          tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))

#            zonalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*cos(angleEdge(iEdge)) &
#              - tangentialBarotropicVel(iEdge)*sin(angleEdge(iEdge))
#            meridionalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*sin(angleEdge(iEdge)) &
#              + tangentialBarotropicVel(iEdge)*cos(angleEdge(iEdge))


println(sum(anl)) # This must be unity with some numerical errors

############################################################
############################################################

# Output
out_file_name = "./output.nc"
ds_out = Dataset(out_file_name,"c")

defDim(ds_out,"nCells",nCells)

nc_var = defVar(ds_out,"lonCell",Float64,("nCells",))
nc_var[:] = lonCell[:]
nc_var = defVar(ds_out,"latCell",Float64,("nCells",))
nc_var[:] = latCell[:]
nc_var = defVar(ds_out,"init",Float64,("nCells",))
nc_var[:] = init[:]
nc_var = defVar(ds_out,"u",Float64,("nCells",))
nc_var[:] = u[:]
nc_var = defVar(ds_out,"v",Float64,("nCells",))
nc_var[:] = v[:]
nc_var = defVar(ds_out,"normalVelocity",Float64,("nCells",))
nc_var[:] = normalVelocity[:]

close(ds_out)

#print(latCell)

###########################################################
