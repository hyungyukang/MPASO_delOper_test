using NCDatasets
#using Plots
using PyPlot
import OffsetArrays.OffsetArray
using SHTOOLS

#nxgrid = "2562"
#
#    if ( nxgrid == "2562" )
#         # dx = 480 km ~ 4deg
#         const WN = 30
#elseif ( nxgrid == "10242" )
#         # dx = 240 km ~ 2deg
#         const WN = 60
#elseif ( nxgrid == "40962" )
#         # dx = 120 km ~ 1deg
#         const WN = 120
#elseif ( nxgrid == "655362" )
#         # dx = 60 km ~ 1deg
#         const WN = 240
#elseif ( nxgrid == "2621442" )
#         # dx = 30 km ~ 1deg
#         const WN = 480
#elseif ( nxgrid == "10485762" )
#         # dx = 15 km ~ 1deg
#         const WN = 960
#end

###########################################################

#mesh_file_name = "../meshes/x1."*nxgrid*".grid.nc"   # Miniweather output file name
#mesh_file_name = "./x1."*nxgrid*".grid.nc"   # Miniweather output file name

mesh_file_name = "../../meshes/x1.2562.grid.nc"   # Miniweather output file name
const WN = 10

#mesh_file_name = "../../meshes/oQU240.nc"   # Miniweather output file name
#const WN = 60

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
edgeSignOnVertex = zeros(Float64,maxEdges, nVertices)
kiteIndexOnCell = zeros(Integer,maxEdges,nCells)

###########################################################
for iCell in 1:nCells
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        iVertex = verticesOnCell[i, iCell]
        
        if ( iCell == cellsOnEdge[1,iEdge])
           edgeSignOnCell[i,iCell] = -1.0
        else
           edgeSignOnCell[i,iCell] =  1.0
        end

        for j in 1:vertexDegree
           if (cellsOnVertex[j,iVertex] == iCell)
              kiteIndexOnCell[i,iCell] = j
           end
        end

    end
end 

#for iVertex in 1:nVertices
#    for i in 1:vertexDegree
#      iEdge = edgesOnVertex[i, iVertex]
#      if ( iEdge == 0 )
#         edgesOnVertex[i,iVertex] = nEdges+1
#      end
#    end
#end

for iVertex in 1:nVertices
    for i in 1:vertexDegree
      iEdge = edgesOnVertex[i, iVertex]

      if (iVertex == verticesOnEdge[1,iEdge])
         edgeSignOnVertex[i,iVertex] = -1.0
      else
         edgeSignOnVertex[i,iVertex] =  1.0
      end
   end
end
 
for iEdge in 1:nEdges
    areaEdge[iEdge] = dcEdge[iEdge] * dvEdge[iEdge]
end

###########################################################
###########################################################

# Allocation array
_alp = zeros(Float64,WN+1,WN+1,nEdges)
alp  = OffsetArray(_alp,0:WN,0:WN,1:nEdges)
_dalp = zeros(Float64,WN+1,WN+1,nEdges)
dalp  = OffsetArray(_dalp,0:WN,0:WN,1:nEdges)
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


# First few Associated Legendre Polynomials (ALP)
for iEdge in 1:nEdges
    z = sin(latEdge[iEdge])
    x = 1.0-z^2

    alp[ 0,0,iEdge] = 1.0

    alp[ 0,1,iEdge] = z
    alp[ 1,1,iEdge] = -sqrt(x)
#   alp[-1,1,iEdge] = -0.5*alp[1,1,iEdge]

    alp[ 0,2,iEdge] = 0.5*(3.0*z^2-1.0)
    alp[ 1,2,iEdge] = -3.0*z*sqrt(x)
    alp[ 2,2,iEdge] =  3.0*(x)
#   alp[-2,2,iEdge] = (1.0/24.0)*alp[2,2,iEdge]
#   alp[-1,2,iEdge] =-(1.0/6.0)*alp[1,2,iEdge]

    alp[ 0,3,iEdge] = 0.5*(5.0*z^3-3.0*z)
    alp[ 1,3,iEdge] = 1.5*(1.0-5.0*z^2)*sqrt(x)
    alp[ 2,3,iEdge] = 15.0*z*(x)
    alp[ 3,3,iEdge] =-15.0*x^(3.0/2.0)
#   alp[-3,3,iEdge] =-(1.0/720.0)*alp[3,3,iEdge]
#   alp[-2,3,iEdge] = (1.0/120.0)*alp[2,3,iEdge]
#   alp[-1,3,iEdge] =-(1.0/12.0)*alp[1,3,iEdge]
end


# Recurrence
for iEdge in 1:nEdges

    z = sin(latEdge[iEdge])
    a = sqrt(1.0-z^2)

    # m = 0,1,2
    for m in 0:2
        for n in 4:WN
            alp[m,n,iEdge] = ( z*(2.0*n-1.0)*alp[m,n-1,iEdge]
                                    -(n+m-1)*alp[m,n-2,iEdge]) / (n-m)
        end
    end

    for n in 4:WN-1
        alp[n,n  ,iEdge] =  -(2.0*(n-1)+1.0)*a*alp[n-1,n-1,iEdge]
        alp[n,n+1,iEdge] = z*(2.0* n   +1.0)  *alp[n,n,iEdge]
    end
    alp[WN,WN,iEdge] =  -(2.0*(WN-1)+1.0)*a*alp[WN-1,WN-1,iEdge]

    for m in 3:WN-2
        for n in m+2:WN
            alp[m,n,iEdge] = ( z*(2.0*n-1.0)*alp[m,n-1,iEdge]
                                    -(n+m-1)*alp[m,n-2,iEdge]) / (n-m)
        end
    end

    for n in 0:WN
        for m in 0:n
            alp[m,n,iEdge] = norm[m,n] * alp[m,n,iEdge] * sqrt(2.0)
        end
    end

end 

# Derivation of ALP
for iEdge in 1:nEdges
    z = sin(latEdge[iEdge])
    a = sqrt(1.0-z^2)

    m = 0
    for n in m+1:WN
        dalp[m,n,iEdge] = ( n*z*alp[m,n,iEdge]-(n+m)*alp[m,n-1,iEdge] ) / (-(a^2))
    end

    for n in 1:WN
        for m in 1:n
            dalp[m,n,iEdge] = ( -(n+m)*(n-m+1)*a*alp[m-1,n,iEdge]
                                -m*z*alp[m,n,iEdge] ) / (-(a^2))
        end
    end
end

# Save the last mode
init = zeros(Float64, nEdges)
anl = zeros(Float64, nEdges)
diff = zeros(Float64, nEdges)
num = zeros(Float64, nEdges)
u = zeros(Float64, nEdges)
v = zeros(Float64, nEdges)
u_anl = zeros(Float64, nEdges)
v_anl = zeros(Float64, nEdges)
normalVelocity = zeros(Float64, nEdges)
tangentialVelocity = zeros(Float64, nEdges)
divergence = zeros(Float64, nCells)
divergenceVert = zeros(Float64, nVertices)
divergenceEdge = zeros(Float64, nEdges)
relativeVorticity = zeros(Float64,nVertices)
relaiveVorticityCell = zeros(Float64,nCells)
relativeVorticityEdge = zeros(Float64,nEdges)
del2u = zeros(Float64, nEdges)
Vx = zeros(Float64,nEdges)
Vy = zeros(Float64,nEdges)
Vz = zeros(Float64,nEdges)

###############################################################################################
#m = Integer(WN/2)
m = 0
n = 2
for iEdge in 1:nEdges
    a =      exp(im*m*lonEdge[iEdge]) *  alp[m,n,iEdge] 
    b = im*m*exp(im*m*lonEdge[iEdge]) *  alp[m,n,iEdge]
    c =      exp(im*m*lonEdge[iEdge]) * ( dalp[m,n,iEdge] * cos(latEdge[iEdge])
                                         + alp[m,n,iEdge] * sin(latEdge[iEdge]))
    #c =      exp(im*m*lonEdge[iEdge]) * ( dalp[m,n,iEdge])

    init[iEdge] = a.re

    u_anl[iEdge] = c.re
    v_anl[iEdge] = b.re

    # To check probability
    anl[iEdge] = (abs(a)^2.0) * areaEdge[iEdge] / 4.0

#   diff[iEdge] = u[iEdge]
end
println("int|Y| = ",sum(anl))
###############################################################################################

# u,v -> normalVelocity, tangentialVelocity
for iEdge in 1:nEdges
        normalVelocity[iEdge] =  u_anl[iEdge]*cos(angleEdge[iEdge]) + v_anl[iEdge]*sin(angleEdge[iEdge])
    tangentialVelocity[iEdge] = -u_anl[iEdge]*sin(angleEdge[iEdge]) + v_anl[iEdge]*cos(angleEdge[iEdge])
end

#####################################
# Divergence - way 1: Default
#for iCell in 1:nCells
#    div = 0.0
#    for i in 1:nEdgesOnCell[iCell]
#        iEdge = edgesOnCell[i,iCell]
#        # div
#        r_tmp = dvEdge[iEdge]*normalVelocity[iEdge]
#        div = div + edgeSignOnCell[i,iCell] * r_tmp
#    end
#    divergence[iCell] = div / areaCell[iCell]
#end
#####################################

## Divergence - way 2: Divergence at vertices
#for iVertex in 1:nVertices
#    invAreaTri1 = 1.0 / areaTriangle[iVertex]
#    divergenceVert[iVertex] = 0.0
#    for i in 1:vertexDegree
#        iEdge = edgesOnVertex[i,iVertex]
#        r_tmp = dcEdge[iEdge] * tangentialVelocity[iEdge]
#        divergenceVert[iVertex] = divergenceVert[iVertex] - edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
#
##       r_tmp = dcEdge[iEdge] * normalVelocity[iEdge]
##       divergenceVert[iVertex] = divergenceVert[iVertex] + edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
#    end
#end 
#
## Interpolation - Vertex => Cell
#for iCell in 1:nCells
#    invAreaCell1 = 1.0 / areaCell[iCell]
#    divergence[iCell] = 0.0
#    for i in 1:nEdgesOnCell[iCell]
#        j = kiteIndexOnCell[i,iCell]
#        iVertex = verticesOnCell[i,iCell]
#        divergence[iCell] =divergence[iCell] +
#                           kiteAreasOnVertex[j,iVertex] * divergenceVert[iVertex]*invAreaCell1
#    end
#end
#
## Interpolation - Vertex => Edge
#for iEdge in 1:nEdges
#    iVertex1 = verticesOnEdge[1,iEdge]
#    iVertex2 = verticesOnEdge[2,iEdge]
#    divergenceEdge[iEdge] = (divergenceVert[iVertex1]+divergenceVert[iVertex2])*0.5
#    #divergenceEdge[iEdge] = (divergenceVert[iVertex1]+divergenceVert[iVertex2])*dcEdge[iEdge]
#    #divergenceEdge[iEdge] = dcEdge[iEdge]*dvEdge[iEdge]
#    #divergenceEdge[iEdge] = dcEdge[iEdge]
#   
#end

#println(divergence[1:50])


## Relative Vorticity
#for iVertex in 1:nVertices
#    invAreaTri1 = 1.0 / areaTriangle[iVertex]
#    relativeVorticity[iVertex] = 0.0
#    for i in 1:vertexDegree
#        iEdge = edgesOnVertex[i,iVertex]
#        r_tmp = dcEdge[iEdge] * normalVelocity[iEdge]
#        relativeVorticity[iVertex] = relativeVorticity[iVertex] + edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
#    end
#end 

###################################
## Del2
#for iEdge in 1:nEdges
#   cell1 = cellsOnEdge[1,iEdge]
#   cell2 = cellsOnEdge[2,iEdge]
#   vertex1 = verticesOnEdge[1,iEdge]
#   vertex2 = verticesOnEdge[2,iEdge]
#
#   dcEdgeInv = 1.0 / dcEdge[iEdge]
#   dvEdgeInv = 1.0 / dvEdge[iEdge]
#
#   #uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv - (relativeVorticity[vertex2]-relativeVorticity[vertex1])*dvEdgeInv
#   #uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv 
#
#   #println((divergence[cell2] - divergence[cell1])*dcEdgeInv,',',(relativeVorticity[vertex2]-relativeVorticity[vertex1])*dvEdgeInv)
#   uDiff =-(relativeVorticity[vertex2]-relativeVorticity[vertex1])*dvEdgeInv
#   #uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv
#
#   del2u[iEdge] = uDiff
#   #println((divergence[cell2] - divergence[cell1])*dcEdgeInv,',',(relativeVorticity[vertex2]-relativeVorticity[vertex1])*dvEdgeInv)
#   #println(del2u[iEdge])
#end
#
## Tangential velocity from normalVelocity
#for iEdge in 1:nEdges
#    tangentialVelocity[iEdge] = 0.0
#    for i in 1:nEdgesOnEdge[iEdge]
#        eoe = edgesOnEdge[i,iEdge]
#        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * del2u[eoe]
#    end
#end
#
#
## normalVelocity, tangentialVelocity -> u,v
for iEdge in 1:nEdges
#    u[iEdge] = del2u[iEdge] * cos(angleEdge[iEdge]) - tangentialVelocity[iEdge]*sin(angleEdge[iEdge])
#    v[iEdge] = del2u[iEdge] * sin(angleEdge[iEdge]) + tangentialVelocity[iEdge]*cos(angleEdge[iEdge])
     u[iEdge] = 0.0
     v[iEdge] = 0.0
end

################################################3
## Vcart
#for iEdge in 1:nEdges
##   sinlon=sin(lonEdge[iEdge])
##   coslon=cos(lonEdge[iEdge])
##   sinlat=sin(latEdge[iEdge])
##   coslat=cos(latEdge[iEdge])
#
##   Vx[iEdge] = -n*(n+1)*(-sinlon * u_anl[iEdge] - sinlat*coslon*v_anl[iEdge])
##   Vy[iEdge] = -n*(n+1)*( coslon * u_anl[iEdge] - sinlat*sinlon*v_anl[iEdge])
##   Vz[iEdge] = -n*(n+1)*(                                coslat*v_anl[iEdge])
#
##   u_anl[iEdge] = -sinlon*Vx[iEdge] + coslon*Vy[iEdge]
##   v_anl[iEdge] = -sinlat*coslon*Vx[iEdge] - sinlat*sinlon*Vy[iEdge] + coslat*Vz[iEdge]
#
#    u_anl[iEdge] = -n*(n+1)*u_anl[iEdge]
#    v_anl[iEdge] = -n*(n+1)*v_anl[iEdge]
#end
################################################3

## Vcart -> UV
#
##
#for iEdge in 1:20
#    println(iEdge,',',u[iEdge],',',u_anl[iEdge],',', u_anl[iEdge]/u[iEdge])
#    #println(iEdge,',',v[iEdge])
#end


## Tangential velocity from normalVelocity
#for iEdge in 1:nEdges
#    for i in 1:nEdgesOnEdge[iEdge]
#        eoe = edgesOnEdge[i,iEdge]
#        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * init[eoe]
#    end
#end
#
## Vector reconstruction
## normalVelocity, tangentialVelocity -> u,v
#for iEdge in 1:nEdges
#    u[iEdge] = init[iEdge] * cos(angleEdge[iEdge]) - tangentialVelocity[iEdge]*sin(angleEdge[iEdge])
#    v[iEdge] = init[iEdge] * sin(angleEdge[iEdge]) + tangentialVelocity[iEdge]*cos(angleEdge[iEdge])
#end
#
## UV -> Vxyz
#for iEdge in 1:nEdges
#    sinlon = sin(lonEdge[iEdge])
#    sinlat = sin(latEdge[iEdge])
#    coslon = cos(lonEdge[iEdge])
#    coslat = cos(latEdge[iEdge])
#    
#    
#    Vx[iEdge] = -sinlon * u[iEdge] -sinlat*coslon*v[iEdge]
#    Vy[iEdge] =  coslon * u[iEdge] -sinlat*sinlon*v[iEdge]
#    Vz[iEdge] =                     coslat       *v[iEdge]
#
#    Vx[iEdge] = -n*(n+1) * Vx[iEdge]
#    Vy[iEdge] = -n*(n+1) * Vy[iEdge]
#    Vz[iEdge] = -n*(n+1) * Vz[iEdge]
#end
#
#
## Vxyz -> UV
#for iEdge in 1:nEdges
#    sinlon = sin(lonEdge[iEdge])
#    sinlat = sin(latEdge[iEdge])
#    coslon = cos(lonEdge[iEdge])
#    coslat = cos(latEdge[iEdge])
#
#    u[iEdge] = -sinlon       *Vx[iEdge] +coslon       *Vy[iEdge]
#    v[iEdge] = -sinlat*coslon*Vx[iEdge] -sinlat*sinlon*Vy[iEdge] +coslat*Vz[iEdge]
#end
#
## UV -> norm, tan
#for iEdge in 1:nEdges
#         normalVelocity[iEdge] =  u[iEdge]*cos(angleEdge[iEdge]) + v[iEdge]*sin(angleEdge[iEdge])
#     tangentialVelocity[iEdge] = -u[iEdge]*sin(angleEdge[iEdge]) + v[iEdge]*cos(angleEdge[iEdge])
#end



# Laplacian
## u,v -> normalVelocity, tangentialVelocity
#for iEdge in 1:nEdges
#        normalVelocity[iEdge] =  u[iEdge]*cos(angleEdge[iEdge]) + v[iEdge]*sin(angleEdge[iEdge])
#    tangentialVelocity[iEdge] = -u[iEdge]*sin(angleEdge[iEdge]) + v[iEdge]*cos(angleEdge[iEdge])
#end
#
## normalVelocity, tangentialVelocity -> u,v
#for iEdge in 1:nEdges
#    u[iEdge] = normalVelocity[iEdge] * cos(angleEdge[iEdge]) - tangentialVelocity[iEdge]*sin(angleEdge[iEdge])
#    v[iEdge] = normalVelocity[iEdge] * sin(angleEdge[iEdge]) + tangentialVelocity[iEdge]*cos(angleEdge[iEdge])
#end


#for iEdge in 1:50
#    println(iEdge,',',del2u[iEdge],',',anl[iEdge],',',normalVelocity[iEdge])
#    println(iEdge,',',del2u[iEdge],',',anl[iEdge],',',init[iEdge])
#    println(iEdge,',',dcEdge[iEdge],',',dvEdge[iEdge])
#    println(iEdge,',',divergence[iEdge],',',relativeVorticity[iEdge])
#end





## Interpolation - Vertex => Cell
#for iCell in 1:nCells
#    invAreaCell1 = 1.0 / areaCell[iCell]
#    for i in 1:nEdgesOnCell[iCell]
#        j = kiteIndexOnCell[i,iCell]
#        iVertex = verticesOnCell[i,iCell]
#        relativeVorticityCell[iCell] = relativeVorticityCell[iCell]
#                                     + kiteAreasOnVertex[j,iVertex] * relativeVorticity[iVertex]*invAreaCell1
#    end
#end

###########################################################

# Interpolation - Cell => Edge
#for iEdge in 1:nEdges
#    cell1 = cellsOnEdge[1,iEdge]
#    cell2 = cellsOnEdge[2,iEdge]
#    relativeVorticityEdge[iEdge] = ( relativeVorticityCell[cell1]
#                                    +relativeVorticityCell[cell2])/2.0
#    anl[iEdge] =  +(n*(n+1)) * init[iEdge]
#    diff[iEdge] = abs(relativeVorticityEdge[iEdge] - anl[iEdge])
#end
#
#println("Rel Vol error = ",sum(diff))

#l2norm = sqrt(asum/bsum)
#println("L2 norm error (Cell => Edge) =" , l2norm)

###########################################################

#          normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
#          tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))

#            zonalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*cos(angleEdge(iEdge)) &
#              - tangentialBarotropicVel(iEdge)*sin(angleEdge(iEdge))
#            meridionalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*sin(angleEdge(iEdge)) &
#              + tangentialBarotropicVel(iEdge)*cos(angleEdge(iEdge))

#println(normalVelocity[1:100])
#println(divergence[1:100])

#println(sum(diff)) # This must be unity with some numerical errors

############################################################
############################################################

# Output
out_file_name = "./output.nc"
ds_out = Dataset(out_file_name,"c")

defDim(ds_out,"nCells",nCells)
defDim(ds_out,"nEdges",nEdges)
defDim(ds_out,"nVertices",nVertices)

nc_var = defVar(ds_out,"lonCell",Float64,("nCells",))
nc_var[:] = lonCell[:]
nc_var = defVar(ds_out,"latCell",Float64,("nCells",))
nc_var[:] = latCell[:]
nc_var = defVar(ds_out,"lonEdge",Float64,("nEdges",))
nc_var[:] = lonEdge[:]
nc_var = defVar(ds_out,"latEdge",Float64,("nEdges",))
nc_var[:] = latEdge[:]
nc_var = defVar(ds_out,"lonVertex",Float64,("nVertices",))
nc_var[:] = lonVertex[:]
nc_var = defVar(ds_out,"latVertex",Float64,("nVertices",))
nc_var[:] = latVertex[:]
nc_var = defVar(ds_out,"init",Float64,("nEdges",))
nc_var[:] = init[:]
nc_var = defVar(ds_out,"u",Float64,("nEdges",))
nc_var[:] = u[:]
nc_var = defVar(ds_out,"v",Float64,("nEdges",))
nc_var[:] = v[:]
nc_var = defVar(ds_out,"uAnl",Float64,("nEdges",))
nc_var[:] = u_anl[:]
nc_var = defVar(ds_out,"vAnl",Float64,("nEdges",))
nc_var[:] = v_anl[:]
nc_var = defVar(ds_out,"normalVelocity",Float64,("nEdges",))
nc_var[:] = normalVelocity[:]
nc_var = defVar(ds_out,"divergence",Float64,("nCells",))
nc_var[:] = divergence[:]
nc_var = defVar(ds_out,"divergenceVert",Float64,("nVertices",))
nc_var[:] = divergenceVert[:]
nc_var = defVar(ds_out,"divergenceEdge",Float64,("nEdges",))
nc_var[:] = divergenceEdge[:]
nc_var = defVar(ds_out,"anl",Float64,("nEdges",))
nc_var[:] = anl[:]
nc_var = defVar(ds_out,"del2u",Float64,("nEdges",))
nc_var[:] = del2u[:]

close(ds_out)

###########################################################
