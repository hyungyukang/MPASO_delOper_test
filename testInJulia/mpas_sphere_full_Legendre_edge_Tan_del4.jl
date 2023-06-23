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

#mesh_file_name = "../../meshes/oQU240.nc"   # Miniweather output file name

#mesh_file_name = "../../meshes/x1.2562.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.10242.grid.nc"   # Miniweather output file name
mesh_file_name = "../../meshes/x1.40962.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.163842.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.655362.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.2621442.grid.nc"   # Miniweather output file name
#mesh_file_name = "../../meshes/x1.10485762.grid.nc"   # Miniweather output file name

const WN = 5

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

for iVertex in 1:nVertices
    for i in 1:vertexDegree
      iEdge = edgesOnVertex[i, iVertex]

      if ( iEdge > 0 ) 

      if (iVertex == verticesOnEdge[1,iEdge])
         edgeSignOnVertex[i,iVertex] = -1.0
      else
         edgeSignOnVertex[i,iVertex] =  1.0
      end

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
    x = sin(latEdge[iEdge])
    z = 1.0-x^2

    alp[ 0,0,iEdge] = 1.0

    alp[ 0,1,iEdge] = x
    alp[ 1,1,iEdge] = -sqrt(z)

    alp[ 0,2,iEdge] = 0.5*(3.0*x^2-1.0)
    alp[ 1,2,iEdge] = -3.0*x*sqrt(z)
    alp[ 2,2,iEdge] =  3.0*(z)

    alp[ 0,3,iEdge] = 0.5*(5.0*x^3-3.0*x)
    alp[ 1,3,iEdge] = 1.5*(1.0-5.0*x^2)*sqrt(z)
    alp[ 2,3,iEdge] = 15.0*x*(z)
    alp[ 3,3,iEdge] =-15.0*z^(1.5)

    alp[ 0,4,iEdge] = 0.125*(35.0*x^4-30.0*x^2+3.0)
    alp[ 1,4,iEdge] = -2.5*(7.0*x^3-3.0*x)*sqrt(z)
    alp[ 2,4,iEdge] =  7.5*(7.0*x^2-1.0)*z
    alp[ 3,4,iEdge] = -105.0*x*z^1.5
    alp[ 4,4,iEdge] =  105.0*z^2
end


# Recurrence
#for iEdge in 1:nEdges
#
#    z = sin(latEdge[iEdge])
#    a = sqrt(1.0-z^2)
#
#    # m = 0,1,2
#    for m in 0:2
#        for n in 4:WN
#            alp[m,n,iEdge] = ( z*(2.0*n-1.0)*alp[m,n-1,iEdge]
#                                    -(n+m-1)*alp[m,n-2,iEdge]) / (n-m)
#        end
#    end
#
#    for n in 4:WN-1
#        alp[n,n  ,iEdge] =  -(2.0*(n-1)+1.0)*a*alp[n-1,n-1,iEdge]
#        alp[n,n+1,iEdge] = z*(2.0* n   +1.0)  *alp[n,n,iEdge]
#    end
#    alp[WN,WN,iEdge] =  -(2.0*(WN-1)+1.0)*a*alp[WN-1,WN-1,iEdge]
#
#    for m in 3:WN-2
#        for n in m+2:WN
#            alp[m,n,iEdge] = ( z*(2.0*n-1.0)*alp[m,n-1,iEdge]
#                                    -(n+m-1)*alp[m,n-2,iEdge]) / (n-m)
#        end
#    end
#
#    for n in 0:WN
#        for m in 0:n
#            alp[m,n,iEdge] = norm[m,n] * alp[m,n,iEdge] * sqrt(2.0)
#        end
#    end
#
#end 

# Derivation of ALP
for iEdge in 1:nEdges
    x = sin(latEdge[iEdge])
    z = 1.0-x^2

    m = 0
    for n in m+1:WN
        dalp[m,n,iEdge] = ( n*x*alp[m,n,iEdge]-(n+m)*alp[m,n-1,iEdge] ) / (-z)
    end

    for n in 1:WN
        for m in 1:n
            dalp[m,n,iEdge] = ( -(n+m)*(n-m+1)*sqrt(z)*alp[m-1,n,iEdge]
                                -m*x*alp[m,n,iEdge] ) / (-z)
             
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
uLap_anl = zeros(Float64, nEdges)
vLap_anl = zeros(Float64, nEdges)
vx_anl = zeros(Float64, nEdges)
vy_anl = zeros(Float64, nEdges)
vz_anl = zeros(Float64, nEdges)
normalVelocity = zeros(Float64, nEdges)
normalVelocityAnl = zeros(Float64, nEdges)
tangentialVelocityAnl = zeros(Float64, nEdges)
tangentialVelocity = zeros(Float64, nEdges)
tangentialVelocityTemp = zeros(Float64, nEdges)
divergence = zeros(Float64, nCells)
divergenceVert = zeros(Float64, nVertices)
divergenceEdge = zeros(Float64, nEdges)
relativeVorticity = zeros(Float64,nVertices)
relativeVorticityCell = zeros(Float64,nCells)
relativeVorticityEdge = zeros(Float64,nEdges)
relativeVorticityAnl = zeros(Float64,nVertices)
relativeVorticityAnlEdge = zeros(Float64,nEdges)
del2u = zeros(Float64, nEdges)
Vx = zeros(Float64,nEdges)
Vy = zeros(Float64,nEdges)
Vz = zeros(Float64,nEdges)

###############################################################################################
#m = Integer(WN/2)
m = 2
n = 3
#lap = -n * (n+1)
lap = (-n * (n+1))^2

for iEdge in 1:nEdges

    a =  exp(im*m*lonEdge[iEdge]) *  alp[m,n,iEdge] 
    b = m*a*im / cos(latEdge[iEdge])
    c =  exp(im*m*lonEdge[iEdge]) * ( dalp[m,n,iEdge])

    init[iEdge] = a.re

    u_anl[iEdge] = -c.re * cos(latEdge[iEdge]) 
    v_anl[iEdge] =  b.re

    uLap_anl[iEdge] = -lap * c.re * cos(latEdge[iEdge])
    vLap_anl[iEdge] =  lap * b.re

    relativeVorticityAnlEdge[iEdge] = lap * a.re

    # To check probability
    anl[iEdge] = (abs(a)^2.0) * areaEdge[iEdge] / 4.0

#   diff[iEdge] = u[iEdge]
end
println("int|Y| = ",sum(anl))
###############################################################################################

# Vxyz -> UV
#for iEdge in 1:nEdges
#    sinlon = sin(lonEdge[iEdge])
#    sinlat = sin(latEdge[iEdge])
#    coslon = cos(lonEdge[iEdge])
#    coslat = cos(latEdge[iEdge])
#
#    u_anl[iEdge] = -sinlon       *vx_anl[iEdge] +coslon       *vy_anl[iEdge]
#    v_anl[iEdge] = -sinlat*coslon*vx_anl[iEdge] -sinlat*sinlon*vy_anl[iEdge] +coslat*vz_anl[iEdge]
#    
#    lap = -n * (n+1) #* 1.24321792
#
#    uLap_anl[iEdge] = -sinlon       *lap*vx_anl[iEdge] +coslon       *lap*vy_anl[iEdge]
#    #vLap_anl[iEdge] = -sinlat*coslon*lap*vx_anl[iEdge] -sinlat*sinlon*lap*vy_anl[iEdge] +coslat*lap*vz_anl[iEdge]
#    vLap_anl[iEdge] = -sinlat*coslon*lap*vx_anl[iEdge] -sinlat*sinlon*lap*vy_anl[iEdge] +coslat*lap*vz_anl[iEdge]
#end

# u,v -> normalVelocity, tangentialVelocity
for iEdge in 1:nEdges
#   sinlon = sin(lonEdge[iEdge])
#   sinlat = sin(latEdge[iEdge])
#   coslon = cos(lonEdge[iEdge])
#   coslat = cos(latEdge[iEdge])

#   clon_clat = cos(lonEdge[iEdge]) / cos(latEdge[iEdge])
#   slon_clat = sin(lonEdge[iEdge]) / cos(latEdge[iEdge])
#   slon_clat_clat = (sin(lonEdge[iEdge])/cos(latEdge[iEdge])) / cos(latEdge[iEdge])
#   clon_clat_clat = (cos(lonEdge[iEdge])/cos(latEdge[iEdge])) / cos(latEdge[iEdge])

        normalVelocity[iEdge] =  u_anl[iEdge]*cos(angleEdge[iEdge]) + v_anl[iEdge]*sin(angleEdge[iEdge])
    tangentialVelocity[iEdge] = -u_anl[iEdge]*sin(angleEdge[iEdge]) + v_anl[iEdge]*cos(angleEdge[iEdge])

        normalVelocityAnl[iEdge] =  uLap_anl[iEdge]*cos(angleEdge[iEdge]) + vLap_anl[iEdge]*sin(angleEdge[iEdge])
    tangentialVelocityAnl[iEdge] = -uLap_anl[iEdge]*sin(angleEdge[iEdge]) + vLap_anl[iEdge]*cos(angleEdge[iEdge])

end

#####################################
# Divergence - way 1: Default
for iCell in 1:nCells
    divergence[iCell] = 0.0
    invAreaCell = 1.0 / areaCell[iCell]
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        r_tmp = dvEdge[iEdge]*normalVelocity[iEdge]
        divergence[iCell] = divergence[iCell] - edgeSignOnCell[i,iCell] * r_tmp * invAreaCell
    end
end

# Relative Vorticity
for iVertex in 1:nVertices
    invAreaTri1 = 1.0 / areaTriangle[iVertex]
    relativeVorticity[iVertex] = 0.0
    for i in 1:vertexDegree
        iEdge = edgesOnVertex[i,iVertex]
        r_tmp = dcEdge[iEdge] * normalVelocity[iEdge]
        relativeVorticity[iVertex] = relativeVorticity[iVertex] + edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
    end
end 
#
#####################################
## Del2
for iEdge in 1:nEdges
   cell1 = cellsOnEdge[1,iEdge]
   cell2 = cellsOnEdge[2,iEdge]
   vertex1 = verticesOnEdge[1,iEdge]
   vertex2 = verticesOnEdge[2,iEdge]

   dcEdgeInv = 1.0 / dcEdge[iEdge]
   dvEdgeInv = 1.0 / dvEdge[iEdge]

   uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv - (relativeVorticity[vertex2] - relativeVorticity[vertex1])*dvEdgeInv

   #uDiff = - (relativeVorticity[vertex2] - relativeVorticity[vertex1])*dvEdgeInv
   del2u[iEdge] = uDiff

end

#####################################
# Tangential velocity from normalVelocity
for iEdge in 1:nEdges
    tangentialVelocity[iEdge] = 0.0
    for i in 1:nEdgesOnEdge[iEdge]
        eoe = edgesOnEdge[i,iEdge]
        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * del2u[eoe]
    end
end
###########################################################
###########################################################
###########################################################

# divVert from tangential velocity
for iVertex in 1:nVertices
    invAreaTri1 = 1.0 / areaTriangle[iVertex]
    divergenceVert[iVertex] = 0.0
    for i in 1:vertexDegree
        iEdge = edgesOnVertex[i,iVertex]

        if ( iEdge > 0 )
        r_tmp = dcEdge[iEdge] * tangentialVelocity[iEdge]
        divergenceVert[iVertex] = divergenceVert[iVertex] + edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
        end
    end
end 

# RelVorCell from tangential velocity
for iCell in 1:nCells
    relativeVorticityCell[iCell] = 0.0
    invAreaCell = 1.0 / areaCell[iCell]
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        r_tmp = dvEdge[iEdge]*tangentialVelocity[iEdge]
        relativeVorticityCell[iCell] = relativeVorticityCell[iCell] - edgeSignOnCell[i,iCell] * r_tmp * invAreaCell
    end
end

# Del2 - ver 2
for iEdge in 1:nEdges
    cell1 = cellsOnEdge[1,iEdge]
    cell2 = cellsOnEdge[2,iEdge]
    vertex1 = verticesOnEdge[1,iEdge]
    vertex2 = verticesOnEdge[2,iEdge]
 
    dcEdgeInv = 1.0 / dcEdge[iEdge]
    dvEdgeInv = 1.0 / dvEdge[iEdge]
 
    uDiff =(divergenceVert[vertex2]-divergenceVert[vertex1])*dvEdgeInv -(relativeVorticityCell[cell2] - relativeVorticityCell[cell1]) * dcEdgeInv

    del2u[iEdge] = uDiff

end

# Tangential velocity from normalVelocity
for iEdge in 1:nEdges
    tangentialVelocity[iEdge] = 0.0
    for i in 1:nEdgesOnEdge[iEdge]
        eoe = edgesOnEdge[i,iEdge]
        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * del2u[eoe]
    end
end

#------------------------
for iEdge in 1:nEdges
    del2u[iEdge] = tangentialVelocity[iEdge]
end

# Tangential velocity from normalVelocity
for iEdge in 1:nEdges
    tangentialVelocity[iEdge] = 0.0
    for i in 1:nEdgesOnEdge[iEdge]
        eoe = edgesOnEdge[i,iEdge]
        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * del2u[eoe]
    end
end


#
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



#------------------------------------
##Cell -> Vertex
#for iVertex in 1:nVertices
#    invAreaTri1 = 1.0 / areaTriangle[iVertex]
#    areaSum = 0.0 
#    divergenceVert[iVertex] = 0.0
#    for i in 1:vertexDegree
#        iCell = cellsOnVertex[i,iVertex]
#        areaSum = areaSum + kiteAreasOnVertex[i,iVertex]
#        divergenceVert[iVertex] = divergenceVert[iVertex] + kiteAreasOnVertex[i,iVertex] * divergence[iCell]
#    end
#    divergenceVert[iVertex] = divergenceVert[iVertex] / areaSum
#end 
#
##------------------------------------
##Vertex -> Cell
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

#####################################
##------------------------------------
##Vertex -> Cell
#for iCell in 1:nCells
#    invAreaCell1 = 1.0 / areaCell[iCell]
#    relativeVorticityCell[iCell] = 0.0
#    for i in 1:nEdgesOnCell[iCell]
#        j = kiteIndexOnCell[i,iCell]
#        iVertex = verticesOnCell[i,iCell]
#        relativeVorticityCell[iCell] =relativeVorticityCell[iCell] +
#                           kiteAreasOnVertex[j,iVertex] * relativeVorticity[iVertex]*invAreaCell1
#    end
#end

####------------------------------------
## Cell -> Vertex
#for iVertex in 1:nVertices
#    invAreaTri1 = 1.0 / areaTriangle[iVertex]
#    areaSum = 0.0 
#    relativeVorticity[iVertex] = 0.0
#    for i in 1:vertexDegree
#        iCell = cellsOnVertex[i,iVertex]
#        areaSum = areaSum + kiteAreasOnVertex[i,iVertex]
#        if ( iCell > 0 )
#        relativeVorticity[iVertex] = relativeVorticity[iVertex] + kiteAreasOnVertex[i,iVertex] * relativeVorticityCell[iCell]
#        end
#    end
#    relativeVorticity[iVertex] = relativeVorticity[iVertex] / areaSum
#end 

##------------------------------------
# Vertex -> Edge
for iEdge in 1:nEdges
   vertex1 = verticesOnEdge[1,iEdge]
   vertex2 = verticesOnEdge[2,iEdge]
   relativeVorticityEdge[iEdge] = 0.5*(relativeVorticity[vertex1]+relativeVorticity[vertex2])
end


#####################################
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
#   uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv - (relativeVorticity[vertex2] - relativeVorticity[vertex1])*dvEdgeInv
#
#   #uDiff = - (relativeVorticity[vertex2] - relativeVorticity[vertex1])*dvEdgeInv
#   del2u[iEdge] = uDiff
#
#end

## Del2 - ver 2
#for iEdge in 1:nEdges
#    cell1 = cellsOnEdge[1,iEdge]
#    cell2 = cellsOnEdge[2,iEdge]
#    vertex1 = verticesOnEdge[1,iEdge]
#    vertex2 = verticesOnEdge[2,iEdge]
# 
#    dcEdgeInv = 1.0 / dcEdge[iEdge]
#    dvEdgeInv = 1.0 / dvEdge[iEdge]
# 
#    #uDiff = -(relativeVorticityCell[cell2] - relativeVorticityCell[cell1]) * dcEdgeInv
#    uDiff =(divergenceVert[vertex2]-divergenceVert[vertex1])*dvEdgeInv -(relativeVorticityCell[cell2] - relativeVorticityCell[cell1]) * dcEdgeInv
#    del2u[iEdge] = uDiff
#
##   #uDiff = (divergence[cell2] - divergence[cell1])*dcEdgeInv - (relativeVorticity[vertex2] - relativeVorticity[vertex1])*dvEdgeInv
#end

#####################################

#
#for iEdge in 1:nEdges
#    tangentialVelocity[iEdge] = 0.0
#    for i in 1:nEdgesOnEdge[iEdge]
#        eoe = edgesOnEdge[i,iEdge]
#        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * del2u[eoe]
#    end
#end

#####################################
#
#

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



## normalVelocity, tangentialVelocity -> u,v
for iEdge in 1:nEdges
#   clon_clat = cos(lonEdge[iEdge]) / cos(latEdge[iEdge])
#   slon_clat = sin(lonEdge[iEdge]) / cos(latEdge[iEdge])
#   slon_clat_clat = (sin(lonEdge[iEdge])/cos(latEdge[iEdge])) / cos(latEdge[iEdge])
#   clon_clat_clat = (cos(lonEdge[iEdge])/cos(latEdge[iEdge])) / cos(latEdge[iEdge])
#   v[iEdge] =  clon_clat * del2u[iEdge] - slon_clat * tangentialVelocity[iEdge]
#   u[iEdge] = -slon_clat_clat * del2u[iEdge] - clon_clat_clat * tangentialVelocity[iEdge]
    dvEdgeInv = 1.0 / dvEdge[iEdge]

    u[iEdge] = del2u[iEdge] * cos(angleEdge[iEdge]) - tangentialVelocity[iEdge]*sin(angleEdge[iEdge])
    v[iEdge] = del2u[iEdge] * sin(angleEdge[iEdge]) + tangentialVelocity[iEdge]*cos(angleEdge[iEdge])

    #u[iEdge] = u[iEdge] / dvEdgeInv
    #uLap_anl[iEdge] = uLap_anl[iEdge] / dvEdgeInv

end


# RMSE
anl = zeros(Float64, nEdges)
#Udiff = zeros(Float64, nEdges)
#Vdiff = zeros(Float64, nEdges)

global div0 = 0.0
global div1 = 0.0
for iCell in 1:nCells
    global div0 = div0 + (divergence[iCell]- 0.0)^2
end
#println("div   "," ",sqrt(div0)/nCells)

#println(minimum(divergence),' ',maximum(divergence))

global Usum0 = 0.0
global Usum1 = 0.0
global Vsum0 = 0.0
global Vsum1 = 0.0
global rvor0 = 0.0
global rvor1 = 0.0

global norm0 = 0.0
global norm1 = 0.0
global tanv0 = 0.0
global tanv1 = 0.0
for iEdge in 1:nEdges
#   if ( latEdge[iEdge] * 180.0/pi > -60.0 && latEdge[iEdge]*180.0/pi < 60.0 )
       global Usum0 = Usum0 + (u[iEdge] - uLap_anl[iEdge]) ^2.0
       global Usum1 = Usum1 + (           uLap_anl[iEdge]) ^2.0
       global Vsum0 = Vsum0 + (v[iEdge] - vLap_anl[iEdge]) ^2.0
       global Vsum1 = Vsum1 + (           vLap_anl[iEdge]) ^2.0
       global norm0 = norm0 + (del2u[iEdge] - normalVelocityAnl[iEdge]) ^2.0
       global norm1 = norm1 + (               normalVelocityAnl[iEdge]) ^2.0
       global tanv0 = tanv0 + (tangentialVelocity[iEdge] - tangentialVelocityAnl[iEdge]) ^2.0
       global tanv1 = tanv1 + (               tangentialVelocityAnl[iEdge]) ^2.0
   
       global rvor0 = rvor0 + (relativeVorticityEdge[iEdge] - relativeVorticityAnlEdge[iEdge]) ^2.0
       global rvor1 = rvor1 + (                               relativeVorticityAnlEdge[iEdge]) ^2.0
#   end
end
#println(minimum(relativeVorticityEdge),' ',maximum(relativeVorticityEdge))

println("NormV  "," ",sqrt(norm0/norm1))
println("TangV  "," ",sqrt(tanv0/tanv1))
println("Usph   "," ",sqrt(Usum0/Usum1))
println("Vsph   "," ",sqrt(Vsum0/Vsum1))
println("div    "," ",sqrt(div0/nCells))
println("relVor "," ",sqrt(rvor0/nEdges))
#println("relVor "," ",sqrt(rvor0/rvor1))




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
nc_var = defVar(ds_out,"uLapAnl",Float64,("nEdges",))
nc_var[:] = uLap_anl[:]
nc_var = defVar(ds_out,"vLapAnl",Float64,("nEdges",))
nc_var[:] = vLap_anl[:]
nc_var = defVar(ds_out,"normalVelocity",Float64,("nEdges",))
nc_var[:] = normalVelocity[:]
nc_var = defVar(ds_out,"divergence",Float64,("nCells",))
nc_var[:] = divergence[:]
nc_var = defVar(ds_out,"relativeVorticity",Float64,("nVertices",))
nc_var[:] = relativeVorticity[:]
nc_var = defVar(ds_out,"relativeVorticityAnlEdge",Float64,("nEdges",))
nc_var[:] = relativeVorticityAnlEdge[:]
nc_var = defVar(ds_out,"relativeVorticityEdge",Float64,("nEdges",))
nc_var[:] = relativeVorticityEdge[:]
nc_var = defVar(ds_out,"divergenceVert",Float64,("nVertices",))
nc_var[:] = divergenceVert[:]
nc_var = defVar(ds_out,"divergenceEdge",Float64,("nEdges",))
nc_var[:] = divergenceEdge[:]
nc_var = defVar(ds_out,"anl",Float64,("nEdges",))
nc_var[:] = anl[:]
nc_var = defVar(ds_out,"del2u",Float64,("nEdges",))
nc_var[:] = del2u[:]
nc_var = defVar(ds_out,"normalVelocityAnl",Float64,("nEdges",))
nc_var[:] = normalVelocityAnl[:]
nc_var = defVar(ds_out,"tangentialVelocityAnl",Float64,("nEdges",))
nc_var[:] = tangentialVelocityAnl[:]

close(ds_out)

###########################################################
