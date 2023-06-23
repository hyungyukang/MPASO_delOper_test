using NCDatasets
#using Plots
using PyPlot
import OffsetArrays.OffsetArray
using SHTOOLS

nxgrid = "10242"

    if ( nxgrid == "2562" )
         # dx = 480 km ~ 4deg
         const WN = 40
         const WN1 = WN + 1
elseif ( nxgrid == "10242" )
         # dx = 240 km ~ 2deg
         const WN = 80
         const WN1 = WN + 1
elseif ( nxgrid == "40962" )
         # dx = 120 km ~ 1deg
         const WN = 120
         const WN1 = WN + 1
elseif ( nxgrid == "655362" )
         # dx = 60 km ~ 1deg
         const WN = 240
         const WN1 = WN + 1
elseif ( nxgrid == "2621442" )
         # dx = 30 km ~ 1deg
         const WN = 480
         const WN1 = WN + 1
elseif ( nxgrid == "10485762" )
         # dx = 15 km ~ 1deg
         const WN = 960
         const WN1 = WN + 1
end

###########################################################

mesh_file_name = "../../meshes/x1."*nxgrid*".grid.nc"   # Miniweather output file name
#mesh_file_name = "./x1."*nxgrid*".grid.nc"   # Miniweather output file name
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

nCells = ds.dim["nCells"]
nEdges = ds.dim["nEdges"]
nVertices = ds.dim["nVertices"]
maxEdges = ds.dim["maxEdges"]
maxEdges2 = ds.dim["maxEdges2"]
vertexDegree = ds.dim["vertexDegree"]

edgeSignOnCell = zeros(Float64, maxEdges,nCells)
areaEdge = zeros(Float64, nEdges)

###########################################################
 
for iCell in 1:nCells
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        
        if ( iCell == cellsOnEdge[1,iEdge])
           edgeSignOnCell[i,iCell] = -1.0
        else
           edgeSignOnCell[i,iCell] =  1.0
        end

    end
    #areaCell[iCell] = areaCell[iCell] * 6371000^2
end 

for iEdge in 1:nEdges
    areaEdge[iEdge] = dcEdge[iEdge] * dvEdge[iEdge]/2.0
#    #areaEdge[iEdge] = dcEdge[iEdge] * (dvEdge[iEdge]+dv2Edge[iEdge]+dv1Edge[iEdge])/4.0
end

###########################################################
###########################################################

# Allocation array
_alp = zeros(Float64,WN+1,WN+1,nCells)
alp  = OffsetArray(_alp,0:WN,0:WN,1:nCells)
_dalp = zeros(Float64,WN+1,WN+1,nCells)
dalp  = OffsetArray(_dalp,0:WN,0:WN,1:nCells)
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
for iCell in 1:nCells
    z = sin(latCell[iCell])
    x = 1.0-z^2

    alp[ 0,0,iCell] = 1.0

    alp[ 0,1,iCell] = z
    alp[ 1,1,iCell] = -sqrt(x)

    alp[ 0,2,iCell] = 0.5*(3.0*z^2-1.0)
    alp[ 1,2,iCell] = -3.0*z*sqrt(x)
    alp[ 2,2,iCell] =  3.0*(x)

    alp[ 0,3,iCell] = 0.5*(5.0*z^3-3.0*z)
    alp[ 1,3,iCell] = 1.5*(1.0-5.0*z^2)*sqrt(x)
    alp[ 2,3,iCell] = 15.0*z*(x)
    alp[ 3,3,iCell] =-15.0*x^(3.0/2.0)
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

    for n = 0:WN
        for m = 0:n
            alp[m,n,iCell] = norm[m,n] * alp[m,n,iCell] #* sqrt(2.0)
        end
    end

end 


#########################################################

# Derivation of ALP
for iCell in 1:nCells
    z = sin(latCell[iCell])
    a = sqrt(1.0-z^2)

    m = 0
    for n in m+1:WN
#       dalp[m,n,iCell] = ( n*z*alp[m,n,iCell]-(n+m)*alp[m,n-1,iCell] ) / ((-a))

        # Derivative w.r.t [theta] not [x]
        dalp[m,n,iCell] =  (n*z*alp[m,n,iCell] - (n+m)*alp[m,n-1,iCell])
    end

    for n in 1:WN
        for m in 1:n
#           dalp[m,n,iCell] = ( -(n+m)*(n-m+1)*a*alp[m-1,n,iCell]
#                               -m*z*alp[m,n,iCell] ) / ((-a))

            # Derivative w.r.t [theta] not [x]
            dalp[m,n,iCell] =  ( -(n+m)*(n-m+1)*a*alp[m-1,n,iCell]
                                -m*z*alp[m,n,iCell] )
        end
    end
end

# Save the last mode
init = zeros(Float64, nCells)
anl = zeros(Float64, nCells)
diff = zeros(Float64, nCells)
diffEdge1 = zeros(Float64, nEdges)
diffEdge2 = zeros(Float64, nEdges)
num = zeros(Float64, nCells)
uCell = zeros(Float64, nCells)
vCell = zeros(Float64, nCells)
uEdge = zeros(Float64, nEdges)
vEdge = zeros(Float64, nEdges)
uEdge_anl = zeros(Float64, nEdges)
vEdge_anl = zeros(Float64, nEdges)
normalVelocity = zeros(Float64, nEdges)
tangentialVelocity = zeros(Float64, nEdges)
normalVelocity_anl = zeros(Float64, nEdges)
angentialVelocity_anl = zeros(Float64, nEdges)


#m = Integer(WN/2)
#m = 5
#n = 10
#m = 0
for n in 5:5:80

    m = n# floor(Int8,n/2.0)

for iCell in 1:nCells
    a =      exp(im*m*lonCell[iCell]) * alp[m,n,iCell] 
    b = im*m*exp(im*m*lonCell[iCell]) * alp[m,n,iCell] 
    c =      exp(im*m*lonCell[iCell]) * ( dalp[m,n,iCell] 
                                         -sin(latCell[iCell])* alp[m,n,iCell])
    init[iCell] = a.re

    
    #uCell[iCell] = -c.re 
    #vCell[iCell] =  b.re

    uCell[iCell] = a.re

    # To check probability
    anl[iCell] = (abs(a)^2.0) * areaCell[iCell]

end
#println("Check integ(|Ynm|^2) : ",sum(anl)) # This must be unity with some numerical errors


# Divergence
for iCell in 1:nCells
    div = 0
    for i in 1:nEdgesOnCell[iCell]
        iEdge = edgesOnCell[i,iCell]
        cell1 = cellsOnEdge[1,iEdge]
        cell2 = cellsOnEdge[2,iEdge]

        # grad
        grad = - (init[cell2] - init[cell1]) / dcEdge[iEdge]
 
        # div
        num[iCell] = num[iCell] + edgeSignOnCell[i,iCell] * grad * dvEdge[iEdge]
    end

    num[iCell] = num[iCell] / areaCell[iCell]

    anl[iCell] = -n*(n+1)*init[iCell]
    diff[iCell] = (num[iCell]-anl[iCell])^2
    anl[iCell] = anl[iCell]^2
    
end

sum1 = sum(diff)
sum2 = sum(anl)

println(m,',',n,',',"Error = ",sqrt(sum1/sum2))

for iCell in 1:nCells
    anl[iCell] = 0.0
    num[iCell] = 0.0
end

end










## Interpolation - Cell => Edge
#for iEdge in 1:nEdges
#    cell1 = cellsOnEdge[1,iEdge]
#    cell2 = cellsOnEdge[2,iEdge]
#    uEdge_anl[iEdge] = (uCell[cell1]+uCell[cell2])/2.0
#    vEdge_anl[iEdge] = (vCell[cell1]+vCell[cell2])/2.0
#end
#
## u,v -> normalVelocity, tangentialVelocity
#for iEdge in 1:nEdges
#        normalVelocity_anl[iEdge] =  uEdge_anl[iEdge]*cos(angleEdge[iEdge]) + vEdge_anl[iEdge]*sin(angleEdge[iEdge])
#    tangentialVelocity_anl[iEdge] = -uEdge_anl[iEdge]*sin(angleEdge[iEdge]) + vEdge_anl[iEdge]*cos(angleEdge[iEdge])
#end
#
## Divergence
#for iCell in 1:nCells
#    div = 0.0
#    for i in 1:nEdgesOnCell[iCell]
#        iEdge = edgesOnCell[i,iCell]
#        # div
#        #r_tmp = dvEdge[iEdge]*normalVelocity[iEdge]*invAreaCell1
#        r_tmp = dvEdge[iEdge]*normalVelocity_anl[iEdge]
#        div = div - edgeSignOnCell[i,iCell] * r_tmp
#    end
#    divergence[iCell] = div / areaCell[iCell]
#end
#
## Relative Vorticity
#for iVertex in 1:nVertices
#    invAreaTri1 = 1.0 / areaTriangle[iVertex]
#    for i in 1:vertexDegree
#        iEdge = edgesOnVertex[i,iVertex]
#        #r_tmp = dcEdge[iEdge] * normalVelocity[iEdge]
#        r_tmp = dcEdge[iEdge] * normalVelocity_anl[iEdge]
#        relativeVorticity[iVertex] = relativeVorticity[iVertex] + edgeSignOnVertex[i,iVertex] * r_tmp * invAreaTri1
#    end
#end
#
#
#
## MPAS grad
#for iEdge in 1:nEdges
#    cell1 = cellsOnEdge[1,iEdge]
#    cell2 = cellsOnEdge[2,iEdge]
#    grad = - (init[cell2] - init[cell1]) / dcEdge[iEdge]
#    normalVelocity[iEdge] = grad / sqrt(3.0)
#end
#
#for iEdge in 1:nEdges
#    for i in 1:nEdgesOnEdge[iEdge]
#        eoe = edgesOnEdge[i,iEdge]
#        tangentialVelocity[iEdge] = tangentialVelocity[iEdge] + weightsOnEdge[i,iEdge] * normalVelocity[eoe]
#    end
#end
#
##print(tangentialVelocity[1:100])
#
## u,v -> normalVelocity, tangentialVelocity
#
###for iEdge in 1:nEdges
##for iEdge in 1:50
##        normalVelocity_anl[iEdge] =  uEdge[iEdge]*cos(angleEdge[iEdge]) + vEdge[iEdge]*sin(angleEdge[iEdge])
##    tangentialVelocity_anl[iEdge] = -uEdge[iEdge]*sin(angleEdge[iEdge]) + vEdge[iEdge]*cos(angleEdge[iEdge])
##
##    println(iEdge,',',normalVelocity[iEdge],',',normalVelocity_anl[iEdge])
##    diffEdge[iEdge] = (normalVelocity[iEdge]-normalVelocity_anl[iEdge])
##end
#
#
### normalVelocity, tangentialVelocity -> u,v
#for iEdge in 1:nEdges
#    # Zonal com => V spherical
#    vEdge[iEdge] = normalVelocity[iEdge] * cos(angleEdge[iEdge]) - tangentialVelocity[iEdge]*sin(angleEdge[iEdge])
#    
#    # Meridional comp => -U spherical
#    uEdge[iEdge] =-(normalVelocity[iEdge] * sin(angleEdge[iEdge]) + tangentialVelocity[iEdge]*cos(angleEdge[iEdge])) /
#                  (cos(latEdge[iEdge]) )#* sqrt(3.0))
#                                           # 1/root(3) for U wind based on Gaussman (2011)
#    diffEdge1[iEdge] = (uEdge[iEdge] - uEdge_anl[iEdge])^2
#    diffEdge2[iEdge] = uEdge_anl[iEdge]^2
#end
#
#for iEdge in 1:50
#    println(iEdge,',',uEdge[iEdge],',',uEdge_anl[iEdge])
#end
#
#sum1 = sum(diffEdge1)
#sum2 = sum(diffEdge2)
#
#println("U-wind error =", sqrt(sum1/sum2))


#####################################################################

#          normalVelocity(:,iEdge) = velocityZonal*cos(angleEdge(iEdge)) + velocityMeridional*sin(angleEdge(iEdge))
#          tangentialVelocity(:,iEdge) = -velocityZonal*sin(angleEdge(iEdge)) + velocityMeridional*cos(angleEdge(iEdge))

#            zonalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*cos(angleEdge(iEdge)) &
#              - tangentialBarotropicVel(iEdge)*sin(angleEdge(iEdge))
#            meridionalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*sin(angleEdge(iEdge)) &
#              + tangentialBarotropicVel(iEdge)*cos(angleEdge(iEdge))



############################################################
############################################################

# Output

#out_file_name = "./output.nc"
#ds_out = Dataset(out_file_name,"c")
#
#defDim(ds_out,"nCells",nCells)
#defDim(ds_out,"nEdges",nEdges)
#
#nc_var = defVar(ds_out,"lonCell",Float64,("nCells",))
#nc_var[:] = lonCell[:]
#nc_var = defVar(ds_out,"latCell",Float64,("nCells",))
#nc_var[:] = latCell[:]
#nc_var = defVar(ds_out,"init",Float64,("nCells",))
#nc_var[:] = init[:]
#nc_var = defVar(ds_out,"uEdge",Float64,("nEdges",))
#nc_var[:] = uEdge[:]
#nc_var = defVar(ds_out,"vEdge",Float64,("nEdges",))
#nc_var[:] = vEdge[:]
#nc_var = defVar(ds_out,"uEdge_anl",Float64,("nEdges",))
#nc_var[:] = uEdge_anl[:]
#nc_var = defVar(ds_out,"vEdge_anl",Float64,("nEdges",))
#nc_var[:] = vEdge_anl[:]
#nc_var = defVar(ds_out,"normalVelocity",Float64,("nEdges",))
#nc_var[:] = normalVelocity[:]
#
#close(ds_out)

#print(latCell)

###########################################################
