#!/bin/sh

#----
julia mpas_sphere_full_Legendre_edge_del4.jl
python display_scatter_4.py
mv fig_scatter.png fig_scatter_del4_default.png

#----
julia mpas_sphere_full_Legendre_edge_del4_vorintp.jl
python display_scatter_4.py
mv fig_scatter.png fig_scatter_del4_vorintp.png

#----
julia mpas_sphere_full_Legendre_edge_Tan1_del4.jl
python display_scatter_4.py
mv fig_scatter.png fig_scatter_del4_tan1.png

#----
julia mpas_sphere_full_Legendre_edge_Tan2_del4.jl
python display_scatter_4.py
mv fig_scatter.png fig_scatter_del4_tan2.png


#julia mpas_sphere_full_Legendre_edge_UsphTest.jl
#julia mpas_sphere_full_Legendre_edge_Vxyz.jl
#julia mpas_sphere_full_Legendre_edge_UsphTest_slopeLimiter.jl
#julia mpas_sphere_full_Legendre_edge_UsphTest_extAvg.jl
#python display_scatter.py

#python display_scatter_4.py
#display fig_scatter.png
