################################################################################
# NAVIER_STOKES
#
# The navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NAVIER_STOKES)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	navier_stokes_base.cpp
			upwind_interface.cpp
			upwind.cpp
			
			
			incompressible/incompressible_navier_stokes_base.cpp	
			
			incompressible/fv1/navier_stokes_fv1.cpp
			incompressible/fe/navier_stokes_fe.cpp
			incompressible/fvcr/navier_stokes_fvcr.cpp
			incompressible/fv/navier_stokes_fv.cpp

			incompressible/bnd/no_normal_stress_outflow_base.cpp
			incompressible/fv1/bnd/no_normal_stress_outflow_fv1.cpp
			incompressible/fv/bnd/no_normal_stress_outflow_fv.cpp
			incompressible/fvcr/bnd/no_normal_stress_outflow_fvcr.cpp

			incompressible/fvcr/bnd/symmetric_boundary_fvcr.cpp
			incompressible/fv1/bnd/symmetric_boundary_fv1.cpp

			incompressible/fvcr/cr_reorder.cpp

			incompressible/fv1/stabilization.cpp
			
			incompressible/fv1/register_fv1.cpp
			incompressible/fv/register_fv.cpp
			incompressible/fvcr/register_fvcr.cpp
			incompressible/fe/register_fe.cpp	
			
			incompressible/register_incompressible.cpp
			
			
			compressible/compressible_navier_stokes_base.cpp
			
			compressible/fv1/compressible_navier_stokes_fv1.cpp
			compressible/fv1/register_fv1.cpp
			
			compressible/register_compressible.cpp
			
			
			navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(NavierStokes SHARED ${SOURCES})
	target_link_libraries (NavierStokes ug4)
endif(buildEmbeddedPlugins)
