################################################################################
# INCOMPRESSIBLE_NAVIER_STOKES
#
# The incompressible-navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(INCOMPRESSIBLE_NAVIER_STOKES)

# include the definitions and dependencies for ug-plugins.
include("../../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	../navier_stokes_base.cpp
			../upwind_interface.cpp
			../upwind.cpp
			
			
			incompressible_navier_stokes_base.cpp	
			
			fv1/navier_stokes_fv1.cpp
			fe/navier_stokes_fe.cpp
			fvcr/navier_stokes_fvcr.cpp
			fv/navier_stokes_fv.cpp

			bnd/no_normal_stress_outflow_base.cpp
			fv1/bnd/no_normal_stress_outflow_fv1.cpp
			fv/bnd/no_normal_stress_outflow_fv.cpp
			fvcr/bnd/no_normal_stress_outflow_fvcr.cpp

			fvcr/bnd/symmetric_boundary_fvcr.cpp
			fv1/bnd/symmetric_boundary_fv1.cpp

			fvcr/cr_reorder.cpp

			fv1/stabilization.cpp
			
			fv1/register_fv1.cpp
			fv/register_fv.cpp
			fvcr/register_fvcr.cpp
			fe/register_fe.cpp	
			
			../register_navier_stokes.cpp
			incompressible_navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes/incompressible ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(IncompressibleNavierStokes SHARED ${SOURCES})
	target_link_libraries (IncompressibleNavierStokes ug4)
endif(buildEmbeddedPlugins)