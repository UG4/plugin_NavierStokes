################################################################################
# NAVIER_STOKES
#
# The navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NAVIER_STOKES)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	navier_stokes_base.cpp
			fv1/navier_stokes_fv1.cpp
			fe/navier_stokes_fe.cpp
			fvcr/navier_stokes_fvcr.cpp
			fv/navier_stokes_fv.cpp

			bnd/no_normal_stress_outflow_base.cpp
			fv1/bnd/no_normal_stress_outflow_fv1.cpp
			fvcr/bnd/no_normal_stress_outflow_fvcr.cpp

			fvcr/bnd/symmetric_boundary.cpp
			fvcr/cr_reorder.cpp
			
			upwind.cpp
			
			fv1/register_fv1.cpp
			fv/register_fv.cpp
			fvcr/register_fvcr.cpp
			fe/register_fe.cpp
			navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(NavierStokes SHARED ${SOURCES})
	target_link_libraries (NavierStokes ug4)
endif(buildEmbeddedPlugins)
