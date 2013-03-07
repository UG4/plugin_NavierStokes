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
set(SOURCES	navier_stokes.cpp
			upwind.cpp
			fvcr/cr_reorder.cpp
			bnd/no_normal_stress_outflow.cpp
			bnd/symmetric_boundary.cpp
			navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(NavierStokes SHARED ${SOURCES})
	target_link_libraries (NavierStokes ug4)
endif(buildEmbeddedPlugins)
