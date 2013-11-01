################################################################################
# COMPRESSIBLE_NAVIER_STOKES
#
# The compressible-navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(COMPRESSIBLE_NAVIER_STOKES)

# include the definitions and dependencies for ug-plugins.
include("../../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	../navier_stokes_base.cpp
			../upwind_interface.cpp
			../upwind.cpp			
			
			
			compressible_navier_stokes_base.cpp
			
			fv1/compressible_navier_stokes_fv1.cpp
			fv1/register_fv1.cpp
			
			register_compressible.cpp
			
			
			../navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes/compressible ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(CompressibleNavierStokes SHARED ${SOURCES})
	target_link_libraries (CompressibleNavierStokes ug4)
endif(buildEmbeddedPlugins)
