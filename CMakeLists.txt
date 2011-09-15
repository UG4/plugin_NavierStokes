################################################################################
# NAVIER_STOKES_PLUGIN
#
# The navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(NAVIER_STOKES_PLUGIN)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	navier_stokes_fv1.cpp
			navier_stokes_plugin.cpp)


add_library(NavierStokesPlugin SHARED ${SOURCES})
target_link_libraries (NavierStokesPlugin ug4)
