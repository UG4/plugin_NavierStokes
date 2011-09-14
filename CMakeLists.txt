# created by Sebastian Reiter

cmake_minimum_required(VERSION 2.6)

################################################
# SAMPLE_PLUGIN
################################################

project(NAVIER_STOKES_PLUGIN)

# Plugins only work in a shared environment
#set(UG_SHARED ON)
#set(UG_PLUGINS ON)

# include the definitions and dependencies from ug4
include("../../../ug_cmake_includes.txt")

# the plugin shall be written to the plugins folder.
message(STATUS "Info: Plugin will be written to bin/plugins")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/plugins)

# sources
set(SOURCES	navier_stokes_fv1.cpp
			navier_stokes_plugin.cpp)

add_library(NavierStokesPlugin SHARED ${SOURCES})
target_link_libraries (NavierStokesPlugin ug4)