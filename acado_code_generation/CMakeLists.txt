#
# Project settings
#

# Minimum required version of cmake
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

# Project name and programming languages used
PROJECT( code_generation_ACADO_for_f110_nmpc CXX )

# CMake module(s) path
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR} )

#
# Prerequisites
#
FIND_PACKAGE( ACADO REQUIRED )

#
# Include directories
#
INCLUDE_DIRECTORIES( . ${ACADO_INCLUDE_DIRS} )

FILE(GLOB SOURCES "src/formulation.cpp")

#
# Build an executable
#
ADD_EXECUTABLE( main ${SOURCES} )

# Linux and OS X
TARGET_LINK_LIBRARIES( main ${ACADO_SHARED_LIBRARIES} )
