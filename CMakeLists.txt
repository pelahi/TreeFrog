# CMakeLists.txt for treefrog
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

cmake_minimum_required(VERSION 3.0)

# We have the version here
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION TF_VERSION)
string(REPLACE "." ";" VERSION_LIST ${TF_VERSION})
list(GET VERSION_LIST 0 TF_VERSION_MAJOR)
list(GET VERSION_LIST 1 TF_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(TF_LANGS C CXX)
else()
	set(TF_LANGS CXX)
endif()

project(TreeFrog VERSION ${TF_VERSION} LANGUAGES ${TF_LANGS})
set(CMAKE_CXX_STANDARD 11)

# Command-line treefrog-specific options
#
# Options users can give on the command line via -D
#
macro(tf_option optname optdesc status)
        option(TF_${optname} "${optdesc}" "${status}")
endmacro()

# Input formats
tf_option(HDF5 "Attempt to include HDF5 support in TreeFrog " ON)

# MPI options
tf_option(MPI               "Attempt to include MPI support in TreeFrog" ON)

tf_option(TF_BUILD_VELOCIRAPTOR "Whether VELOCiraptor should be built together with treefrog or not" OFF)

set(TF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(TF_DEFINES "")
set(TF_LIBS "")
set(TF_CXX_FLAGS "")
set(TF_LINK_FLAGS "")

# Include nbodylib, which is definitely needed, and check its version
add_subdirectory(NBodylib NBodylib)
if (NBODYLIB_VERSION VERSION_LESS "1.10")
	message(FATAL_ERROR "NBodyLib version ${NBODYLIB_VERSION} unsupported, treefrog requires >= 1.10")
endif()

list(APPEND TF_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
list(APPEND TF_DEFINES ${NBODYLIB_DEFINES})
list(APPEND TF_LIBS ${NBODYLIB_LIBS})
list(APPEND TF_CXX_FLAGS ${NBODYLIB_CXX_FLAGS})
list(APPEND TF_LINK_FLAGS ${NBODYLIB_LINK_FLAGS})
list(APPEND TF_DOC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/doc")

# Before had option to build velociraptor as a submodule but keeping the repos separate

list(INSERT TF_DOC_DIRS 0 ${NBODYLIB_DOC_DIRS})

# cleanup duplicates and convert to single string
list(REMOVE_DUPLICATES TF_CXX_FLAGS)
string(REPLACE ";" " " TF_CXX_FLAGS "${TF_CXX_FLAGS}")
list(REMOVE_DUPLICATES TF_LINK_FLAGS)
string(REPLACE ";" " " TF_LINK_FLAGS "${TF_LINK_FLAGS}")

include_directories(${TF_INCLUDE_DIRS})

# This provides us with the velociraptor library and the stf binary
add_subdirectory(src)

# Export the include directories, if necessary
# If building on our own, add the "doc" target
if (_export)
        set(TF_VERSION "${TF_VERSION_MAJOR}.${TF_VERSION_MINOR}" PARENT_SCOPE)
        set(TF_VERSION_MAJOR "${TF_VERSION_MAJOR}" PARENT_SCOPE)
        set(TF_VERSION_MINOR "${TF_VERSION_MINOR}" PARENT_SCOPE)
        set(TF_INCLUDE_DIRS "${TF_INCLUDE_DIRS}" PARENT_SCOPE)
        set(TF_DEFINES "${TF_DEFINES}" PARENT_SCOPE)
        set(TF_CXX_FLAGS "${TF_CXX_FLAGS}" PARENT_SCOPE)
        set(TF_LINK_FLAGS "${TF_LINK_FLAGS}" PARENT_SCOPE)
        set(TF_LIBS "${TF_LIBS}" PARENT_SCOPE)
        set(TF_DOC_DIRS "${TF_DOC_DIRS}" PARENT_SCOPE)
else()
        try_add_doc_target("${TF_DOC_DIRS}")
endif()

