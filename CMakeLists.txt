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

#Define a varible that can be used in TreeFrog to know the version
add_definitions( -DTREEFROGVERSION=${TF_VERSION})

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

# Precison options
tf_option(UNSIGNED_IDS "Use unsigned ID's in TreeFrog" OFF)
tf_option(LONG_INT "Use long ints ID's in TreeFrog, this will be set by NBODY_LONG_INT if not explicity set" ON)

#ID related options
tf_option(HALOIDNOTINDEX "If the haloID does not correspond to an index which it is located in the snapshot, if so this should be switched on and a index ordered, temporally unique ID will be outputted along with original ID" OFF)

#Other options
tf_option(BUILD_VELOCIRAPTOR "Whether VELOCIraptor should be built together with treefrog or not" OFF)


set(TF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(TF_DEFINES "")
set(TF_LIBS "")
set(TF_CXX_FLAGS "")
set(TF_LINK_FLAGS "")


#
# How we find HDF5 and set it up
#

macro(find_hdf5)
        # FindHDF5 needs an environment variable, oddly, unlike
        # most other packages that use normal cmake variables
        if (HDF5_ROOT)
                set(ENV{HDF5_ROOT} ${HDF5_ROOT})
        endif()
        find_package(HDF5 COMPONENTS CXX)
        if (HDF5_FOUND)
                list(APPEND TF_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
                list(APPEND TF_LIBS ${HDF5_LIBRARIES})
                list(APPEND TF_DEFINES USEHDF)
                set(TF_HAS_HDF5 Yes)
                if (NOT HDF5_VERSION VERSION_LESS "1.10.1")
                        list(APPEND TF_DEFINES "HDF5_NEWER_THAN_1_10_0")
                endif()
        endif()
endmacro()

#
# How we find MPI and set it up
#
macro(find_mpi)
        find_package(MPI)
        if (MPI_FOUND)
                list(APPEND TF_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
                list(APPEND TF_LIBS ${MPI_CXX_LIBRARIES})
                list(APPEND TF_CXX_FLAGS ${MPI_CXX_FLAGS})
                list(APPEND TF_LINK_FLAGS ${MPI_CXX_FLAGS})
                list(APPEND TF_DEFINES USEMPI)
                set(TF_HAS_MPI Yes)
        endif()
endmacro()

#
# Make sure we have the git submodules we need
#
macro(ensure_git_submodules)
        if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/NBodylib/CMakeLists.txt" OR NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/velociraptor_python_tools.py")
                find_package(Git QUIET)
                if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
                        # Update submodules as needed
                        message(STATUS "Updating NBodylib and tools submodule")
                        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
                                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
                        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
                        endif()
                else()
                        message(FATAL_ERROR "Cannot get NBodylib submodule or tools submodule automatically.
                                  Make sure you get your submodules")
                endif()
        endif()
endmacro()

#
# Go and find our dependencies, depending on whether
# we actually need them or not
#

ensure_git_submodules()

set(TF_HAS_HDF5 No)
if (TF_HDF5)
        find_hdf5()
endif()

set(TF_HAS_MPI No)
if (TF_MPI)
        find_mpi()
endif()


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

#
# Find what has been defined in NBodylib thats relevant to TreeFrog
#
list (FIND TF_DEFINES USEOPENMP _index)
if (${_index} GREATER -1)
      set(TF_HAS_OPENMP Yes)
else()
      set(TF_HAS_OPENMP No)
endif()

#Only if the treefrog longID's have been explicity set
if(NOT TF_LONG_INT)
        list (FIND TF_DEFINES LONGINT _index)
        if (${_index} GREATER -1)
              set(TF_HAS_LONG_INT Yes)
        else()
              set(TF_HAS_LONG_INT No)
        endif()
else()
        set(TF_HAS_LONG_INT Yes)
endif()

# Lets set which ID type TreeFrog will be using from, long int, long int unsigned, unsigned int and int
if(TF_HAS_LONG_INT)
        if(TF_UNSIGNED_IDS)
                set(TF_HAS_UNSIGNED_IDS Yes)
                list(APPEND TF_DEFINES TREEFROGLONGUIDS)
        else()
                set(TF_HAS_UNSIGNED_IDS No)
                list(APPEND TF_DEFINES TREEFROGLONGIDS)
        endif()
else()
        if(TF_UNSIGNED_IDS)
                set(TF_HAS_UNSIGNED_IDS Yes)
        else()
                set(TF_HAS_UNSIGNED_IDS No)
                list(APPEND TF_DEFINES TREEFROGINTIDS)
        endif()
endif()

if(TF_HALOIDNOTINDEX)
    set(TF_HAS_HALOIDNOTINDEX Yes)
    list(APPEND TF_DEFINES HALOIDNOTINDEX)
else()
    set(TF_HAS_HALOIDNOTINDEX No)
endif()

#
# Tell the world what what we are doing
#
macro(tf_report feature)

        # Output feature name and underscore it in the next line
        message("\n${feature}")
        string(REGEX REPLACE "." "-" _underscores ${feature})
        message("${_underscores}\n")

        set(_args "${ARGN}")
        list(LENGTH _args _nargs)
        math(EXPR _nargs "${_nargs} - 1")
        foreach(_idx RANGE 0 ${_nargs} 2)

                # Items in the list come with a message first, then the variable name
                list(GET _args ${_idx} _msg)
                math(EXPR _idx2 "${_idx} + 1")
                list(GET _args ${_idx2} _varname)

                # We try to keep things up to 80 cols
                string(LENGTH ${_msg} _len)
                math(EXPR _nspaces "75 - ${_len}")
                string(RANDOM LENGTH ${_nspaces} _spaces)
                string(REGEX REPLACE "." " " _spaces "${_spaces}")
                string(CONCAT _msg "${_msg}" ${_spaces})
                message(" ${_msg} ${TF_HAS_${_varname}}")
        endforeach()
endmacro()

message("\nTreeFrog successfully configured with the following settings:")
tf_report("Dependencies"
          "OpenMP" OPENMP
          "MPI support" MPI)
tf_report("File formats"
          "HDF5" HDF5)
tf_report("Types"
          "Long particle IDs" LONG_INT
          "Unsigned particles IDs" UNSIGNED_IDS)
tf_report("ID related"
          "Halo ID does not correspond to index within snapshot" HALOIDNOTINDEX)
# tf_report("Others"
#           "VELOCIraptor" TF_BUILD_VELOCIRAPTOR)
message("")
message("Compilation")
message("-----------")
message("")
message(" Include dirs: ${TF_INCLUDE_DIRS}")
message(" Macros defined: ${TF_DEFINES}")
message(" Libs: ${TF_LIBS}")
message(" C++ flags: ${TF_CXX_FLAGS}")
message(" Link flags: ${TF_LINK_FLAGS}")
message("")

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
