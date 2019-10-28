# CMakeLists.txt for treefrog
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#

set(TF_SOURCES
    allvars.cxx
    crosscheck.cxx
    idroutines.cxx
    io.cxx
    mpiroutines.cxx
    otherio.cxx
    stfio.cxx
    ui.cxx
    utilities.cxx)

#build the halo merger tree library.
add_library(hmt STATIC ${TF_SOURCES})
target_compile_definitions(hmt PRIVATE ${TF_DEFINES})
if (TF_CXX_FLAGS)
	set_target_properties(hmt PROPERTIES COMPILE_FLAGS ${TF_CXX_FLAGS})
endif()
set_target_properties(hmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Bringing NBodyLib objects into the velociraptor library
target_link_libraries(hmt nbodylib_iface)

#build the treefrog executable.
add_executable(treefrog main.cxx)
target_link_libraries(treefrog hmt ${TF_LIBS})
target_compile_definitions(treefrog PRIVATE ${TF_DEFINES})

if (TF_LINK_FLAGS)
	set_target_properties(treefrog PROPERTIES LINK_FLAGS ${TF_LINK_FLAGS})
endif()
if (TF_CXX_FLAGS)
    set_target_properties(treefrog PROPERTIES COMPILE_FLAGS ${TF_CXX_FLAGS})
endif()
set_target_properties(treefrog PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
