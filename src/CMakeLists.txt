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

set(_treefrog_sources
    allvars.cxx
    crosscheck.cxx
    idroutines.cxx
    io.cxx
    main.cxx
    mpiroutines.cxx
    otherio.cxx
    stfio.cxx
    ui.cxx
    utilities.cxx)

add_executable(treefrog ${_treefrog_sources})
target_link_libraries(treefrog nbodylib ${TF_LIBS})
target_compile_definitions(treefrog PRIVATE ${TF_DEFINES})
if (TF_CXX_FLAGS)
	set_target_properties(treefrog PROPERTIES COMPILE_FLAGS ${TF_CXX_FLAGS})
endif()
if (TF_LINK_FLAGS)
	set_target_properties(treefrog PROPERTIES LINK_FLAGS ${TF_LINK_FLAGS})
endif()
set_target_properties(treefrog PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
