cmake_minimum_required (VERSION 3.6)

project (MPUI)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})

set(LIB_DIR ${PROJECT_SOURCE_DIR}/libs)

# GLM
include_directories(${LIB_DIR}/glm)

# Nanoflann
include_directories(${LIB_DIR}/Nanoflann)

# Eigen
include_directories(${LIB_DIR}/eigen)

# Folders considered for compilation
add_subdirectory(src)

####################################
####################################
####################################
# First we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# Check if Doxygen is installed
find_package(Doxygen)

if (DOXYGEN_FOUND)
    # Set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc_doxygen)

    # Note the option ALL which allows to build the docs together with the application
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)