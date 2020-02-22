cmake_minimum_required(VERSION 3.8)

# build Losc documentation file at configure time.
# first we can indicate the documentation build as an option and set it to ON by
# default
option(BUILD_DOC "Build documentation for Losc library" ON)

if (BUILD_DOC)
    if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
        # check if Doxygen is installed
        find_package(Doxygen)
        if (DOXYGEN_FOUND)
            # set input and output files
            set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
            set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

            # request to configure the file
            configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
            message(STATUS "Doxygen build started")

            # note the option ALL which allows to build the docs together with the
            # application
            add_custom_target(doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )
        else (DOXYGEN_FOUND)
          message("Doxygen need to be installed to generate the doxygen documentation")
        endif (DOXYGEN_FOUND)
    else ()
        message(STATUS "Documentations need to be build under release compliation.")
    endif ()
else (BUILD_DOC)
   message("NOte: If you want to build documentation, use `-DBUILD_DOC=On -DCMAKE_BUILD_TYPE=Release`. "
    "Doxygen will be required.")
endif (BUILD_DOC)

# compiler setting: use intel C/C++ compiler.
find_program(ICC icc)
if (ICC)
    message(STATUS "icc found: ${ICC}")
    set(CMAKE_C_COMPILER "${ICC}" CACHE FILEPATH "Intel C compiler" FORCE)
else()
    message(WARNING "icc not found.")
endif()
find_program(ICPC icpc)
if (ICPC)
    message(STATUS "icpc found: ${ICPC}")
    set(CMAKE_CXX_COMPILER "${ICPC}" CACHE FILEPATH "Intel C++ compiler" FORCE)
else()
    message(WARNING "icpc not found.")
endif()

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(project_name losc)
project(${project_name} C CXX)

enable_testing()

file(GLOB Headers
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
file(GLOB Sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

# enable openmp or not.
if (OPENMP)
    message(STATUS "trying to use openmp.")
    add_compile_definitions(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "openmp for CXX is found and enabled.")
    else()
        message(WARNING "Fail to find openmp for CXX. Will compile WITHOUT openmp for CXX.")
    endif()
else()
    message(STATUS "openmp for CXX is not used. If you want to enable openmp, use `cmake -DOPENMP=on`.")
endif()

find_package (Eigen3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 library found.")
else()
    message(WARNING "Can not found library Eigen3, which will be used for testing.\
        Test code will be skipped until Eigen3 is installed.")
    return()
endif()

# blas library
find_library(
    BLAS_LIBRARY
    NAMES blas
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "libblas: ${BLAS_LIBRARY}")
else()
    message(FATAL_ERROR "libblas not found")
endif()

# lapack library
find_library(
    LAPACK_LIBRARY
    NAMES lapack
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "liblapack: ${LAPACK_LIBRARY}")
else()
    message(FATAL_ERROR "liblapack not found")
endif()


add_library(
    ${project_name} SHARED
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name}
    PUBLIC
    ${BLAS_LIBRARY}
    ${LAPAK_LIBRARY}
    Eigen3::Eigen)

if (OPENMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(
            ${project_name}
            PUBLIC
            OpenMP::OpenMP_CXX)
    endif()
endif()

target_include_directories(
    ${project_name}
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(test)
