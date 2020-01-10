cmake_minimum_required(VERSION 3.8)

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

set(project_name localization)
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

add_library(
    ${project_name} SHARED
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name}
    PUBLIC
    matrix)

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
