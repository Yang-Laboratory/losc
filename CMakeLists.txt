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

# set name of this project.
set(project_name losc)

project(${project_name} C CXX)

file(GLOB Headers
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB Sources
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

enable_testing()

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# build google test library
#add_subdirectory(googletest)
if(IS_DIRECTORY $ENV{HOME}/local/googletest)
    message(STATUS "found googletest source folder.")
    add_subdirectory($ENV{HOME}/local/googletest $ENV{HOME}/local/googletest/build)
else()
    message(FATAL_ERROR "googletest source folder not found.
        Run `git clone https://github.com/google/googletest.git` under path `$ENV{HOME}/local/`")
endif()

# build matrix library
add_subdirectory(matrix)

# build localization library
add_subdirectory(localization)

add_library(
    ${project_name} SHARED
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name}
    PUBLIC
    localization)
