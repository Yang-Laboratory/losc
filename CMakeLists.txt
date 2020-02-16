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

## set name of this project.
#set(project_name losc)

#project(${project_name} C CXX)

#file(GLOB Headers
#    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
#    "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
#file(GLOB Sources
#    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
#    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
#    "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

enable_testing()

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# copied from google-test readme.md. See
# https://github.com/google/googletest/blob/master/googletest/README.md.
# {{{
# Download and unpack googletest at configure time
configure_file(cmake/cmake_gtest.in ${CMAKE_CURRENT_BINARY_DIR}/download/googletest/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/googletest )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/googletest )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/googletest)

# copied from google-test readme.md. See
# https://github.com/google/googletest/blob/master/googletest/README.md.
# {{{
# Download and unpack matrix at configure time
configure_file(cmake/cmake_matrix.in ${CMAKE_CURRENT_BINARY_DIR}/download/matrix/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/matrix)
if(result)
  message(FATAL_ERROR "CMake step for libmatrix failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/matrix)
if(result)
  message(FATAL_ERROR "Build step for libmatrix failed: ${result}")
endif()
# Add libmatrix directly to our build. This defines the `matrix` targets.
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/matrix)

# build localization library
add_subdirectory(losc)

#add_library(
#    ${project_name} SHARED
#    ${Sources}
#    ${Headers})
#
#target_link_libraries(
#    ${project_name}
#    PUBLIC
#    localization)
