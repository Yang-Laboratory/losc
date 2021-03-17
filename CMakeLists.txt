cmake_minimum_required(VERSION 3.8)

# names of projects
set(PROJECT_LOSC losc)
project(${PROJECT_LOSC} C CXX)

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Download pybind11 library
set(PYBIND11_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/pybind11")
if(NOT EXISTS "${PYBIND11_SOURCE_DIR}")
    message(STATUS "Downloading external library: pybind11.")
    configure_file(cmake/download_pybind11.in download/pybind11/CMakeLists.txt)
    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/pybind11"
    )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/pybind11"
    )
endif()

# find Eigen3 library.
find_package (Eigen3 NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 library found.")
else()
    # Download and build Eigen3
    message(STATUS "Downloading external library: Eigen3.")
    set(EXTERNAL_INSTALL_LOCATION "${CMAKE_SOURCE_DIR}/external")
    configure_file(cmake/download_eigen3.in download/eigen3/CMakeLists.txt)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/eigen3"
    )
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/eigen3"
    )
    # find Eigen3 again
    find_package (Eigen3 REQUIRED NO_MODULE)
    if(Eigen3_FOUND)
        message(STATUS "External Eigen3 library is installed.")
    else()
        message(FATAL_ERROR "Eigen3 library not found. Please install Eigen3 library.")
    endif()
endif()

# ==> build pybind11 <==
add_subdirectory("${PYBIND11_SOURCE_DIR}")

# ==> build source <==
add_subdirectory(src)

# ==> build tests <==
option(BUILD_TEST "Build test or not." OFF)
if (BUILD_TEST)
    enable_testing()
    # Download googletest library
    set(GOOGLE_TEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/googletest")
    if(NOT EXISTS "${GOOGLE_TEST_SOURCE_DIR}")
        message(STATUS "Downloading external library: googletest.")
        configure_file(cmake/download_gtest.in download/googletest/CMakeLists.txt)
        execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/googletest"
        )
        execute_process(COMMAND "${CMAKE_COMMAND}" --build .
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/googletest"
        )
    endif()

    # build google test
    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory("${GOOGLE_TEST_SOURCE_DIR}")

    # build test source
    add_subdirectory(tests)
else()
    message(STATUS "Note: Tests are disabled. Use \"-DBUILD_TEST=On\" to build tests.")
endif()

# ==> build docs <==
option(BUILD_DOC "Build documentation for Losc library" OFF)
if (BUILD_DOC)
    add_subdirectory(docs_source)
else (BUILD_DOC)
    message("Note: If you want to build documentation, use `-DBUILD_DOC=On"
            "Doxygen will be required.")
endif (BUILD_DOC)
