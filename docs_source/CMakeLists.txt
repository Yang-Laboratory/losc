# To use Sphinx + Breathe + Doxygen to generate documentations
# for C/C++ project. Refer to this blog for an example:
# https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/.
#
# I think doxygen generates the documentations for C/C++ in a better look
# compared to sphinx. So I switch back to using doxygen.

# Add the cmake folder, so the FindSphinx module can be found.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# require doxygen and sphinx is installed.
find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# some variables used to configure sphinx source files.
set(LOSC_DOXYGEN_HTML_OUTPUT "")
set(C_LOSC_DOXYGEN_HTML_OUTPUT "")

# ==> Doxygen <==
macro(build_doxygen module)
    set(DOXYGEN_INPUT_DIR ${CMAKE_SOURCE_DIR}/include/${module})
    set(DOXYGEN_SOURCE_ROOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen/${module})
    set(DOXYGEN_OUTPUT_DIR ${DOXYGEN_SOURCE_ROOT_DIR})

    set(DOXYGEN_HTML_OUTPUT_DIR ${DOXYGEN_SOURCE_ROOT_DIR}/html)
    set(DOXYGEN_HTML_OUTPUT ${DOXYGEN_HTML_OUTPUT_DIR}/index.html)

    if ("${module}" STREQUAL "losc")
        set(DOXYGEN_PROJECT_NAME "\"Localized Orbital Scaling Correction (LOSC) C++ library\"")
        set(LOSC_DOXYGEN_HTML_OUTPUT ${DOXYGEN_HTML_OUTPUT})
    elseif ("${module}" STREQUAL "c_losc")
        set(DOXYGEN_PROJECT_NAME "\"Localized Orbital Scaling Correction (LOSC) C library\"")
        set(C_LOSC_DOXYGEN_HTML_OUTPUT ${DOXYGEN_HTML_OUTPUT})
    else()
        message(FATAL_ERROR "Unknown module to build docs: ${module}")
    endif()
    message(STATUS "Build docs for module: ${module}")

    # create some directories that Doxygen won't create for us
    file(MAKE_DIRECTORY ${DOXYGEN_SOURCE_ROOT_DIR})
    file(MAKE_DIRECTORY ${DOXYGEN_HTML_OUTPUT_DIR})

    # configure the configuration file for doxygen.
    file(GLOB DOXYGEN_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/*)
    file(COPY ${DOXYGEN_SOURCE_FILES} DESTINATION ${DOXYGEN_SOURCE_ROOT_DIR})
    set(DOXYFILE_IN ${DOXYGEN_SOURCE_ROOT_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${DOXYGEN_SOURCE_ROOT_DIR}/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    # print some informational message.
    message(STATUS "${module}: doxygen configuration file: ${DOXYFILE_OUT}")
    message(STATUS "${module}: DOXYGEN_INPUT_DIR: ${DOXYGEN_INPUT_DIR}")
    message(STATUS "${module}: DOXYGEN_SOURCE_ROOT_DIR: ${DOXYGEN_SOURCE_ROOT_DIR}")
    message(STATUS "${module}: DOXYGEN_HTML_OUTPUT: ${DOXYGEN_HTML_OUTPUT}")

    # build documentations for losc and c_losc modules with doxygen.
    add_custom_command(OUTPUT ${DOXYGEN_HTML_OUTPUT}
                       DEPENDS
                       ${DOXYGEN_INPUT_DIR}
                       ${DOXYFILE_OUT}
                       losc c_losc py_losc_core
                       ${DOXYGEN_SOURCE_FILES}
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                       MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                       COMMENT "Generating docs with doxygen"
                       )

    # create a target for doxygen docs.
    add_custom_target(docs_doxygen_${module} ALL DEPENDS ${DOXYGEN_HTML_OUTPUT})
endmacro()

# build docs for losc module
build_doxygen("losc")

# build docs for c_losc module
build_doxygen("c_losc")

message(STATUS "LOSC_DOXYGEN_HTML_OUTPUT: ${LOSC_DOXYGEN_HTML_OUTPUT}")
message(STATUS "C_LOSC_DOXYGEN_HTML_OUTPUT: ${C_LOSC_DOXYGEN_HTML_OUTPUT}")

# ==> Sphinx <==
set(SPHINX_BUILD_ROOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_BUILD_DIR ${SPHINX_BUILD_ROOT_DIR}/build)
set(SPHINX_SOURCE_DIR ${SPHINX_BUILD_ROOT_DIR}/source)
set(SPHINX_HTML_OUTPUT ${SPHINX_BUILD_DIR}/index.html)

# copy sphinx source files to building directory
file(MAKE_DIRECTORY ${SPHINX_BUILD_DIR})
file(GLOB SPHINX_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/*)
file(COPY ${SPHINX_SOURCE_FILES} DESTINATION ${SPHINX_SOURCE_DIR})

# configure sphinx source files
configure_file(${SPHINX_SOURCE_DIR}/losc.rst.in ${SPHINX_SOURCE_DIR}/losc.rst @ONLY)
configure_file(${SPHINX_SOURCE_DIR}/c_losc.rst.in ${SPHINX_SOURCE_DIR}/c_losc.rst @ONLY)

# print some informational message.
message(STATUS "SPHINX_SOURCE_DIR: ${SPHINX_SOURCE_DIR}")
message(STATUS "SPHINX_BUILD_DIR: ${SPHINX_BUILD_DIR}")
message(STATUS "SPHINX_HTML_OUTPUT: ${SPHINX_HTML_OUTPUT}")

# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated
add_custom_command(OUTPUT ${SPHINX_HTML_OUTPUT}
                   COMMAND
                   ${SPHINX_EXECUTABLE} -b html
                   # Tell Breathe where to find the Doxygen output
                   #-Dbreathe_projects.losc=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE_DIR} ${SPHINX_BUILD_DIR}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${SPHINX_SOURCE_FILES}
                   losc
                   c_losc
                   py_losc_core
                   MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py
                   COMMENT "Generating documentation with Sphinx")

# create a target for sphinx docs.
add_custom_target(dosc_sphinx ALL DEPENDS ${SPHINX_HTML_OUTPUT})

message(STATUS "Note: Enable to build the documentation.\n"
               "   After the building process, see documentation in:\n"
               "   ${SPHINX_HTML_OUTPUT}")
