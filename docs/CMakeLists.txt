# Use Sphinx + Breathe + Doxygen to generate documentations for C/C++ project.
# Refer to this blog for an example:
# https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/.

# Add the cmake folder, so the FindSphinx module can be found.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# require doxygen and sphinx is installed.
find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# ==> Doxygen <==
set(DOXYGEN_INPUT_DIR ${CMAKE_SOURCE_DIR}/include)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)

set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in)
set(DOXYFILE_OUT ${DOXYGEN_OUTPUT_DIR}/Doxyfile)
set(DOXYGEN_HTML_OUTPUT_DIR ${DOXYGEN_OUTPUT_DIR}/html)
set(DOXYGEN_HTML_OUTPUT ${DOXYGEN_OUTPUT_DIR}/html/index.html)

# create some directories that Doxygen won't create for us
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
file(MAKE_DIRECTORY ${DOXYGEN_HTML_OUTPUT_DIR})

# print some informational message.
message(STATUS "doxygen configuration file: ${DOXYFILE_OUT}")
message(STATUS "DOXYGEN_INPUT_DIR: ${DOXYGEN_INPUT_DIR}")
message(STATUS "DOXYGEN_OUTPUT_DIR: ${DOXYGEN_OUTPUT_DIR}")

# configure the configuration file for doxygen.
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# build documentations for losc and c_losc modules with doxygen.
add_custom_command(OUTPUT ${DOXYGEN_HTML_OUTPUT}
                   DEPENDS ${DOXYGEN_INPUT_DIR} losc c_losc py_losc_core
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs with doxygen"
                   )

# create a target for doxygen docs.
add_custom_target(docs_doxygen ALL DEPENDS ${DOXYGEN_HTML_OUTPUT})

# ==> Sphinx <==
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx/build)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

# print some informational message.
message(STATUS "SPHINX_SOURCE: ${SPHINX_SOURCE}")
message(STATUS "SPHINX_BUILD: ${SPHINX_BUILD}")

file(GLOB SPHINX_SOURCE_ALL_FILES ${SPHINX_SOURCE}/*.rst)
file(MAKE_DIRECTORY ${SPHINX_BUILD}) #Doxygen won't create this for us

# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND
                   ${SPHINX_EXECUTABLE} -b html
                   # Tell Breathe where to find the Doxygen output
                   -Dbreathe_projects.losc=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${SPHINX_SOURCE_ALL_FILES}
                   ${DOXYGEN_INDEX_FILE}
                   losc
                   c_losc
                   py_losc_core
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx")

# create a target for sphinx docs.
add_custom_target(dosc_sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

message(STATUS "Note: Enable to build the documentation.\n"
               "      After the building process, see documentation in:\n"
               "      ${SPHINX_INDEX_FILE}")
