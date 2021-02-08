set(PROJECT_NAME ${PROJECT_LOSC})

# build Losc documentation file at configure time.
# first we can indicate the documentation build as an option and set it to ON by
# default
option(BUILD_DOC "Build documentation for Losc library" OFF)
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
   message("Note: If you want to build documentation, use `-DBUILD_DOC=On -DCMAKE_BUILD_TYPE=Release`. "
    "Doxygen will be required.")
endif (BUILD_DOC)

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

# sources
file(GLOB SOURCES "*.cpp")
set(LOSC_PUBLIC_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# ==> build losc library <==
add_library(
    ${PROJECT_NAME} SHARED
    ${SOURCES})

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    Eigen3::Eigen)

if (OPENMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(
            ${PROJECT_NAME}
            PUBLIC
            OpenMP::OpenMP_CXX)
    endif()
endif()

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${LOSC_PUBLIC_INCLUDE_DIR})

# ==> install losc library <==
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY "${LOSC_PUBLIC_INCLUDE_DIR}/losc" DESTINATION include)
