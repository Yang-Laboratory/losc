set(PROJECT_NAME ${PROJECT_LOSC})

option(OPENMP "Enable openmp" On)
# enable openmp or not.
if (OPENMP)
    message(STATUS "Openmp is enabled.")
    add_compile_definitions(_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "Found openmp for CXX.")
    else()
        message(FATAL_ERROR "Fail to find openmp for CXX. Abort compilation.")
    endif()
else()
    message(STATUS "Openmp is disabled, use `cmake -DOPENMP=on` to turn on.")
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
    target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC
        OpenMP::OpenMP_CXX)
endif()

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${LOSC_PUBLIC_INCLUDE_DIR})

# ==> install losc library <==
install(TARGETS ${PROJECT_NAME} DESTINATION liblosc/losc)
install(DIRECTORY "${LOSC_PUBLIC_INCLUDE_DIR}/losc" DESTINATION liblosc/losc/include)
