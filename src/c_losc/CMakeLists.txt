set(PROJECT_NAME c_losc)
project(${PROJECT_NAME} C CXX)

# sources
file(GLOB SOURCES "*.cpp")
set(C_LOSC_PUBLIC_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# ==> build c_losc library <==
add_library(
    ${PROJECT_NAME} SHARED
    ${SOURCES})

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    Eigen3::Eigen
    ${PROJECT_LOSC})

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${C_LOSC_PUBLIC_INCLUDE_DIR})

# ==> install losc library <==
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY "${C_LOSC_PUBLIC_INCLUDE_DIR}/c_losc" DESTINATION include)
