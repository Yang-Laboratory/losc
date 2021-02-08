set(PROJECT_NAME py_losc_core)
project(${PROJECT_NAME} C CXX)

# sources
file(GLOB SOURCES "*.cpp")

pybind11_add_module(
    ${PROJECT_NAME} SHARED
    ${SOURCES})

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_LOSC})

# ==> install py_losc module <==
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
